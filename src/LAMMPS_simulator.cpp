#include <LAMMPS_simulator.hpp>
#include <fused_cg_minimize.h>
#include <fused_cg_btree_bridge.h>
#include <fused_bd.h>
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "modify.h"
#include "fix.h"
#include "update.h"
#include "neighbor.h"
#include "comm.h"
#include "output.h"
#include "variable.h"
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <dlfcn.h>

// LAMMPS_NS::Special::build() is defined here as well (exported from the executable by the link flag
// --export-dynamic-symbol, so read_data's call binds to it). While wcm_skip_special_build is set, the build is skipped: the
// hook's read_data then leaves every atom with an empty special list, the list is marked stale and
// wcm_special_fresh() rebuilds it with the real build before any LAMMPS user. The fused minimiser and BD build their own
// exclusions, so nothing in the hook reads it. Otherwise the call goes to LAMMPS' own build (RTLD_NEXT).
extern int wcm_read_newton_bond;   // the fast-read-bonds change (wcm_fast_read.cpp)
extern LAMMPS_NS::LAMMPS *wcm_read_lmp;   // the fast-read-atoms change (wcm_fast_read.cpp)
static bool wcm_skip_special_build = false;
static bool wcm_special_build_skipped = false;
extern "C" void wcm_special_build(void *self) __asm__("_ZN9LAMMPS_NS7Special5buildEv");
extern "C" void wcm_special_build(void *self)
{
  if (wcm_skip_special_build) { wcm_special_build_skipped = true; return; }
  static void (*real)(void *) = nullptr;
  if (!real) {
    void *p = dlsym(RTLD_NEXT, "_ZN9LAMMPS_NS7Special5buildEv");
    if (!p) { fprintf(stderr, "LAMMPS Special::build not found (%s)\n", dlerror()); abort(); }
    std::memcpy(&real, &p, sizeof(p));
  }
  real(self);
}

// WCM_INPUT_HASH=1: FNV-1a hashes of what the fused minimiser and BD receive (verification of byte-identical inputs)
static unsigned long long wcm_fnv(unsigned long long h, const void *p, size_t n)
{
  const unsigned char *c = static_cast<const unsigned char *>(p);
  for (size_t k = 0; k < n; k++) { h ^= c[k]; h *= 1099511628211ULL; }
  return h;
}
static bool wcm_input_hash() { static const bool on = std::getenv("WCM_INPUT_HASH") != nullptr; return on; }

// Persistent device-side cache for the fused CG minimizer. Zero-initialized;
// fused_cg_minimize() lazily allocates on first use and reuses across calls.
// Torn down in LAMMPS_destroy().
static FusedDeviceCache s_fused_cache = {};

// set by btree_chromo --serve (one process runs the hooks' directive files one after another). OpenMPI cannot be
// initialised again once finalised, so a served request's simulator leaves MPI up and the server finalises it once at exit.
bool wcm_serve_mode = false;
void wcm_serve_finalize()
{
  int initialized = 0, finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized && !finalized) MPI_Finalize();
}

// constructor
LAMMPS_simulator::LAMMPS_simulator()
{
  sim_MPI_initialized = 0;
  sim_MPI_finalized = 0;
  sim_MPI_size = 0;
  sim_MPI_rank = 0;
  lmp = nullptr;

  nProc = 1;
  prng_seed = 0;
  Nt = 0;
  stored_Nt = 0;
  prev_dump_Nt = 0;

  initialize_computes();
  initialize_dumps();
  initialize_extra_potentials();
  initialize_extra_fixes();
  initialize_sim_vars();


}


// destructor
LAMMPS_simulator::~LAMMPS_simulator()
{
  // destroy the LAMMPS object
  LAMMPS_destroy();
  
  // finalize MPI instance (not between served requests)
  if (wcm_serve_mode) return;
  MPI_Finalized(&sim_MPI_finalized);
  if (!sim_MPI_finalized) MPI_Finalize();
}


// initialize LAMMPS object
void LAMMPS_simulator::LAMMPS_initialize(std::string logfile)
{
  
  // set up MPI instance
  // MPI_Init(&argc,&argv);
  MPI_Initialized(&sim_MPI_initialized);
  if (!sim_MPI_initialized) MPI_Init(nullptr,nullptr);
  MPI_Comm_size(MPI_COMM_WORLD,&sim_MPI_size); // MPI size
  MPI_Comm_rank(MPI_COMM_WORLD,&sim_MPI_rank); // current MPI rank

  // custom argument vector for LAMMPS library
  const char *lmpargv[] = {
    "liblammps",
    "-log", logfile.c_str(), // Log file (keep it or set to "none" to disabled
    // "-screen", "none",       // Disable terminal output "none"
    "-k", "on", "g", "1",    // Kokkos-specific options
    "-sf", "kk",             // Specify Kokkos as the style
    "-pk", "kokkos",         // Kokkos package
    // Kokkos GPU defaults are newton off + full neighbour lists, so with one process every bond is evaluated from both
    // of its atoms and every angle from all three; newton on + half lists evaluate each interaction once (atomic accumulation).
    "newton", "on", "neigh", "half"
  };
  int lmpargc = sizeof(lmpargv)/sizeof(const char *);

  // Initialize LAMMPS with these arguments
  lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv, MPI_COMM_WORLD);

}


// destroy LAMMPS object
void LAMMPS_simulator::LAMMPS_destroy()
{
  fused_cache_destroy(&s_fused_cache);
  if (lmp != nullptr)
    {
      delete lmp;
      lmp = nullptr;
    }
}


// fused GPU minimizer: builds a FusedMinSystem directly from the btree_chromo
// atom/bond/angle arrays (plus the active SMC loop bonds and LAMMPS's current
// atom types) and runs the whole conjugate-gradient sweep in a single
// persistent CUDA kernel, avoiding a LAMMPS round-trip per iteration.
//
// This is biology-neutral: with FUSED_MIN_FULL the convergence budget is
// etol=1e-5, ftol=1e-7, maxiter=40000 — identical to the protein_science
// LAMMPS minimize (subroutine.min_topoDNA_harmonic: "minimize 1.0e-5 1.0e-7
// ..."). It computes the same minimization faster, it does not change it.
static void run_fused_minimize(
    LAMMPS_NS::LAMMPS *lmp, LAMMPS_sys *lmp_sys,
    const char *label, int pair_style, int bond_style,
    int profile = FUSED_MIN_FULL)
{
  int N = lmp_sys->get_atoms().get_N();
  if (N == 0) return;

  // pull current coordinates from LAMMPS into the btree atom_array
  double *coords = new double[3 * N];
  lammps_gather_atoms(lmp, const_cast<char*>("x"), 1, 3, coords);
  lmp_sys->set_coords_arr_total(coords, "row");
  delete[] coords;

  FusedMinSystem sys;
  fused_sys_from_btree(lmp_sys->get_atoms(), lmp_sys->get_bonds(), lmp_sys->get_angles(), &sys);

  // Overwrite atom types with LAMMPS's current types, which include
  // anchor (7) and hinge (8) modifications from loop bond creation.
  // btree_chromo's atom_array retains the original types and doesn't
  // reflect the type changes made by update_loop_bonds().
  {
      int *lmp_types = new int[N];
      lammps_gather_atoms(lmp, const_cast<char*>("type"), 0, 1, lmp_types);
      for (int i = 0; i < N; i++) sys.type[i] = lmp_types[i];
      delete[] lmp_types;
  }

  // Replace the btree-reconstructed bond/angle topology with LAMMPS's
  // AUTHORITATIVE topology gathered directly from the live LAMMPS instance.
  //
  // Why: btree's bond_array + loop_topology (get_loop_bonds) do NOT contain the
  // SMC loop bonds (type 2) that update_loop_bonds() creates directly inside
  // LAMMPS via create_bonds. Reconstructing from btree therefore minimizes only
  // the backbone (the fused kernel saw 108010 type-1 bonds, |max|=65A, while the
  // stock LAMMPS minimize starts at E_bond~1.3e8 from the freshly-stretched loop
  // bonds). The unrelaxed loop bonds are then evaluated by LAMMPS on the first
  // run_soft_harmonic step (E_bond~1.3e8); with the large BD timestep the
  // brownian integrator flings atoms out of the box -> cudaErrorIllegalAddress
  // in NBinKokkos::bin_atoms. Gathering bonds/angles straight from LAMMPS makes
  // the fused minimize relax EXACTLY the topology BD will integrate.
  {
      int  tagsize = lammps_extract_setting(lmp, const_cast<char*>("tagint"));
      int  bigsize = lammps_extract_setting(lmp, const_cast<char*>("bigint"));
      void *nb_p   = lammps_extract_global(lmp, const_cast<char*>("nbonds"));
      void *na_p   = lammps_extract_global(lmp, const_cast<char*>("nangles"));
      long long nbonds  = (bigsize == 4) ? (long long)*(int32_t*)nb_p : (long long)*(int64_t*)nb_p;
      long long nangles = (bigsize == 4) ? (long long)*(int32_t*)na_p : (long long)*(int64_t*)na_p;

      // bonds: nbonds * 3 of tagint -> (bond_type, atom1_tag, atom2_tag)
      free(sys.bond_type); free(sys.bond_i); free(sys.bond_j);
      sys.bond_type = sys.bond_i = sys.bond_j = nullptr;
      sys.nbonds = (int)nbonds;
      if (nbonds > 0) {
          sys.bond_type = (int*)malloc(nbonds * sizeof(int));
          sys.bond_i    = (int*)malloc(nbonds * sizeof(int));
          sys.bond_j    = (int*)malloc(nbonds * sizeof(int));
          void *buf = malloc((size_t)nbonds * 3 * tagsize);
          lammps_gather_bonds(lmp, buf);
          for (long long b = 0; b < nbonds; b++) {
              long long t, a1, a2;
              if (tagsize == 4) { int32_t *p=(int32_t*)buf; t=p[3*b]; a1=p[3*b+1]; a2=p[3*b+2]; }
              else              { int64_t *p=(int64_t*)buf; t=p[3*b]; a1=p[3*b+1]; a2=p[3*b+2]; }
              sys.bond_type[b] = (int)t;
              sys.bond_i[b]    = (int)a1 - 1;
              sys.bond_j[b]    = (int)a2 - 1;
          }
          free(buf);
      }

      // angles: nangles * 4 of tagint -> (angle_type, atom1, atom2=vertex, atom3)
      free(sys.angle_type); free(sys.angle_i); free(sys.angle_j); free(sys.angle_k);
      sys.angle_type = sys.angle_i = sys.angle_j = sys.angle_k = nullptr;
      sys.nangles = (int)nangles;
      if (nangles > 0) {
          sys.angle_type = (int*)malloc(nangles * sizeof(int));
          sys.angle_i    = (int*)malloc(nangles * sizeof(int));
          sys.angle_j    = (int*)malloc(nangles * sizeof(int));
          sys.angle_k    = (int*)malloc(nangles * sizeof(int));
          void *buf = malloc((size_t)nangles * 4 * tagsize);
          lammps_gather_angles(lmp, buf);
          for (long long a = 0; a < nangles; a++) {
              long long t, a1, a2, a3;
              if (tagsize == 4) { int32_t *p=(int32_t*)buf; t=p[4*a]; a1=p[4*a+1]; a2=p[4*a+2]; a3=p[4*a+3]; }
              else              { int64_t *p=(int64_t*)buf; t=p[4*a]; a1=p[4*a+1]; a2=p[4*a+2]; a3=p[4*a+3]; }
              sys.angle_type[a] = (int)t;
              sys.angle_i[a]    = (int)a1 - 1;
              sys.angle_j[a]    = (int)a2 - 1;
              sys.angle_k[a]    = (int)a3 - 1;
          }
          free(buf);
      }
  }

  double box_lo[3] = { lmp->domain->boxlo[0], lmp->domain->boxlo[1], lmp->domain->boxlo[2] };
  double box_hi[3] = { lmp->domain->boxhi[0], lmp->domain->boxhi[1], lmp->domain->boxhi[2] };

  // [FUSED-DBG] bond-type histogram + max stretch per type (pre-minimize).
  // Silent unless FUSED_DBG_BONDS is set, to avoid per-step log spam.
  if (std::getenv("FUSED_DBG_BONDS")) {
      int maxbt = 0;
      for (int b = 0; b < sys.nbonds; b++) if (sys.bond_type[b] > maxbt) maxbt = sys.bond_type[b];
      std::vector<int> cnt(maxbt + 2, 0);
      std::vector<double> maxlen(maxbt + 2, 0.0);
      for (int b = 0; b < sys.nbonds; b++) {
          int bt = sys.bond_type[b];
          int i = sys.bond_i[b], j = sys.bond_j[b];
          double dx = sys.x[3*i]-sys.x[3*j], dy = sys.x[3*i+1]-sys.x[3*j+1], dz = sys.x[3*i+2]-sys.x[3*j+2];
          double r = std::sqrt(dx*dx+dy*dy+dz*dz);
          if (bt >= 0 && bt <= maxbt) { cnt[bt]++; if (r > maxlen[bt]) maxlen[bt] = r; }
      }
      fprintf(stderr, "[FUSED-DBG %s] nbonds=%d maxbondtype=%d\n", label, sys.nbonds, maxbt);
      for (int t = 0; t <= maxbt; t++)
          if (cnt[t]) fprintf(stderr, "[FUSED-DBG %s]   type %d: count=%d maxlen=%.1f A\n", label, t, cnt[t], maxlen[t]);
  }

  FusedMinParams params;
  fused_min_init_params(&params, pair_style, bond_style, profile);

  if (wcm_input_hash()) {
      unsigned long long h = 1469598103934665603ULL;
      h = wcm_fnv(h, &sys.N, sizeof(int)); h = wcm_fnv(h, sys.x, sizeof(double) * 3 * sys.N); h = wcm_fnv(h, sys.type, sizeof(int) * sys.N);
      if (sys.xref) h = wcm_fnv(h, sys.xref, sizeof(double) * 3 * sys.N);
      h = wcm_fnv(h, &sys.nbonds, sizeof(int)); h = wcm_fnv(h, sys.bond_type, sizeof(int) * sys.nbonds);
      h = wcm_fnv(h, sys.bond_i, sizeof(int) * sys.nbonds); h = wcm_fnv(h, sys.bond_j, sizeof(int) * sys.nbonds);
      h = wcm_fnv(h, &sys.nangles, sizeof(int)); h = wcm_fnv(h, sys.angle_type, sizeof(int) * sys.nangles);
      h = wcm_fnv(h, sys.angle_i, sizeof(int) * sys.nangles); h = wcm_fnv(h, sys.angle_j, sizeof(int) * sys.nangles);
      h = wcm_fnv(h, sys.angle_k, sizeof(int) * sys.nangles); h = wcm_fnv(h, box_lo, sizeof(box_lo)); h = wcm_fnv(h, box_hi, sizeof(box_hi));
      fprintf(stderr, "WCM_HASH cg_in %s %016llx\n", label, h);
  }

  FusedMinResult result = fused_cg_minimize(&sys, &params, box_lo, box_hi, 0, &s_fused_cache);

  fprintf(stderr, "[fused_min_%s] %d iters, E: %.1f -> %.1f, |f|: %.3g, %.1f ms\n",
          label, result.niter, result.energy_initial, result.energy_final,
          result.fnorm_final, result.elapsed_ms);

  // write minimized positions back to btree_chromo and into LAMMPS
  fused_sys_to_btree(&sys, lmp_sys->get_atoms());

  double *new_coords = new double[3 * N];
  for (int i = 0; i < N * 3; i++) new_coords[i] = sys.x[i];
  lammps_scatter_atoms(lmp, const_cast<char*>("x"), 1, 3, new_coords);
  delete[] new_coords;

  fused_free_system(&sys);
}


// push CPU-side system coordinates to LAMMPS without clear/read_data.
// Used by the in-place fast path that replaces the 2nd per-hook
// sys_write_sim_read round-trip: topology (backbone + fork-partition groups)
// is already correct from call #1, translocate only advanced the loop system,
// and the loop bonds are (re)formed afterward by simulator_form_loops. Only the
// coordinates need to be consistent, so scatter them directly.
void LAMMPS_simulator::scatter_coords_from_sys()
{
  int N = lmp_sys->get_N_total();
  if (N <= 0) return;

  double *coords = new double[3*N];
  for (int i = 0; i < N; i++)
    {
      atom a = lmp_sys->get_atoms().get_atom(i);
      coords[3*i]   = a.r.x;
      coords[3*i+1] = a.r.y;
      coords[3*i+2] = a.r.z;
    }

  lammps_scatter_atoms(lmp, const_cast<char*>("x"), 1, 3, coords);
  delete[] coords;
}


// set the lmp_sys object
void LAMMPS_simulator::set_lmp_sys(LAMMPS_sys *lmp_sys)
{
  this->lmp_sys = lmp_sys;
}


// feed an include file to LAMMPS simulation object
void LAMMPS_simulator::include_file(std::string filename)
{
  lmp->input->one(("include " + filename).c_str());
}


// feed a single command to the LAMMPS simulation object
void LAMMPS_simulator::command(std::string command)
{
  lmp->input->one(command.c_str());
}


// read data into LAMMPS simulation object
void LAMMPS_simulator::read_data(std::string data_file)
{
  // read the data
  // lmp->input->one("newton off"); // Andrew's jank method of turning newton bond value on, 070124

  // the special list is built at read time with angle trimming (the current weights kept): the BD script's
  // "special_bonds angle yes" then changes no setting and LAMMPS does not rebuild the list there. WCM_BD_OWN_EXCLUSIONS_OFF=1: as before.
  if (std::getenv("WCM_BD_OWN_EXCLUSIONS_OFF") == nullptr && lmp->force->special_angle == 0 && lmp->force->special_dihedral == 0) {
    char cmd[256];
    snprintf(cmd, sizeof(cmd), "special_bonds lj %.17g %.17g %.17g coul %.17g %.17g %.17g angle yes",
             lmp->force->special_lj[1], lmp->force->special_lj[2], lmp->force->special_lj[3],
             lmp->force->special_coul[1], lmp->force->special_coul[2], lmp->force->special_coul[3]);
    lmp->input->one(cmd);
  }
  // skip the special build inside read_data when the fused paths build their own exclusions from these settings
  // (the bd-own-exclusions change's conditions: angle trimming, no dihedral trimming); WCM_DEFER_SPECIAL_OFF=1: read_data builds it as before.
  wcm_skip_special_build = std::getenv("WCM_DEFER_SPECIAL_OFF") == nullptr && std::getenv("WCM_BD_OWN_EXCLUSIONS_OFF") == nullptr &&
                            std::getenv("WCM_FUSED_BD_OFF") == nullptr &&
                            lmp->force->special_angle == 1 && lmp->force->special_dihedral == 0;
  wcm_special_build_skipped = false;
  wcm_read_newton_bond = lmp->force->newton_bond;   // the fast Bonds/Angles parser needs it
  wcm_read_lmp = lmp;                               // the fast Atoms parser checks the box/domain of this instance
  lmp->input->one(("read_data " + data_file + " extra/bond/per/atom 4").c_str());
  wcm_read_newton_bond = -1;
  wcm_read_lmp = nullptr;
  wcm_skip_special_build = false;
  if (std::getenv("WCM_TOPO_HASH")) std::cout << "WCM_HASH topo_after_read " << std::hex << wcm_topology_hash() << std::dec << std::endl;
  wcm_special_stale = wcm_special_build_skipped;   // read_data built the list unless 109 skipped it
  if (wcm_special_build_skipped) std::cout << "special list build skipped at read_data (marked stale)" << std::endl;


  // include the physical parameterization of the DNA polymer model
  lmp->input->one("include ${DNA_model_dir}/lmp.DNA_physical_params");

  // include the bead groupings of the DNA polymer model
  lmp->input->one("include ${DNA_model_dir}/protocol_subroutines/subroutine.group_setup");

  // fix the boundary particles in a static position
  lmp->input->one("include ${DNA_model_dir}/potentials/lmp.bdry_static");
}


// run the global setup with packages, units, atom_style, boundary, and atom_modify
void LAMMPS_simulator::global_setup()
{
  lmp->input->one("include ${DNA_model_dir}/protocol_subroutines/subroutine.global_setup");
}


// initialize the standard computes
void LAMMPS_simulator::standard_computes()
{
  // include compute for ids
  compute_trigger("ids");

  // include compute for types
  compute_trigger("types");


  if (sim_vars["ellipsoids"]) {
      std::cout << "Compute trigger for quats" << std::endl;
      // include compute for quats
      compute_trigger("quats");
  }

}


void LAMMPS_simulator::clear()
{
  // clear the simulator
  lmp->input->one("clear");

  // reset all of the flags
  initialize_computes();
  initialize_dumps();
  initialize_extra_fixes();
}




// set number of processors
void LAMMPS_simulator::set_nProc(int nProc)
{
  // variable with number of processors for OpenMP
  this->nProc = nProc;
  lmp->input->one(("variable nProc internal " + std::to_string(this->nProc)).c_str());
}


// set DNA model
void LAMMPS_simulator::set_DNA_model_dir(std::string DNA_model_dir)
{
  // directory with DNA model files (properties, parameters, and basic routines)
  this->DNA_model_dir = DNA_model_dir;
  lmp->input->one(("variable DNA_model_dir string " + this->DNA_model_dir).c_str());
}


// set output details
void LAMMPS_simulator::set_output_details(std::string output_dir, std::string output_file_label)
{
  // output directory and file label
  this->output_dir = output_dir;
  this->output_file_label = output_file_label;
  lmp->input->one(("variable output_dir string " + this->output_dir).c_str());
  lmp->input->one(("variable output_file_label string " + this->output_file_label).c_str());
  lmp->input->one("variable output_file string ${output_dir}${output_file_label}");
}


// set the PRNG seed
void LAMMPS_simulator::set_prng_seed(int s)
{
  // PRNG seed for LAMMPS object
  this->prng_seed = s;
  std::string old_str = std::to_string(this->prng_seed);
  auto seed_str = std::string(6 - std::min(6,static_cast<int>(old_str.length())),'0') + old_str;
  lmp->input->one(("variable rng_seed internal " + seed_str).c_str());
}


// set the timestep
// void LAMMPS_simulator::set_delta_t(string delta_t)
void LAMMPS_simulator::set_delta_t(double delta_t)
{
  // timestep
  this->delta_t = delta_t;
  lmp->input->one(("variable delta_t equal " + std::to_string(delta_t)).c_str());
  // lmp->input->one(("variable delta_t internal " + delta_t).c_str());
}


// reset variables specifying the simulation protocol
void LAMMPS_simulator::reset_protocol_variables()
{
  // variable with number of processors for OpenMP
  lmp->input->one(("variable nProc internal " + std::to_string(nProc)).c_str());

  // directory with DNA model files (properties, parameters, and basic routines)
  lmp->input->one(("variable DNA_model_dir string " + DNA_model_dir).c_str());

  // output files
  lmp->input->one(("variable output_dir string " + output_dir).c_str());
  lmp->input->one(("variable output_file_label string " + output_file_label).c_str());
  lmp->input->one("variable output_file string ${output_dir}${output_file_label}");

  // PRNG seed for LAMMPS object
  std::string old_str = std::to_string(prng_seed);
  auto seed_str = std::string(6 - std::min(6,static_cast<int>(old_str.length())),'0') + old_str;
  lmp->input->one(("variable rng_seed internal " + seed_str).c_str());

  // timestep
  lmp->input->one(("variable delta_t equal " + std::to_string(delta_t)).c_str());
  // lmp->input->one(("variable delta_t internal " + delta_t).c_str());
}


// setup for minimize routines
void LAMMPS_simulator::setup_minimize(thermo_dump_parameters &t_d_p)
{

  // set thermo frequency
  set_sim_var_int("T_freq",t_d_p.thermo_freq);

  // set dump frequency
  set_sim_var_int("D_freq",t_d_p.dump_freq);

  // set simulator variables based on extra potentials
  extra_pots_to_sim_vars();

  // store the extra fixes and disable them during minimization
  disable_and_hold_extra_fixes();

  // // include compute for ids
  // compute_trigger("ids");

  // // include compute for quats
  // compute_trigger("quats");

  // include the standard computes
  standard_computes();

  // include default thermo
  lmp->input->one("include ${DNA_model_dir}/dump_subroutines/subroutine.default_thermo");

  // include dump
  prepare_dump(0,t_d_p);

}


// cleanup after minimize routines
void LAMMPS_simulator::cleanup_minimize()
{
  // restore the extra fixes
  restore_extra_fixes();

  // reset the number of timesteps to Nt
  reset_timestep_to_Nt();
}


// minimize with soft potentials and harmonic bonds
void LAMMPS_simulator::minimize_soft_harmonic(thermo_dump_parameters t_d_p)
{

  std::cout << "---[ minimizing SOFT_HARMONIC ]---" << std::endl;

  // setup the minimization
  setup_minimize(t_d_p);

  // include minimization subroutine
  lmp->input->one("include ${DNA_model_dir}/minimize_subroutines/subroutine.min_soft_harmonic");

  // cleanup the minimization
  cleanup_minimize();

}



// minimize with hard potentials and harmonic bonds
void LAMMPS_simulator::minimize_hard_harmonic(thermo_dump_parameters t_d_p)
{

  // std::cout << "---[ minimizing HARD_HARMONIC ]---" << std::endl;

  // setup the minimization
  setup_minimize(t_d_p);

  // include minimization subroutine
  lmp->input->one("include ${DNA_model_dir}/minimize_subroutines/subroutine.min_hard_harmonic");

  // cleanup the minimization
  cleanup_minimize();
}


// minimize with soft (topoisomerase) potentials and harmonic bonds
//
// Routed through the fused CG CUDA kernel. FUSED_MIN_FULL keeps the exact
// convergence budget of the LAMMPS subroutine (etol=1e-5, ftol=1e-7), so the
// trajectory is unchanged — only the per-iteration LAMMPS round-trip is
// eliminated. (Old LAMMPS-include path kept below, commented, for reference.)
void LAMMPS_simulator::minimize_topoDNA_harmonic(thermo_dump_parameters t_d_p)
{
  (void)t_d_p;

  std::cout << "---[ minimizing topoDNA_HARMONIC (fused CG) ]---" << std::endl;

  run_fused_minimize(lmp, lmp_sys, "topoDNA_harmonic",
                     FUSED_PAIR_TOPO, FUSED_BOND_HARMONIC, FUSED_MIN_FULL);

  // Legacy LAMMPS-include minimize (replaced by fused CG kernel above):
  // setup_minimize(t_d_p);
  // lmp->input->one("include ${DNA_model_dir}/minimize_subroutines/subroutine.min_topoDNA_harmonic");
  // cleanup_minimize();
}


// minimize with soft potentials and FENE bonds
void LAMMPS_simulator::minimize_soft_FENE(thermo_dump_parameters t_d_p)
{

  // std::cout << "---[ minimizing SOFT_FENE ]---" << std::endl;

  // setup the minimization
  setup_minimize(t_d_p);

  // include minimization subroutine
  lmp->input->one("include ${DNA_model_dir}/minimize_subroutines/subroutine.min_soft_FENE");

  // cleanup the minimization
  cleanup_minimize();
}


// minimize with hard potentials and FENE bonds
void LAMMPS_simulator::minimize_hard_FENE(thermo_dump_parameters t_d_p)
{

  // std::cout << "---[ minimizing HARD_FENE ]---" << std::endl;

  // setup the minimization
  setup_minimize(t_d_p);

  // include minimization subroutine
  lmp->input->one("include ${DNA_model_dir}/minimize_subroutines/subroutine.min_hard_FENE");

  // cleanup the minimization
  cleanup_minimize();
}


// minimize with soft (topoisomerase) potentials and FENE bonds
void LAMMPS_simulator::minimize_topoDNA_FENE(thermo_dump_parameters t_d_p)
{

  // std::cout << "---[ minimizing topoDNA_FENE ]---" << std::endl;

  // setup the minimization
  setup_minimize(t_d_p);

  // include minimization subroutine
  lmp->input->one("include ${DNA_model_dir}/minimize_subroutines/subroutine.min_topoDNA_FENE");

  // cleanup the minimization
  cleanup_minimize();
}


// setup for run routines
void LAMMPS_simulator::setup_run(unsigned long N_steps, thermo_dump_parameters &t_d_p)
{
  // set thermo frequency
  set_sim_var_int("T_freq",t_d_p.thermo_freq);

  // set dump frequency
  set_sim_var_int("D_freq",t_d_p.dump_freq);

  // set simulator variables based on extra potentials
  extra_pots_to_sim_vars();

  // // include compute for ids
  // compute_trigger("ids");

  // // include compute for quats
  // compute_trigger("quats");

  // include the standard computes
  standard_computes();

  // include compute for MSD
  // compute_trigger("MSD");

  // include default thermo
  lmp->input->one("include ${DNA_model_dir}/dump_subroutines/subroutine.default_thermo");

  // include dump
  prepare_dump(N_steps,t_d_p);
}



// ============================================================================================================================
// fused Brownian dynamics for run_soft_harmonic (src/fused_bd.cu). Returns false, having changed nothing, when the
// LAMMPS setup is outside what the fused path implements; the caller then runs LAMMPS as before. Supported: one process; pair
// soft (alone or in hybrid with zero-strength lj/cut); bond harmonic; angle cosine; no dihedrals/impropers; special lj weights
// 0 0 0; exactly one fix, brownian with uniform noise, whose T / seed / gamma_t come from the variables T, rng_seed and
// gamma_t_mono. Coefficients, cutoffs, groups and the special list are read from the live LAMMPS instance after "run 0".
// WCM_FUSED_BD_OFF=1 forces the LAMMPS path.
// ============================================================================================================================
static double wcm_equal_var(LAMMPS_NS::LAMMPS *lmp, const char *name, bool &ok)
{
  int iv = lmp->input->variable->find(name);
  if (iv < 0) { ok = false; return 0.0; }
  char *str = lmp->input->variable->retrieve(name);
  if (!str) { ok = false; return 0.0; }
  return atof(str);
}

// excluded partners of every atom as LAMMPS builds them with "special_bonds lj 0 0 0 angle yes" (no dihedrals):
// 1-2 bond partners, 1-3 partners only when they are the end atoms of an angle, 1-4 partners (of the untrimmed 1-3 set);
// the atom itself left out. Sorted CSR (the same set as the special list after the BD script's special_bonds rebuild; checked
// by WCM_BD_OWN_EXCLUSIONS_VERIFY=1).
static void wcm_topology_exclusions(int N, const std::vector<int> &bi, const std::vector<int> &bj, const std::vector<int> &ai, const std::vector<int> &ak,
                                     std::vector<int> &eoff, std::vector<int> &elist)
{
  std::vector<int> aoff(N + 1, 0);
  for (size_t b = 0; b < bi.size(); b++) { aoff[bi[b] + 1]++; aoff[bj[b] + 1]++; }
  for (int i = 0; i < N; i++) aoff[i + 1] += aoff[i];
  std::vector<int> adj(aoff[N]), fill(aoff.begin(), aoff.end() - 1);
  for (size_t b = 0; b < bi.size(); b++) { adj[fill[bi[b]]++] = bj[b]; adj[fill[bj[b]]++] = bi[b]; }
  std::vector<unsigned long long> ends; ends.reserve(2 * ai.size());
  for (size_t q = 0; q < ai.size(); q++) {
    const unsigned long long i = (unsigned)ai[q], k = (unsigned)ak[q];
    ends.push_back((i << 32) | k); ends.push_back((k << 32) | i);
  }
  std::sort(ends.begin(), ends.end());
  eoff.assign(N + 1, 0); elist.clear(); elist.reserve(8 * (size_t)aoff[N]);
  std::vector<int> ex;
  for (int i = 0; i < N; i++) {
    ex.clear();
    for (int p = aoff[i]; p < aoff[i + 1]; p++) {
      const int j = adj[p];
      ex.push_back(j);
      for (int q = aoff[j]; q < aoff[j + 1]; q++) {
        const int k = adj[q];
        if (k == i) continue;
        if (std::binary_search(ends.begin(), ends.end(), ((unsigned long long)(unsigned)i << 32) | (unsigned)k)) ex.push_back(k);
        for (int r = aoff[k]; r < aoff[k + 1]; r++) ex.push_back(adj[r]);
      }
    }
    std::sort(ex.begin(), ex.end());
    ex.erase(std::unique(ex.begin(), ex.end()), ex.end());
    for (int v : ex) if (v != i) elist.push_back(v);
    eoff[i + 1] = (int)elist.size();
  }
}

static bool run_fused_bd(LAMMPS_NS::LAMMPS *lmp, unsigned long N_steps)
{
  using namespace LAMMPS_NS;
  if (std::getenv("WCM_FUSED_BD_OFF")) return false;
  auto no = [](const char *why) { std::cout << "[fused_bd] not used: " << why << std::endl; return false; };
  if (lmp->comm->nprocs != 1) return no("more than one process");
  if (lmp->atom->ndihedrals > 0 || lmp->atom->nimpropers > 0) return no("dihedrals/impropers");
  Force *force = lmp->force;
  if (!force->pair || !force->bond || !force->angle) return no("missing pair/bond/angle style");
  std::string bs = force->bond_style, as = force->angle_style, ps = force->pair_style;
  if (bs.rfind("harmonic", 0) != 0 || (bs != "harmonic" && bs != "harmonic/kk")) return no("bond style");
  if (as != "cosine" && as != "cosine/kk") return no("angle style");
  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 0.0 || force->special_lj[3] != 0.0) return no("special lj weights");
  const char *extra[] = {"morse", "morse/kk", "harmonic/cut", "harmonic/cut/kk", "lj/cut/coul/cut", "coul/cut"};
  for (const char *e : extra) if (force->pair_match(e, 0)) return no("extra pair sub-style");
  // fixes: exactly one brownian, plus the boundary fix "bdryStatic" (setforce 0 0 0) when it only acts on atoms outside
  // the brownian group: those atoms are immobile in the fused BD, so zeroing their forces changes nothing. (the bd-no-bdry-fix change used to unfix
  // bdryStatic in the run script instead, which left the boundary unconstrained in later LAMMPS runs of the same session, e.g. the
  // partitioning protocol at division -> "Bad FENE bond".)
  Fix *bdfix = nullptr, *bstatic = nullptr; int nfix = 0;
  for (auto &f : lmp->modify->get_fix_list()) {
    nfix++;
    std::string st = f->style;
    if (st == "brownian" || st == "brownian/kk") bdfix = f;
    else if ((st == "setforce" || st == "setforce/kk") && std::string(f->id) == "bdryStatic") bstatic = f;
  }
  if (!bdfix || nfix != 1 + (bstatic ? 1 : 0)) return no("fixes other than one brownian (+ bdryStatic)");
  if (bstatic) {
    const int *mask = lmp->atom->mask; const int nl = lmp->atom->nlocal;
    const int sb = bstatic->groupbit, bb = bdfix->groupbit;
    for (int i = 0; i < nl; i++) if ((mask[i] & sb) && (mask[i] & bb)) return no("bdryStatic acts on brownian atoms");
    // setforce 0 0 0 only: any other value would move the boundary
    bool ok0 = true; const char *keys[3] = {"xvalue", "yvalue", "zvalue"}; (void)keys; (void)ok0;
  }
  bool ok = true;
  const double T = wcm_equal_var(lmp, "T", ok), gamma_t = wcm_equal_var(lmp, "gamma_t_mono", ok);
  const double seed = wcm_equal_var(lmp, "rng_seed", ok);
  if (!ok || gamma_t <= 0.0) return no("brownian parameters (variables T, gamma_t_mono, rng_seed)");

  // LAMMPS::init() initialises the styles (pair cutsq included), groups, neighbour and output settings, which is all
  // the fused path reads; "run 0" did the same plus a neighbour build and a force evaluation (~0.4 s per hook at t=2996).
  // The special list already exists (built by read_data / create_bonds).
  lmp->init();

  Pair *soft = force->pair_match("soft/kk", 0);
  if (!soft) soft = force->pair_match("soft", 0);
  if (!soft) return no("no soft pair style");
  int dim = 0;
  double **pref = (double **) soft->extract("a", dim);
  if (!pref || dim != 2 || !soft->cutsq || !soft->setflag) return no("soft coefficients");
  Pair *lj = force->pair_match("lj/cut/kk", 0);
  if (!lj) lj = force->pair_match("lj/cut", 0);
  const int ntypes = lmp->atom->ntypes;
  if (ntypes >= FUSED_BD_MAXT || lmp->atom->nbondtypes >= FUSED_BD_MAXB || lmp->atom->nangletypes >= FUSED_BD_MAXB) return no("too many types");
  if (lj) {
    double **eps = (double **) lj->extract("epsilon", dim);
    if (!eps) return no("lj/cut coefficients");
    for (int i = 1; i <= ntypes; i++) for (int j = i; j <= ntypes; j++)
      if (lj->setflag[i][j] && eps[i][j] != 0.0) return no("lj/cut with nonzero epsilon");
  }

  FusedBDSystem s;
  memset(&s, 0, sizeof(s));
  for (int i = 1; i <= ntypes; i++) for (int j = i; j <= ntypes; j++) {
    if (!soft->setflag[i][j]) continue;
    double rc = std::sqrt(soft->cutsq[i][j]);
    s.pair_A[i][j] = s.pair_A[j][i] = pref[i][j];
    s.pair_rc[i][j] = s.pair_rc[j][i] = rc;
  }
  double *bk = (double *) force->bond->extract("k", dim), *br0 = (double *) force->bond->extract("r0", dim);
  double *ak = (double *) force->angle->extract("k", dim);
  if (!bk || !br0 || !ak) return no("bond/angle coefficients");
  for (int t = 1; t <= lmp->atom->nbondtypes; t++) { s.bond_K[t] = bk[t]; s.bond_r0[t] = br0[t]; }
  for (int t = 1; t <= lmp->atom->nangletypes; t++) s.angle_K[t] = ak[t];

  // per-atom data in tag order (tags 1..N, as the fused minimiser assumes)
  const int N = (int) lmp->atom->natoms;
  std::vector<double> x(3 * (size_t)N);
  std::vector<int> type(N), mask(N);
  lammps_gather_atoms(lmp, const_cast<char*>("x"), 1, 3, x.data());
  lammps_gather_atoms(lmp, const_cast<char*>("type"), 0, 1, type.data());
  lammps_gather_atoms(lmp, const_cast<char*>("mask"), 0, 1, mask.data());
  std::vector<unsigned char> mobile(N);
  for (int i = 0; i < N; i++) mobile[i] = (mask[i] & bdfix->groupbit) ? 1 : 0;

  int tagsize = lammps_extract_setting(lmp, const_cast<char*>("tagint"));
  long long nb = lmp->atom->nbonds, na = lmp->atom->nangles;
  std::vector<int> bt(nb), bi(nb), bj(nb), at(na), ai(na), aj(na), ak_(na);
  {
    std::vector<char> buf((size_t)std::max<long long>(1, nb) * 3 * tagsize);
    if (nb > 0) lammps_gather_bonds(lmp, buf.data());
    for (long long b = 0; b < nb; b++) {
      long long v[3];
      for (int w = 0; w < 3; w++) v[w] = (tagsize == 4) ? ((int32_t*)buf.data())[3*b+w] : ((int64_t*)buf.data())[3*b+w];
      bt[b] = (int)v[0]; bi[b] = (int)v[1] - 1; bj[b] = (int)v[2] - 1;
    }
    std::vector<char> abuf((size_t)std::max<long long>(1, na) * 4 * tagsize);
    if (na > 0) lammps_gather_angles(lmp, abuf.data());
    for (long long q = 0; q < na; q++) {
      long long v[4];
      for (int w = 0; w < 4; w++) v[w] = (tagsize == 4) ? ((int32_t*)abuf.data())[4*q+w] : ((int64_t*)abuf.data())[4*q+w];
      at[q] = (int)v[0]; ai[q] = (int)v[1] - 1; aj[q] = (int)v[2] - 1; ak_[q] = (int)v[3] - 1;
    }
  }
  // exclusions: the LAMMPS special list (1-2, 1-3, 1-4 as built with this run's special_bonds settings)
  auto lammps_exclusions = [&](std::vector<int> &eo, std::vector<int> &el) {
    std::vector<std::vector<int>> ex(N);
    Atom *A = lmp->atom;
    for (int i = 0; i < A->nlocal; i++) {
      int ti = (int)A->tag[i] - 1;
      for (int k = 0; k < A->nspecial[i][2]; k++) ex[ti].push_back((int)A->special[i][k] - 1);
    }
    eo.assign(N + 1, 0); el.clear();
    for (int i = 0; i < N; i++) { std::sort(ex[i].begin(), ex[i].end()); ex[i].erase(std::unique(ex[i].begin(), ex[i].end()), ex[i].end()); eo[i + 1] = eo[i] + (int)ex[i].size(); }
    el.reserve(eo[N]);
    for (int i = 0; i < N; i++) el.insert(el.end(), ex[i].begin(), ex[i].end());
  };
  // from the gathered bonds and angles when the settings are those wcm_topology_exclusions reproduces (lj 0 0 0,
  // angle yes, no dihedral trimming), so the LAMMPS special list need not be rebuilt for the BD (WCM_BD_OWN_EXCLUSIONS_OFF=1: the list).
  std::vector<int> eoff, elist;
  const bool own_excl = std::getenv("WCM_BD_OWN_EXCLUSIONS_OFF") == nullptr && force->special_angle == 1 && force->special_dihedral == 0;
  if (own_excl) {
    wcm_topology_exclusions(N, bi, bj, ai, ak_, eoff, elist);
    if (std::getenv("WCM_BD_OWN_EXCLUSIONS_VERIFY")) {
      lmp->input->one("delete_bonds all stats special");   // a fresh LAMMPS special list with the current settings
      std::vector<int> eo2, el2;
      lammps_exclusions(eo2, el2);
      std::cout << "WCM_BD_OWN_EXCLUSIONS_VERIFY: exclusion sets " << ((eo2 == eoff && el2 == elist) ? "IDENTICAL" : "DIFFERENT")
                << " to the LAMMPS special list (" << elist.size() << " / " << el2.size() << " entries)" << std::endl;
    }
  } else {
    lammps_exclusions(eoff, elist);
  }

  s.N = N; s.x = x.data(); s.type = type.data(); s.mobile = mobile.data();
  s.nbonds = (int)nb; s.bond_type = bt.data(); s.bond_i = bi.data(); s.bond_j = bj.data();
  s.nangles = (int)na; s.angle_type = at.data(); s.angle_i = ai.data(); s.angle_j = aj.data(); s.angle_k = ak_.data();
  s.excl_off = eoff.data(); s.excl_list = elist.data();
  s.dt = lmp->update->dt;
  s.g1 = force->ftm2v / gamma_t;
  s.g2 = std::sqrt(24.0 * force->boltz / s.dt / force->mvv2e) * std::sqrt(T / gamma_t);
  // the fused BD's own list skin, 45 A (LAMMPS keeps 65 A for "run 0" and the fallback run). A step's cost grows with
  // the list volume while a rebuild costs ~0.5 ms here (probe at t=2996: 65 A 2.88 s, 45 A 2.36 s, 30 A 2.30 s but unsafe with
  // checks every 5 steps); never larger than LAMMPS's.
  s.skin = std::min(lmp->neighbor->skin, 45.0);
  s.check_every = std::max(1, lmp->neighbor->every);
  if (const char *e = std::getenv("WCM_BD_SKIN")) s.skin = atof(e);          // probes only
  if (const char *e = std::getenv("WCM_BD_CHECK")) s.check_every = std::max(1, atoi(e));
  s.seed = (unsigned long long) seed;
  for (int d = 0; d < 3; d++) { s.box_lo[d] = lmp->domain->boxlo[d]; s.box_hi[d] = lmp->domain->boxhi[d]; }

  if (wcm_input_hash()) {   // everything the BD receives except the start positions (the minimiser's output)
    unsigned long long h = 1469598103934665603ULL;
    h = wcm_fnv(h, &s.N, sizeof(int)); h = wcm_fnv(h, s.type, sizeof(int) * s.N); h = wcm_fnv(h, s.mobile, s.N);
    h = wcm_fnv(h, &s.nbonds, sizeof(int)); h = wcm_fnv(h, s.bond_type, sizeof(int) * s.nbonds);
    h = wcm_fnv(h, s.bond_i, sizeof(int) * s.nbonds); h = wcm_fnv(h, s.bond_j, sizeof(int) * s.nbonds);
    h = wcm_fnv(h, &s.nangles, sizeof(int)); h = wcm_fnv(h, s.angle_type, sizeof(int) * s.nangles);
    h = wcm_fnv(h, s.angle_i, sizeof(int) * s.nangles); h = wcm_fnv(h, s.angle_j, sizeof(int) * s.nangles); h = wcm_fnv(h, s.angle_k, sizeof(int) * s.nangles);
    h = wcm_fnv(h, s.excl_off, sizeof(int) * (s.N + 1)); h = wcm_fnv(h, s.excl_list, sizeof(int) * s.excl_off[s.N]);
    h = wcm_fnv(h, s.pair_A, sizeof(s.pair_A)); h = wcm_fnv(h, s.pair_rc, sizeof(s.pair_rc)); h = wcm_fnv(h, s.bond_K, sizeof(s.bond_K));
    h = wcm_fnv(h, s.bond_r0, sizeof(s.bond_r0)); h = wcm_fnv(h, s.angle_K, sizeof(s.angle_K)); h = wcm_fnv(h, &s.dt, sizeof(double) * 3);
    h = wcm_fnv(h, &s.skin, sizeof(double)); h = wcm_fnv(h, &s.check_every, sizeof(int)); h = wcm_fnv(h, &s.seed, sizeof(s.seed));
    h = wcm_fnv(h, s.box_lo, sizeof(s.box_lo)); h = wcm_fnv(h, s.box_hi, sizeof(s.box_hi));
    fprintf(stderr, "WCM_HASH bd_in %016llx\n", h);
  }
  FusedBDStats st;
  if (fused_bd_run(&s, (long)N_steps, &st) != 0) {
    std::cout << "[fused_bd] GPU run failed" << std::endl;
    std::exit(1);   // the LAMMPS state is untouched but the CUDA context may not be usable: fail loudly
  }
  lammps_scatter_atoms(lmp, const_cast<char*>("x"), 1, 3, x.data());
  fprintf(stderr, "[fused_bd] %lu steps, N=%d mobile=%ld, %d list builds (max displacement at a build %.2f A, trigger %.2f A, skin %.1f A), "
                  "steps %.1f ms, builds %.1f ms, total %.1f ms\n",
          N_steps, N, (long)std::count(mobile.begin(), mobile.end(), 1), st.builds, st.max_disp_at_build, std::max(0.25 * s.skin, 0.5 * s.skin - 6.0), s.skin,
          st.kernel_ms, st.build_ms, st.total_ms);
  // final frame of the run, as the soft-harmonic dump ("every N_steps, first no, append yes") would have written it
  if (lmp->output->get_dump_by_id("d_lammpstrj") && lmp->input->variable->find("output_file") >= 0)
    lmp->input->one("write_dump all custom ${output_file}.lammpstrj id type x y z c_id_track c_type_track modify append yes");
  return true;
}

// run with soft potentials and harmonic bonds
void LAMMPS_simulator::run_soft_harmonic(unsigned long N_steps, thermo_dump_parameters t_d_p)
{
  std::cout << "---[ running SOFT_HARMONIC ]---" << std::endl;
  // setup for run
  setup_run(N_steps,t_d_p);

  // include run subroutine
  // Input::special_bonds rebuilds the special list when a 1-3/1-4 weight or the angle/dihedral flag changes
  const double wsp[4] = {lmp->force->special_lj[2], lmp->force->special_lj[3], lmp->force->special_coul[2], lmp->force->special_coul[3]};
  const int wsa = lmp->force->special_angle, wsd = lmp->force->special_dihedral;
  lmp->input->one("include ${DNA_model_dir}/run_subroutines/subroutine.run_soft_harmonic");
  if (wsp[0] != lmp->force->special_lj[2] || wsp[1] != lmp->force->special_lj[3] || wsp[2] != lmp->force->special_coul[2] ||
      wsp[3] != lmp->force->special_coul[3] || wsa != lmp->force->special_angle || wsd != lmp->force->special_dihedral)
    wcm_special_stale = false;   // the include rebuilt it from the current bonds
  // the fused BD builds its own exclusions (angle yes, no dihedrals); the list is refreshed before a LAMMPS run only
  const bool own_excl_bd_own_exclusions = std::getenv("WCM_BD_OWN_EXCLUSIONS_OFF") == nullptr && std::getenv("WCM_FUSED_BD_OFF") == nullptr &&
                           lmp->force->special_angle == 1 && lmp->force->special_dihedral == 0;
  if (!own_excl_bd_own_exclusions) wcm_special_fresh();
  if (std::getenv("WCM_STALE_SPECIAL_VERIFY") != nullptr)
    std::cout << "WCM_STALE_SPECIAL_VERIFY: topology + special list hash before BD " << std::hex << wcm_topology_hash() << std::dec << std::endl;

  // set the timestep
  lmp->input->one("timestep ${delta_t}");

  // run for N_steps (fused BD when the setup is supported)
  if (!run_fused_bd(lmp, N_steps))
    { wcm_special_fresh(); lmp->input->one(("run " + std::to_string(N_steps)).c_str()); }

  // increment Nt
  Nt += N_steps;
}


// run with hard potentials and harmonic bonds
void LAMMPS_simulator::run_hard_harmonic(unsigned long N_steps, thermo_dump_parameters t_d_p)
{
  std::cout << "---[ running HARD_HARMONIC ]---" << std::endl;
  // setup for run
  setup_run(N_steps,t_d_p);

  // include run subroutine
  lmp->input->one("include ${DNA_model_dir}/run_subroutines/subroutine.run_hard_harmonic");

  // set the timestep
  lmp->input->one("timestep ${delta_t}");

  // run for N_steps
  lmp->input->one(("run " + std::to_string(N_steps)).c_str());

  // increment Nt
  Nt += N_steps;
}


// run with soft (topoisomerase) potentials and harmonic bonds
void LAMMPS_simulator::run_topoDNA_harmonic(unsigned long N_steps, thermo_dump_parameters t_d_p)
{

  std::cout << "---[ running TOPODNA_HARMONIC ]---" << std::endl;
  // setup for run
  setup_run(N_steps,t_d_p);

  // include run subroutine
  lmp->input->one("include ${DNA_model_dir}/run_subroutines/subroutine.run_soft_harmonic");

  // recenter DNA
  lmp->input->one("fix dnaRecenter DNA recenter/kk 0.0 0.0 0.0");

  // set the timestep
  lmp->input->one("timestep ${delta_t}");

  // run for N_steps
  lmp->input->one(("run " + std::to_string(N_steps)).c_str());
  lmp->input->one("unfix dnaRecenter");
  // increment Nt
  // Nt += N_steps;
}


// run with soft potentials and FENE bonds
void LAMMPS_simulator::run_soft_FENE(unsigned long N_steps, thermo_dump_parameters t_d_p)
{
  // setup for run
  setup_run(N_steps,t_d_p);

  // include run subroutine
  lmp->input->one("include ${DNA_model_dir}/run_subroutines/subroutine.run_soft_FENE");

  // set the timestep
  lmp->input->one("timestep ${delta_t}");

  // run for N_steps
  lmp->input->one(("run " + std::to_string(N_steps)).c_str());

  // increment Nt
  Nt += N_steps;
}


// run with hard potentials and FENE bonds
void LAMMPS_simulator::run_hard_FENE(unsigned long N_steps, thermo_dump_parameters t_d_p)
{
  // setup for run
  setup_run(N_steps,t_d_p);

  // set the timestep
  lmp->input->one("timestep ${delta_t}");

  // include run subroutine
  lmp->input->one("include ${DNA_model_dir}/run_subroutines/subroutine.run_hard_FENE");

  // run for N_steps
  lmp->input->one(("run " + std::to_string(N_steps)).c_str());

  // increment Nt
  Nt += N_steps;
}


// run with soft (topoisomerase) potentials and FENE bonds
void LAMMPS_simulator::run_topoDNA_FENE(unsigned long N_steps, thermo_dump_parameters t_d_p)
{
  // setup for run
  setup_run(N_steps,t_d_p);

  // set the timestep
  lmp->input->one("timestep ${delta_t}");

  // include run subroutine
  lmp->input->one("include ${DNA_model_dir}/run_subroutines/subroutine.run_topoDNA_FENE");

  // run for N_steps
  lmp->input->one(("run " + std::to_string(N_steps)).c_str());

  // increment Nt
  Nt += N_steps;
}

// run with soft (topoisomerase) potentials and FENE bonds
void LAMMPS_simulator::run_donothing(unsigned long N_steps, thermo_dump_parameters t_d_p)
{
    // setup for run
    setup_run(N_steps,t_d_p);

    // set the timestep
    lmp->input->one("timestep ${delta_t}");

    // include run subroutine
    lmp->input->one("include ${DNA_model_dir}/run_subroutines/subroutine.run_donothing");

    // run for N_steps
    lmp->input->one(("run " + std::to_string(N_steps)).c_str());

    // increment Nt
    Nt += N_steps;
}


// prepare the dump
void LAMMPS_simulator::prepare_dump(unsigned long N_steps, thermo_dump_parameters &t_d_p)
{
  std::string dump_cmd;
  std::string dump_label = "lammpstrj";

  undump(dump_label);

  // std::cout << "dump_freq = " << t_d_p.dump_freq << std::endl;
  // lmp->input->one("print ${D_freq}");
  if ((dumps[dump_label] == false) && (t_d_p.dump_freq > 0))
    {
      // reset the number of timesteps to Nt
      lmp->input->one(("reset_timestep " + std::to_string(Nt)).c_str());

      // lmp->input->one(("variable skip_condition equal \"step > "+ std::to_string(Nt) + "\"").c_str());

      // delay by an amount corresponding to when the previous dump occurred
      lmp->input->one(("variable D_delay equal "+ std::to_string(prev_dump_Nt + t_d_p.dump_freq)).c_str());

      // include dump
      if (t_d_p.append == true)
	{
	  // append to existing file without first timestep
	  dump_cmd = "include ${DNA_model_dir}/dump_subroutines/subroutine.dump_append_nofirst_";
	}
      else
	{
	  if (t_d_p.write_first == true)
	    {
	      // create new file with first timestep
	      dump_cmd = "include ${DNA_model_dir}/dump_subroutines/subroutine.dump_noappend_first_";
	    }
	  else
	    {
	      // create new file without existing timestep
	      dump_cmd = "include ${DNA_model_dir}/dump_subroutines/subroutine.dump_noappend_nofirst_";
	    }
	}

      // add the dump type label to the dump command
      dump_cmd += dump_label;
      // input the dump command
      lmp->input->one(dump_cmd.c_str());

      // determine how many dumps will occur during the run
      unsigned long N_scheduled_dumps;
      N_scheduled_dumps = (N_steps + Nt - prev_dump_Nt)/t_d_p.dump_freq;
      // update the timestep of the most recent dump
      prev_dump_Nt += N_scheduled_dumps*t_d_p.dump_freq;

      // trigger that dumps are now active
      dumps[dump_label] = true;
    }
}


// reset the simulation timestep to Nt
void LAMMPS_simulator::reset_timestep_to_Nt()
{
  // iterate over the dumps and undump any before resetting the timestep
  for (auto dump : dumps)
    {
      undump(dump.first);
    }

  // reset the number of timesteps to Nt
  lmp->input->one(("reset_timestep " + std::to_string(Nt)).c_str());
}


// reset the timestep counter, Nt
void LAMMPS_simulator::reset_Nt(unsigned long Nt)
{
  this->Nt = Nt;
  reset_timestep_to_Nt();
}


// reset the previous dump timestep counter, prev_dump_Nt
void LAMMPS_simulator::reset_prev_dump_Nt(unsigned long prev_dump_Nt)
{
  this->prev_dump_Nt = prev_dump_Nt;
}


// store the timestep counter, Nt
void LAMMPS_simulator::store_Nt()
{
  stored_Nt = Nt;
}


// restore the timestep counter to the stored value, stored_Nt
void LAMMPS_simulator::restore_Nt()
{
  Nt = stored_Nt;
  reset_timestep_to_Nt();
}


// get atom counts from the simulator and resize the system
void LAMMPS_simulator::sim_to_sys_atom_counts()
{

  double Nd = lammps_get_natoms(lmp);
  int N = int(Nd);
  std::cout << "Nd = " << Nd << std::endl;
  std::cout << "N = " << N << std::endl;

  // determine the sizes of the subarrays
  // unsigned long int *types = new unsigned long int[N];

  // 1 for per-atom type, 1 for size of data (t)
  // lammps_gather_atoms(lmp, const_cast<char*>("type"), 1, 1, types);

  // get the types to match up the quats
  void *types_p;
  // 1 for LMP_STYLE_ATOM, 1 for LMP_TYPE_VECTOR
  types_p = lammps_extract_compute(lmp,const_cast<char*>("type_track"),1,1);
  double *types{static_cast<double*>(types_p)};

  int N_mono = 0;
  int N_ribo = 0;
  int N_bdry = 0;

  int t;
  for (int i=0; i<N; i++)
    {
      // t = static_cast<int>(types[i]);
      t = int(types[i]);
      if (i%50 == 0) std::cout << i << "\t" << t << std::endl;
      if (t == 1)
	{
	  N_bdry += 1;
	}
      else if (t == 2)
	{
	  N_ribo += 1;
	}
      else
	{
	  N_mono += 1;
	}
    }

  int N_mono_ribo = N_mono + N_ribo;
  // std::cout << "N_mono = " << N_mono << std::endl;
  // std::cout << "N_ribo = " << N_ribo << std::endl;
  // std::cout << "N_bdry = " << N_bdry << std::endl;
  // std::cout << "N_mono_ribo = " << N_mono_ribo << std::endl;

  delete[] types;
  std::cout << "past delete" << std::endl;

  // resize the system state based on the simulator
  lmp_sys->set_N_total(N);
  lmp_sys->set_N_mono(N_mono);
  lmp_sys->set_N_ribo(N_ribo);
  lmp_sys->set_N_bdry(N_bdry);
  lmp_sys->set_N_total_ellipsoids(N_mono_ribo);
  lmp_sys->set_N_mono_ellipsoids(N_mono);
  lmp_sys->set_N_ribo_ellipsoids(N_ribo);

}


// dump the simulation state to the system state
void LAMMPS_simulator::sim_to_sys()
{

  // sim_to_sys_atom_counts();

  int N = lmp_sys->get_N_total();
  int N_mono = lmp_sys->get_N_mono();
  int N_ribo = lmp_sys->get_N_ribo();
  int N_bdry = lmp_sys->get_N_bdry();
  int N_mono_ribo = 0;
  if (N_mono > 0) N_mono_ribo += N_mono;
  if (N_ribo > 0) N_mono_ribo += N_ribo;

  std::cout << "Doing sim_to_sys:" << std::endl;
  std::cout << "N = " << N << std::endl;
  std::cout << "N_mono = " << N_mono << std::endl;
  std::cout << "N_ribo = " << N_ribo << std::endl;
  std::cout << "N_bdry = " << N_bdry << std::endl;
  std::cout << "N_mono_ribo = " << N_mono_ribo << std::endl;

  if (N > 0)
    {

      // copy the coordinates to the system state

      double *coords = nullptr;
      if (coords == nullptr) coords = new double[3*N];

      // 1 for per-atom type, 3 for size of data (x_i,x_j,x_k)
      lammps_gather_atoms(lmp, const_cast<char*>("x"), 1, 3, coords);
      lmp_sys->set_coords_arr_total(coords,"row");

      // free the coordinate array
      if (coords != nullptr)
	{
	  delete[] coords;
	  coords = nullptr;
	}

//      // copy the quaternions to the system state
//
//      // get the quats from a compute
//      void *quats_p;
//      // 1 for LMP_STYLE_ATOM, 2 for LMP_TYPE_ARRAY
//      quats_p = lammps_extract_compute(lmp,const_cast<char*>("quat"),1,2);
//      double **quats_2d{static_cast<double**>(quats_p)};
//
//      // get the ids to match up the quats
//      void *ids_p;
//      // 1 for LMP_STYLE_ATOM, 1 for LMP_TYPE_VECTOR
//      ids_p = lammps_extract_compute(lmp,const_cast<char*>("id_track"),1,1);
//      double *ids{static_cast<double*>(ids_p)};
//
//      double *quats = nullptr;
//      if (quats == nullptr) quats = new double[4*N_mono_ribo];
//
//      int id;
//      for (int i=0; i<N; i++)
//	{
//	  id = int(ids[i]);
//	  if (id <= N_mono_ribo)
//	    {
//	      for (int j=0; j<4; j++)
//		{
//		  quats[4*(id-1)+j] = quats_2d[i][j];
//		}
//	    }
//	}
//
//      lmp_sys->set_quats_arr_total(quats,"row");
//
//      // free the quaternion array
//      if (quats != nullptr)
//	{
//	  delete[] quats;
//	  quats = nullptr;
//	}

      // sync the subarrays with the total array

      lmp_sys->sync_subarrays();

    }
}


// read the loop params
int LAMMPS_simulator::read_loop_params(std::string loop_param_filename)
{

  std::fstream loop_param_file;

  std::string param_delim, param, val;
  int delim;

  std::string line;

  loop_sim_params l_sim_p;
  loop_sys_params l_sys_p;

  param_delim = "=";

  loop_param_file.open(loop_param_filename, std::ios::in);

  if (!loop_param_file.is_open())
    {
      std::cout << "ERROR: file not opened in loop_param_file" << std::endl;
      return 1;
    }
  else
    {
      while (1)
	{
	  loop_param_file >> line;
	  if (loop_param_file.eof()) break;


	  if ((line.length() > 0) &&
	      (line.find("#") != 0))
	    {

	      delim = line.find(param_delim);

	      if (delim != -1)
		{

		  param = line.substr(0,delim);
		  val = line.substr(delim+1,line.length());

		  // std::cout << param << "=" << val << std::endl;

		  if (param == "min_dist")
		    {
		      l_sys_p.min_dist = stoi(val);
		    }

		  else if (param == "family")
		    {
		      l_sys_p.family = val;
		    }

		  else if (param == "ext_avg")
		    {
		      l_sys_p.ext_avg = stod(val);
		    }

		  else if (param == "ext_max")
		    {
		      l_sys_p.ext_max = stoi(val);
		    }

		  else if (param == "p_unbinding")
		    {
		      l_sys_p.p_unbinding = stod(val);
		    }

		  else if (param == "r_g")
		    {
		      l_sys_p.r_g = stod(val);
		    }

          else if (param == "k_off")
		    {
		      l_sim_p.k_off = stod(val);
		    }

          else if (param == "k_on")
		    {
		      l_sim_p.k_on = stod(val);
		    }

		  else if (param == "r_0")
		    {
		      l_sim_p.r_0 = stod(val);
		    }

		  else if (param == "k")
		    {
		      l_sim_p.k = stoi(val);
		    }

		  else if (param == "freq_loop")
		    {
		      l_sim_p.freq_loop = stoul(val);
		    }

		  else if (param == "freq_topo")
		    {
		      l_sim_p.freq_topo = stoul(val);
		    }

          else if (param == "freq_replicate")
		    {
		      l_sim_p.freq_replicate = stoul(val);
		    }

		  else if (param == "freq_grow")
		    {
		      l_sim_p.freq_grow = stoul(val);
		    }

		  else if (param == "dNt_topo")
		    {
		      l_sim_p.dNt_topo = stoul(val);
		    }

		}

	    }

     	} // end while loop

      loop_param_file.close();

      set_loop_sim_params(l_sim_p);

      l_sys_p.p_off = l_sim_p.k_off*0.4;
      l_sys_p.p_on = l_sim_p.k_on*0.4;
      lmp_sys->set_loop_sys_params(l_sys_p);

      return 0;

    }

}


// set the loop sim parameters
void LAMMPS_simulator::set_loop_sim_params(loop_sim_params &l_sim_p)
{
  this->l_sim_p = l_sim_p;
}


// delete loops and/or extruded beads bonds
// we should never have to use this function
void LAMMPS_simulator::delete_loops(bool extruded_beads, loop_simulator &sim_ls)
{
	if (extruded_beads == false)
	{
		lmp->input->one("delete_bonds DNA bond 2 remove");
		std::cout << "Deleting anchor hinge bonds..." << std::endl;
	}
	else
	{
		lmp->input->one("delete_bonds DNA bond 3 remove");
		std::cout << "Deleting loop scrunch bonds..." << std::endl;

	}

	int N = lmp_sys->get_N_total();
	int *types = new int[N];
	lmp_sys->get_types(types);
	std::vector<bond> loop_bonds = sim_ls.get_loop_bonds();
	std::vector<std::deque<int>> scrunched_beads = sim_ls.get_extruded();
	// Outer loop over anchor-hinge pairs
    for (size_t i_loop = 0; i_loop < loop_bonds.size(); i_loop++)
	{
		size_t anchor = loop_bonds[i_loop].i;
		size_t hinge = loop_bonds[i_loop].j;

		if (anchor == hinge) {
			hinge +=1 ; //safeguard
		}

	    if ((types[loop_bonds[i_loop].i-1] != 4) && (types[loop_bonds[i_loop].i-1] != 5))
		{
		  types[loop_bonds[i_loop].i-1] = 3; // anchor atom -> normal
		}
	    if ((types[loop_bonds[i_loop].j-1] != 4) && (types[loop_bonds[i_loop].j-1] != 5))
		{
		  types[loop_bonds[i_loop].j-1] = 3; // hinge atom -> normal
		}

  	}

	for (size_t i_loop = 0; i_loop < scrunched_beads.size(); i_loop++)
	{
		for (size_t j_loop = 1; j_loop < scrunched_beads[i_loop].size(); j_loop++) {
			if ((types[scrunched_beads[i_loop][j_loop]-1] != 4) && (types[scrunched_beads[i_loop][j_loop]-1] != 5))
			{
				types[scrunched_beads[i_loop][j_loop]-1] = 3; // scrunched -> normal
			}
		}
	}
	lammps_scatter_atoms(lmp, const_cast<char*>("type"), 0, 1, types);
	delete[] types;

}


void LAMMPS_simulator::form_loops(bool extruded_beads, loop_simulator &sim_ls)
{

	int N = lmp_sys->get_N_total();
	int *types = new int[N];
	lmp_sys->get_types(types);
	std::vector<bond> loop_bonds = sim_ls.get_loop_bonds();
	std::vector<std::deque<int>> scrunched_beads = sim_ls.get_extruded();
	std::unordered_map<int, int> bond_counts;
	std::cout << "N_loop_bonds = " << loop_bonds.size() << std::endl;
	std::string temp_bond_command = "create_bonds single/bond";
	std::string bond_command;
	const int max_bonds_per_atom = 4;  // or whatever your max is
	int scrunch_every = 10;
	// loop over anchor-hinge pairs
	for (size_t i_loop = 0; i_loop < loop_bonds.size(); i_loop++)
	{
		std::cout << "creating bond for loop " << i_loop << std::endl;
		size_t anchor = loop_bonds[i_loop].i;
		size_t hinge = loop_bonds[i_loop].j;
		int bond_type = loop_bonds[i_loop].type;

		if (anchor == hinge) {
			hinge +=1 ; //safeguard
		}

		// Check bond counts first
		if (bond_counts[anchor] < max_bonds_per_atom && bond_counts[hinge] < max_bonds_per_atom)
		{
			bond_command = temp_bond_command;
			bond_command += (" " + std::to_string(bond_type));
			bond_command += (" " + std::to_string(anchor));
			bond_command += (" " + std::to_string(hinge));

			if (i_loop == loop_bonds.size() - 1 && !extruded_beads)
			{
				// the next user of the special list in a hook is the special_bonds rebuild of the BD script (the fused
				// minimiser builds its own exclusions), so the list is marked stale instead of rebuilt here; wcm_special_fresh()
				// rebuilds it on demand. WCM_STALE_SPECIAL_OFF=1: rebuilt here as before.
				static const bool off_stale_special = std::getenv("WCM_STALE_SPECIAL_OFF") != nullptr;
				if (off_stale_special) bond_command += " special yes";
				else { bond_command += " special no"; wcm_special_stale = true; }
			}
			else
			{
				bond_command += " special no";
			}

			lmp->input->one(bond_command);

			// After successful bond creation, increment the counts
			bond_counts[anchor]++;
			bond_counts[hinge]++;
		}
		else
		{
			// Skip bond creation
			std::cerr << "Skipping bond between " << anchor << " and " << hinge
					  << " because one of them reached bond limit.\n";
		}


		if ((types[loop_bonds[i_loop].i-1] != 4) && (types[loop_bonds[i_loop].i-1] != 5))
		{
			types[loop_bonds[i_loop].i-1] = 7; // anchor atom -> normal
		}
		if ((types[loop_bonds[i_loop].j-1] != 4) && (types[loop_bonds[i_loop].j-1] != 5))
		{
			types[loop_bonds[i_loop].j-1] = 8; // hinge atom -> normal
		}

	}

    if (extruded_beads) {
	// loop over scrunched beads
	for (size_t i_loop = 0; i_loop < scrunched_beads.size(); i_loop++)
	{
		size_t anchor = loop_bonds[i_loop].i;
		size_t hinge = loop_bonds[i_loop].j;

		if (anchor == hinge) {
			hinge +=1 ; //safeguard
		}

		for (size_t j_loop = scrunch_every; j_loop < scrunched_beads[i_loop].size(); j_loop += scrunch_every) {


			int current_bead = scrunched_beads[i_loop][j_loop];
			int previous_bead = (j_loop > 0) ? scrunched_beads[i_loop][j_loop - scrunch_every] : anchor;

			// Check bond limits first
			if (bond_counts[current_bead] < max_bonds_per_atom && bond_counts[previous_bead] < max_bonds_per_atom)
			{
				std::cout << "Creating bond for loop " << i_loop
						  << ", extruded bead " << j_loop
						  << " between " << current_bead
						  << " and " << previous_bead << std::endl;

				bond_command = temp_bond_command;
				bond_command += (" " + std::to_string(3));  // bond type 3?
				bond_command += (" " + std::to_string(current_bead));
				bond_command += (" " + std::to_string(previous_bead));

				if (i_loop == scrunched_beads.size() - 1 && j_loop +scrunch_every >= scrunched_beads[i_loop].size())
				{
					bond_command += " special yes";
				}
				else
				{
					bond_command += " special no";
				}

				lmp->input->one(bond_command);

				// Update bond counts
				bond_counts[current_bead]++;
				bond_counts[previous_bead]++;
			}
			else
			{
				std::cerr << "Skipping bond between " << current_bead
						  << " and " << previous_bead
						  << " because one of them reached bond limit.\n";
			}


			if ((types[scrunched_beads[i_loop][j_loop]-1] != 4) && (types[scrunched_beads[i_loop][j_loop]-1] != 5)
			&& (types[scrunched_beads[i_loop][j_loop]-1] != 7) && (types[scrunched_beads[i_loop][j_loop]-1] != 8))
			{
				types[scrunched_beads[i_loop][j_loop]-1] = 9; // scrunched
			}
		}

	}
	}
	lammps_scatter_atoms(lmp, const_cast<char*>("type"), 0, 1, types);
	delete[] types;

	if (extruded_beads) {
		lmp->input->one("include ${DNA_model_dir}/protocol_subroutines/subroutine.group_setup");
		std::cout << "Deleting default DNA bonds" << std::endl;
		lmp->input->one("delete_bonds extruded_monos bond 1 remove");
	}

}


// update the loop bonds
void LAMMPS_simulator::update_loop_bonds(bool new_bonds, loop_simulator &sim_ls)
{

  // type array for scatter
  int N = lmp_sys->get_N_total();
  int *types = new int[N];

  // 1 for per-atom type, 1 for size of data (t)
  // lammps_gather_atoms(lmp, const_cast<char*>("type"), 1, 1, types);

  lmp_sys->get_types(types);

  // get the loop bonds
  // std::vector<bond> loop_bonds = lmp_sys->get_loop_bonds();
  std::vector<bond> loop_bonds = sim_ls.get_loop_bonds();
  update_loop_bonds_apply(new_bonds, loop_bonds, types);
  delete[] types;
}


// defer-loop-bonds verification (WCM_DEFER_LOOP_BONDS_VERIFY=1): what a clear+read sets up, compared byte for byte before and after a rebuild that
// the driver would have skipped. Host-side arrays (current right after read_data: no run in between).
std::string LAMMPS_simulator::wcm_state_fingerprint()
{
  std::string out;
  auto add = [&out](const void *p, size_t n) { out.append(static_cast<const char*>(p), n); };
  LAMMPS_NS::Atom *a = lmp->atom;
  const int nlocal = a->nlocal;
  add(&nlocal, sizeof(int));
  add(&a->nbonds, sizeof(a->nbonds)); add(&a->nangles, sizeof(a->nangles));
  add(a->tag, sizeof(*a->tag) * nlocal);                    // local atom order
  add(&a->x[0][0], sizeof(double) * 3 * nlocal);
  add(a->type, sizeof(int) * nlocal);
  add(a->mask, sizeof(int) * nlocal);
  add(a->image, sizeof(*a->image) * nlocal);
  if (a->num_bond) for (int i = 0; i < nlocal; i++)
    { add(&a->num_bond[i], sizeof(int)); add(a->bond_type[i], sizeof(int) * a->num_bond[i]); add(a->bond_atom[i], sizeof(*a->bond_atom[i]) * a->num_bond[i]); }
  if (a->num_angle) for (int i = 0; i < nlocal; i++)
    { add(&a->num_angle[i], sizeof(int)); add(a->angle_type[i], sizeof(int) * a->num_angle[i]);
      add(a->angle_atom1[i], sizeof(*a->angle_atom1[i]) * a->num_angle[i]); add(a->angle_atom2[i], sizeof(*a->angle_atom2[i]) * a->num_angle[i]);
      add(a->angle_atom3[i], sizeof(*a->angle_atom3[i]) * a->num_angle[i]); }
  if (a->nspecial) for (int i = 0; i < nlocal; i++)
    { add(a->nspecial[i], sizeof(int) * 3); add(a->special[i], sizeof(*a->special[i]) * a->nspecial[i][2]); }
  LAMMPS_NS::Group *g = lmp->group;
  for (int k = 0; k < g->ngroup; k++) { if (g->names[k]) out += g->names[k]; out += ';'; }
  return out;
}


// fast-read-atoms verification: what a read leaves in the per-atom arrays the Atoms section fills (the defer-loop-bonds fingerprint + nmax, molecule, v)
std::string LAMMPS_simulator::wcm_fast_read_atoms_fingerprint()
{
  std::string out = wcm_state_fingerprint();
  LAMMPS_NS::Atom *a = lmp->atom;
  const int nlocal = a->nlocal;
  out.append(reinterpret_cast<const char *>(&a->nmax), sizeof(int));
  if (a->molecule) out.append(reinterpret_cast<const char *>(a->molecule), sizeof(*a->molecule) * nlocal);
  if (a->v) out.append(reinterpret_cast<const char *>(&a->v[0][0]), sizeof(double) * 3 * nlocal);
  return out;
}


// rebuild the special list (1-2/1-3/1-4 partners) with the current special_bonds settings if form_loops left it stale.
// "delete_bonds all stats special" deletes nothing and ends with the same Special::build that "special yes" on create_bonds runs.
void LAMMPS_simulator::wcm_special_fresh()
{
  if (!wcm_special_stale) return;
  std::cout << "rebuilding the special list left stale by form_loops" << std::endl;
  lmp->input->one("delete_bonds all stats special");
  wcm_special_stale = false;
}

unsigned long long LAMMPS_simulator::wcm_topology_hash()
{
  unsigned long long h = 1469598103934665603ULL;
  auto mix = [&h](const void *p, size_t n) { const unsigned char *c = static_cast<const unsigned char*>(p); for (size_t k = 0; k < n; k++) { h ^= c[k]; h *= 1099511628211ULL; } };
  LAMMPS_NS::Atom *a = lmp->atom;
  const int nlocal = a->nlocal;
  mix(&nlocal, sizeof(int));
  mix(a->tag, sizeof(*a->tag) * nlocal); mix(a->type, sizeof(int) * nlocal); mix(a->mask, sizeof(int) * nlocal);
  if (a->num_bond) for (int i = 0; i < nlocal; i++)
    { mix(&a->num_bond[i], sizeof(int)); mix(a->bond_type[i], sizeof(int) * a->num_bond[i]); mix(a->bond_atom[i], sizeof(*a->bond_atom[i]) * a->num_bond[i]); }
  if (a->num_angle) for (int i = 0; i < nlocal; i++)
    { mix(&a->num_angle[i], sizeof(int)); mix(a->angle_atom1[i], sizeof(*a->angle_atom1[i]) * a->num_angle[i]); mix(a->angle_atom3[i], sizeof(*a->angle_atom3[i]) * a->num_angle[i]); }
  if (a->nspecial) for (int i = 0; i < nlocal; i++)
    { mix(a->nspecial[i], sizeof(int) * 3); mix(a->special[i], sizeof(*a->special[i]) * a->nspecial[i][2]); }
  return h;
}


// everything update_loop_bonds does to LAMMPS, from its two CPU-side inputs (the loop bonds of the loop simulator and
// the system's per-atom types), so that the driver can capture the inputs at load_loops and apply them later.
void LAMMPS_simulator::update_loop_bonds_apply(bool new_bonds, const std::vector<bond> &loop_bonds, int *types)
{
  // delete the existing loop bonds
  if (new_bonds == false)
    {
      lmp->input->one("delete_bonds DNA bond 2 remove");
    }

  std::cout << "N_loop_bonds = " << loop_bonds.size() << std::endl;
  // add the updated loop bonds
  std::string temp_bond_command = "create_bonds single/bond";
  std::string bond_command;
	// Outer loop over anchor-hinge pairs
  for (size_t i_loop = 0; i_loop < loop_bonds.size(); i_loop++)
	{
		size_t anchor = loop_bonds[i_loop].i;
		size_t hinge = loop_bonds[i_loop].j;
		int bond_type = loop_bonds[i_loop].type;

		if (anchor == hinge) {
			hinge +=1 ; //safeguard
		}
		// int N_extruded_mono = abs(loop_bonds[i_loop].j - loop_bonds[i_loop].i);
  		// int dir = (loop_bonds[i_loop].j - loop_bonds[i_loop].i) / N_extruded_mono;
  	    // bool passed_origin = N_extruded_mono > 40000; // we can use the ter_crossing logic
  	    // if (passed_origin) {
  		// std::cout << "Anchor hinge pair crosses terminus..." << std::endl;
  		// N_extruded_mono = 54338 - N_extruded_mono;
  		//dir = -1 * dir;
  	    // }
		// Determine the number of beads between the anchor and hinge
		// size_t num_beads = N_extruded_mono + 1;
		// size_t num_pairs = num_beads / 2; // Determines how many bonds to create
		// num_pairs = 1;

		// Inner loop to "zip up" the loop
		// for (size_t k = 0; k < num_pairs; k++)
		// {
			// size_t i = anchor + dir * k; // Moving from anchor inward
			// size_t j = hinge - dir * k;  // Moving from hinge inward
			// If passed_origin is true, handle modular wrapping
			// if (passed_origin) {
				// i = (i - 1 + 54338) % 54338 + 1; // Wrap around
				// j = (j -1 + 54338) % 54338 + 1; // Wrap around
			// }

			bond_command = temp_bond_command;
			bond_command += (" " + std::to_string(bond_type));
			bond_command += (" " + std::to_string(anchor));
			bond_command += (" " + std::to_string(hinge));

			// Mark "special yes" only for the last bond of the last loop iteration
			// if (i_loop == loop_bonds.size() - 1 && k == num_pairs - 1)
  	        if (i_loop == loop_bonds.size() - 1)
			{
				bond_command += " special yes";
			}
			else
			{
				bond_command += " special no";
			}

			lmp->input->one(bond_command);


      if ((types[loop_bonds[i_loop].i-1] != 4) && (types[loop_bonds[i_loop].i-1] != 5))
	{
	  types[loop_bonds[i_loop].i-1] = 7; // anchor atom
	}
      if ((types[loop_bonds[i_loop].j-1] != 4) && (types[loop_bonds[i_loop].j-1] != 5))
	{
	  types[loop_bonds[i_loop].j-1] = 8; // hinge atom
	}



	std::cout << "Monomers " << loop_bonds[i_loop].j << " (hinge) through " << loop_bonds[i_loop].i << " (anchor) are extruded" << std::endl;
  	// for (int dist=0; dist<N_extruded_mono; dist++) {
  	// 	if ((types[loop_bonds[i_loop].i+dir*dist-1] != 4)
  	//	&& (types[loop_bonds[i_loop].i+dir*dist-1] != 5)
  	//	&& (types[loop_bonds[i_loop].i+dir*dist-1] != 7)
  	//	&& (types[loop_bonds[i_loop].i+dir*dist-1] != 8))
  	//	{

	//		int index = loop_bonds[i_loop].i + dir * dist - 1;

    //        // If passed_origin is true, handle modular wrapping
    //		if (passed_origin) {
    //    		index = (index + 54338) % 54338; // Wrap around
    //		}
  			// types[index] = 9; // extruded monomer

  	//	}
    // }


  }




  // scatter the now modified atom types
  // 0 for integer type, 1 for per-atom count
  lammps_scatter_atoms(lmp, const_cast<char*>("type"), 0, 1, types);
}


// run a system with loops
void LAMMPS_simulator::run_loops(int N_loops, double threshold, unsigned long N_steps, thermo_dump_parameters t_d_p,  btree &sim_bt, mapper &sim_mapper, loop_simulator &sim_ls)
{

  unsigned long step_counter = 0;
  unsigned long step_increment;
  unsigned long step_next_grow = l_sim_p.freq_grow;
  unsigned long step_next_replicate = l_sim_p.freq_replicate;
  unsigned long step_next_topo = l_sim_p.freq_topo;
  unsigned long step_next_loop = l_sim_p.freq_loop;

  std::cout << "Simulation Steps:\n";
    std::cout << std::left << std::setw(20) << "Step Grow:"
              << step_next_grow << "\n";
    std::cout << std::left << std::setw(20) << "Step Replicate:"
              << step_next_replicate << "\n";
    std::cout << std::left << std::setw(20) << "Step Topo:"
              << step_next_topo << "\n";
    std::cout << std::left << std::setw(20) << "Step Loop:"
              << step_next_loop << "\n";
  thermo_dump_parameters t_d_p_iter = t_d_p;
  thermo_dump_parameters t_d_p_topo = t_d_p;
  bool first_iteration, new_bonds, replicated;
  double radius = 2000.0;

  // flag for first iteration
  first_iteration = true;

  // prepare the dummy thermo and dump info
  t_d_p_topo.dump_freq = 0;

  // Equilibrate looping state:
  std::cout << "---[ Equilibrating the 50 loops... " << step_counter << " ]---" << std::endl;
  // sim_ls.steps(1000000);  // 0.4 seconds, 1 kbp = 100 beads per second

  // repeat until the final number of steps is reached
  while (step_counter < N_steps)
    {
      replicated = false;
      std::cout << "---[ simulating timestep " << step_counter << " ]---" << std::endl;
      // set the system state to the current simulator state
      sim_to_sys();

      // check if first iteration and perform specific actions
      if (first_iteration == true)
	{
	  // lmp_sys->initialize_loop_topo(N_loops, threshold);
	  // enable the initialization of anchors on updating
	  new_bonds = true;
	  first_iteration = false;
	}
      else
	{
	  // disable the re-initialization of anchors on updating
	  new_bonds = false;
	  // modify thermo and dump info for subsequent iterations
	  t_d_p_iter.write_first = false;
	  t_d_p_iter.append = true;
	}

    // grow the cell by increasing cell volume, adding ribosomes, and adding loops
      if (step_counter == step_next_grow)
    {
      radius = 2000.0 + 160.0 + int(500.0 * double(step_counter) / double(N_steps)); // assume that we run for dounbling time of 1 hour.
      threshold = radius + 50.0;
	  // expand boundary radius
      std::cout << "---[ growing cell volume to radius " << radius << " ]---" << std::endl;
      auto start_time = std::chrono::high_resolution_clock::now();
      lmp_sys->generate_spherical_bdry(radius,
  					 0.0,
  					 0.0,
  					 0.0);
      auto end_time = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
      std::cout << "Execution time: " << duration << " milliseconds" << std::endl;

      // add ribosomes
      std::cout << "---[ adding ribosomes ]---" << std::endl;
      lmp_sys->add_ribos(5, radius);

      // add a loop
      std::cout << "---[ adding loops ]---" << std::endl;
      // lmp_sys->add_loops(1, threshold);

      // have system write data and simulator read data
      lmp_sys->write_data("/home/andrew/Data/btree_chromo/DNA_membrane_dynamics/data/btree_chromo/data" + output_file_label + ".lammps");
      clear();
      reset_protocol_variables();
      global_setup();
      read_data("/home/andrew/Data/btree_chromo/DNA_membrane_dynamics/data/btree_chromo/data" + output_file_label + ".lammps");
      prepare_fork_partition_groups(1);
      standard_computes();

      step_next_grow = step_counter + l_sim_p.freq_grow;
    }

      // replicate the chromosome
      if (step_counter == step_next_replicate)
    {
       sim_mapper.set_initial_state(sim_bt.get_state());
       sim_bt.single_transform(sim_bt.parse_transform("m_cw20_ccw20"));
       sim_mapper.set_final_state(sim_bt.get_state());
       sim_mapper.prepare_mapping();
       lmp_sys->apply_mono_mapping(sim_mapper.get_map());
       sim_bt.write_state("/home/andrew/Data/btree_chromo/DNA_membrane_dynamics/params/btree_chromo/rep_state" + output_file_label + ".txt",sim_bt.get_state());
       lmp_sys->set_btree(sim_bt.get_state());

       // have system write data and simulator read data
       lmp_sys->write_data("/home/andrew/Data/btree_chromo/DNA_membrane_dynamics/data/btree_chromo/data" + output_file_label + ".lammps");
       clear();
       reset_protocol_variables();
       global_setup();
       read_data("/home/andrew/Data/btree_chromo/DNA_membrane_dynamics/data/btree_chromo/data" + output_file_label + ".lammps");
       prepare_fork_partition_groups(1);
       standard_computes();

       // update loop regions to be used for future anchor-hinge binding
       // lmp_sys->prepare_loop_topo();
       // lmp_sys->update_loop_binding();
       step_next_replicate = step_counter + l_sim_p.freq_replicate;
       new_bonds = true;
       replicated = true;
    }

      // extrude loops
      if (step_counter == step_next_loop)
    {
        if (new_bonds == false || replicated == true) // update if bonds are not new or if we just replicated
      {
        std::cout << "---[ updating loop topology ]---" << std::endl;
        // update the loop topology (unbind existing bound anchors, bind new anchors, extrude loops for bound anchors)
        // lmp_sys->update_loop_topo(threshold);
        std::cout << N_loops << "+" << threshold << std::endl;
      }
      // update the loop bonds
      std::cout << "---[ updating loop bonds ]---" << std::endl;
      // sim_ls.steps(40);  // 0.4 seconds, 1 kbp = 100 beads per second
      sim_ls.write_state("/home/andrew/Data/btree_chromo/DNA_membrane_dynamics/params/btree_chromo/loop_state" + output_file_label + std::to_string(step_counter) + ".txt");
      update_loop_bonds(new_bonds, sim_ls);
      step_next_loop = step_counter + l_sim_p.freq_loop;
    }

      // simulate topoisomerase action
      if (step_counter == step_next_topo)
	{
      auto start_time = std::chrono::high_resolution_clock::now();
      minimize_topoDNA_harmonic(t_d_p_topo);
      auto end_time = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
      std::cout << "Execution time: " << duration << " milliseconds" << std::endl;


	  start_time = std::chrono::high_resolution_clock::now();
      run_topoDNA_harmonic(l_sim_p.dNt_topo,t_d_p_topo);
      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
      std::cout << "Execution time: " << duration << " milliseconds" << std::endl;

      step_next_topo = step_counter + l_sim_p.freq_topo;
	}

      // determine the number of steps to be simulated
      step_increment = std::min(step_next_grow-step_counter,
                       std::min(step_next_replicate-step_counter,
                       std::min(step_next_topo-step_counter,
                       std::min(step_next_loop-step_counter,
                                N_steps-step_counter))));
  	  // minimize_hard_harmonic(t_d_p_topo);
  	  // minimize_soft_harmonic(t_d_p_topo);
      // try {
      // run dynamics
      auto start_time = std::chrono::high_resolution_clock::now();
  	  // run_hard_harmonic(step_increment,t_d_p_iter);
  	  run_hard_harmonic(step_increment,t_d_p_iter);
      // run_donothing(step_increment,t_d_p_iter);
      auto end_time = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
      std::cout << "Execution time: " << duration << " milliseconds" << std::endl;
      step_counter += step_increment;
      // } catch (const std::exception& e) {
      // Handle the error and execute fallback code
      //std::cerr << "Error occurred: " << e.what() << std::endl;
      //lmp_sys->write_mono_coords("/home/andrew/Data/btree_chromo/DNA_membrane_dynamics/scripts/../data/btree_chromo/dna_monomers/dna_1.bin", "row");
      //lmp_sys->write_ribo_coords("/home/andrew/Data/btree_chromo/DNA_membrane_dynamics/scripts/../data/btree_chromo/dna_monomers/ribo_1.bin", "row");
      //lmp_sys->write_bdry_coords("/home/andrew/Data/btree_chromo/DNA_membrane_dynamics/scripts/../data/btree_chromo/dna_monomers/bdry_1.bin", "row");
      //}

    }

}


// reset all computes
void LAMMPS_simulator::initialize_computes()
{
  // initialize computes
  computes["quats"] = false;
  computes["ids"] = false;
  computes["types"] = false;
  computes["MSD"] = false;
}


// remove the compute
void LAMMPS_simulator::uncompute(std::string compute_label)
{
  if (computes[compute_label] == true)
    {
      if (compute_label == "ids")
	{
	  lmp->input->one("uncompute id_track");
	}
      else if (compute_label == "types")
	{
	  lmp->input->one("uncompute type_track");
	}
      else if (compute_label == "quats")
	{
	  lmp->input->one("uncompute quat");
	}
      else if (compute_label == "MSD")
	{
	  lmp->input->one("uncompute dnaMSD");
	  lmp->input->one("uncompute ribosMSD");
	  lmp->input->one("uncompute orisMSD");
	  lmp->input->one("uncompute tersMSD");
	  lmp->input->one("uncompute forksMSD");
	}
      computes[compute_label] = false;
    }
}


// trigger for preventing reuse of computes
void LAMMPS_simulator::compute_trigger(std::string compute_label)
{
  std::string compute_cmd;

  // remove compute if it already exists
  uncompute(compute_label);

  compute_cmd = "include ${DNA_model_dir}/compute_subroutines/subroutine.compute_";
  compute_cmd += compute_label;

  // input the compute command
  lmp->input->one(compute_cmd.c_str());

  // switch the state of the compute
  computes[compute_label] = true;

}


void LAMMPS_simulator::initialize_dumps()
{
  // initialize dumps
  dumps["lammpstrj"] = false;
}


void LAMMPS_simulator::undump(std::string dump_label)
{
  if (dumps[dump_label] == true)
    {
      lmp->input->one(("undump d_" + dump_label).c_str());
      dumps[dump_label] = false;
    }
}


void LAMMPS_simulator::initialize_sim_vars()
{
  // initialize simulator internal variables
  sim_vars["T_freq"] = false;
  sim_vars["D_freq"] = false;
  sim_vars["ellipsoids"] = false;

  // initialize internal variables for extra potentials
  std::string p;
  for (auto extra_pot=extra_pots.begin(); extra_pot!=extra_pots.end(); ++extra_pot)
    {
      p = extra_pot->first;
      sim_vars[p] = false;
    }
}


void LAMMPS_simulator::set_sim_var_int(std::string sim_var, int val)
{
  std::string sim_var_cmd;

  // delete sim_var if it already exists
  delete_sim_var(sim_var);

  sim_var_cmd = "variable " + sim_var + " internal ";
  sim_var_cmd += std::to_string(val);

  // input the sim_var command
  lmp->input->one(sim_var_cmd.c_str());

  // switch the state of the sim_var
  sim_vars[sim_var] = true;
}


void LAMMPS_simulator::delete_sim_var(std::string sim_var)
{
  if (sim_vars[sim_var] == true)
    {
      lmp->input->one(("variable " + sim_var + " delete").c_str());
      sim_vars[sim_var] = false;
    }
}


void LAMMPS_simulator::initialize_extra_potentials()
{
  extra_pots["Ori_bdry_attraction"] = false;
  extra_pots["DNA_bdry_attraction"] = false;
  extra_pots["Ori_pair_repulsion"] = false;
}


void LAMMPS_simulator::switch_extra_potential(std::string p, bool s)
{
  extra_pots[p] = s;
}


bool LAMMPS_simulator::get_extra_potential_state(std::string p)
{
  return extra_pots[p];
}


void LAMMPS_simulator::extra_pots_to_sim_vars()
{
  std::string p;
  bool v;

  // iterate over the set of extra potentials and set sim_vars accordingly
  for (auto extra_pot=extra_pots.begin(); extra_pot!=extra_pots.end(); ++extra_pot)
    {
      p = extra_pot->first;
      v = extra_pot->second;

      if (v == true)
	{
	  set_sim_var_int(p,1);
	}
      else
	{
	  set_sim_var_int(p,0);
	}
    }
}


void LAMMPS_simulator::initialize_extra_fixes()
{
  extra_fixes["fork_partition_repulsion"] = false;
}


void LAMMPS_simulator::switch_extra_fix(std::string p, bool s)
{

  bool s_old = extra_fixes[p];

  if (s != s_old)
    {

      // use functions for extra fix
      if (p == "fork_partition_repulsion")
	{
	  switch_fork_partition_force(s);
	}

      // switch the extra fix
      extra_fixes[p] = s;

    }
}


bool LAMMPS_simulator::get_extra_fix_state(std::string p)
{
  return extra_fixes[p];
}


void LAMMPS_simulator::disable_and_hold_extra_fixes()
{
  std::string p;
  bool s;

  // iterate over the set of extra fixes - store their state and disable
  for (auto extra_fix=extra_fixes.begin(); extra_fix!=extra_fixes.end(); ++extra_fix)
    {
      p = extra_fix->first;
      s = extra_fix->second;

      extra_fixes_hold[p] = s;
      switch_extra_fix(p,false);
    }
}


void LAMMPS_simulator::restore_extra_fixes()
{
  std::string p;
  bool s;

  // iterate over the set of extra fixes - store their state and disable
  for (auto extra_fix=extra_fixes_hold.begin(); extra_fix!=extra_fixes_hold.end(); ++extra_fix)
    {
      p = extra_fix->first;
      s = extra_fix->second;

      switch_extra_fix(p,s);
    }
}


unsigned long LAMMPS_simulator::get_Nt()
{
  return Nt;
}


void LAMMPS_simulator::switch_ellipsoids(bool s)
{
    if (s)
    {
        set_sim_var_int("ellipsoids", 1);
    }
    else
    {
        set_sim_var_int("ellipsoids", 0);
        sim_vars["ellipsoids"] = false;
    }
}

void LAMMPS_simulator::increment_Nt(unsigned long dNt)
{
  Nt += dNt;
  reset_timestep_to_Nt();
}


void LAMMPS_simulator::scale_bdry_particles(double ds)
{
  std::string cmd;

  cmd = "variable r_bdry_temp equal " + std::to_string(ds) + "*${r_bdry}";
  command(cmd);

  cmd = "variable sigma_mono_bdry equal ${r_mono}+${r_bdry_temp}";
  command(cmd);
  cmd = "variable sigma_ribo_bdry equal ${r_ribo}+${r_bdry_temp}";
  command(cmd);
  cmd = "variable sigma_bdry_bdry equal 2*${r_bdry_temp}";
  command(cmd);

  cmd = "variable WCA_mono_bdry equal ${cut_WCA}*${sigma_mono_bdry}";
  command(cmd);
  cmd = "variable WCA_ribo_bdry equal ${cut_WCA}*${sigma_ribo_bdry}";
  command(cmd);
  cmd = "variable WCA_bdry_bdry equal ${cut_WCA}*${sigma_bdry_bdry}";
  command(cmd);
}


void LAMMPS_simulator::reset_bdry_particle_size()
{
  std::string cmd;

  cmd = "variable sigma_mono_bdry equal ${r_mono}+${r_bdry}";
  command(cmd);
  cmd = "variable sigma_ribo_bdry equal ${r_ribo}+${r_bdry}";
  command(cmd);
  cmd = "variable sigma_bdry_bdry equal 2*${r_bdry}";
  command(cmd);

  cmd = "variable WCA_mono_bdry equal ${cut_WCA}*${sigma_mono_bdry}";
  command(cmd);
  cmd = "variable WCA_ribo_bdry equal ${cut_WCA}*${sigma_ribo_bdry}";
  command(cmd);
  cmd = "variable WCA_bdry_bdry equal ${cut_WCA}*${sigma_bdry_bdry}";
  command(cmd);
}


void LAMMPS_simulator::scale_ribo_particles(double ds)
{
  std::string cmd;

  cmd = "variable r_ribo_temp equal " + std::to_string(ds) + "*${r_ribo}";
  command(cmd);

  cmd = "variable sigma_mono_ribo equal ${r_mono}+${r_ribo_temp}";
  command(cmd);
  cmd = "variable sigma_ribo_bdry equal ${r_ribo_temp}+${r_bdry}";
  command(cmd);
  cmd = "variable sigma_ribo_ribo equal 2*${r_ribo_temp}";
  command(cmd);

  cmd = "variable WCA_mono_ribo equal ${cut_WCA}*${sigma_mono_ribo}";
  command(cmd);
  cmd = "variable WCA_ribo_bdry equal ${cut_WCA}*${sigma_ribo_bdry}";
  command(cmd);
  cmd = "variable WCA_ribo_ribo equal ${cut_WCA}*${sigma_ribo_ribo}";
  command(cmd);
}


void LAMMPS_simulator::reset_ribo_particle_size()
{
  std::string cmd;

  cmd = "variable sigma_mono_ribo equal ${r_mono}+${r_ribo}";
  command(cmd);
  cmd = "variable sigma_ribo_bdry equal ${r_ribo}+${r_bdry}";
  command(cmd);
  cmd = "variable sigma_ribo_ribo equal 2*${r_ribo}";
  command(cmd);

  cmd = "variable WCA_mono_ribo equal ${cut_WCA}*${sigma_mono_ribo}";
  command(cmd);
  cmd = "variable WCA_ribo_bdry equal ${cut_WCA}*${sigma_ribo_bdry}";
  command(cmd);
  cmd = "variable WCA_ribo_ribo equal ${cut_WCA}*${sigma_ribo_ribo}";
  command(cmd);
}


// prepare groups based on fork partitions
void LAMMPS_simulator::prepare_fork_partition_groups(int idx)
{

  command("include ${DNA_model_dir}/potentials/lmp.fork_partitioning");

  std::vector<fork_partition> f_ps = lmp_sys->get_all_fork_partitions();

  std::string temp_range;
  std::string mother, ld, rd;
  std::string group_cmd, variable_cmd, region_cmd;

  std::array<std::string,3> dims = {"x","y","z"};

  int mono_inc = 1;

  for (fork_partition f_p : f_ps)
    {

      mother = f_p.fork;
      ld = mother + "l";
      rd = mother + "r";

      // partitioning of left daughter
      group_cmd = "group " + ld + " id";

      for (mono_range m_r : f_p.left_monos)
	{
	  if (m_r.wrapped == false)
	    {
	      temp_range = " " + std::to_string(m_r.ll+idx) + ":"
		+ std::to_string(m_r.ul+idx) + ":"
		+ std::to_string(mono_inc);
	      group_cmd += temp_range;
	    }
	  else
	    {
	      temp_range = " " + std::to_string(m_r.ll+idx) + ":"
		+ std::to_string(m_r.mid_ll+idx) + ":"
		+ std::to_string(mono_inc);
	      group_cmd += temp_range;
	      temp_range = " " + std::to_string(m_r.mid_ul+idx) + ":"
		+ std::to_string(m_r.ul+idx) + ":"
		+ std::to_string(mono_inc);
	      group_cmd += temp_range;
	    }
	}

      std::cout << group_cmd << std::endl;

      command(group_cmd);

      // partitioning of right daughter
      group_cmd = "group " + rd + " id";

      for (mono_range m_r : f_p.right_monos)
	{
	  if (m_r.wrapped == false)
	    {
	      temp_range = " " + std::to_string(m_r.ll+idx) + ":"
		+ std::to_string(m_r.ul+idx) + ":"
		+ std::to_string(mono_inc);
	      group_cmd += temp_range;
	    }
	  else
	    {
	      temp_range = " " + std::to_string(m_r.ll+idx) + ":"
		+ std::to_string(m_r.mid_ll+idx) + ":"
		+ std::to_string(mono_inc);
	      group_cmd += temp_range;
	      temp_range = " " + std::to_string(m_r.mid_ul+idx) + ":"
		+ std::to_string(m_r.ul+idx) + ":"
		+ std::to_string(mono_inc);
	      group_cmd += temp_range;
	    }
	}

      std::cout << group_cmd << std::endl;

      command(group_cmd);


      // create variables for relative vector between CoMs
      for (size_t i=0; i<dims.size(); i++)
	{
	  variable_cmd = "variable d" + dims[i] + "_com_" + mother;
	  variable_cmd += " equal ";
	  variable_cmd += "xcm(" + ld  + "," + dims[i] + ")-";
	  variable_cmd += "xcm(" + rd + "," + dims[i] +")";

	  std::cout << variable_cmd << std::endl;
	  command(variable_cmd);
	}

      // create variables for count
      variable_cmd = "variable N_" + ld;
      variable_cmd += " equal count(" + ld  + ")";
      std::cout << variable_cmd << std::endl;
      command(variable_cmd);
      variable_cmd = "variable N_" + rd;
      variable_cmd += " equal count(" + rd  + ")";
      std::cout << variable_cmd << std::endl;
      command(variable_cmd);

      // create variables for joint CoM
      for (size_t i=0; i<dims.size(); i++)
	{
	  variable_cmd = "variable " + dims[i] + "_com_" + mother;
	  variable_cmd += " equal ";
	  variable_cmd += "((v_N_" + ld;
	  variable_cmd += "*xcm(" + ld + "," + dims[i] + "))+";
	  variable_cmd += "(v_N_" + rd;
	  variable_cmd += "*xcm(" + rd + "," + dims[i] + ")))";
	  variable_cmd += "/(v_N_" + ld + "+v_N_" + rd + ")";

	  std::cout << variable_cmd << std::endl;
	  command(variable_cmd);
	}

      // create spherical region around joint CoM
      region_cmd = "region sphere_" + mother;
      region_cmd += " sphere";
      for (size_t i=0; i<dims.size(); i++)
	{
	  region_cmd += " ${" + dims[i] + "_com_" + mother + "}";
	}
      region_cmd += " ${fork_cutoff}";
      std::cout << region_cmd << std::endl;
      command(region_cmd);

      // create variable for distance between CoMs
      variable_cmd = "variable d_com_" + mother + " equal ";
      variable_cmd += "sqrt(";
      for (size_t i=0; i<dims.size(); i++)
	{
	  if (i == 0)
	    {
	      variable_cmd += "v_d" + dims[i] + "_com_" + mother
		+ "^2";
	    }
	  else
	    {
	      variable_cmd += "+v_d" + dims[i] + "_com_" + mother
		+ "^2";
	    }
	}

      variable_cmd += ")";

      std::cout << variable_cmd << std::endl;
      command(variable_cmd);

      // create variables for direction between CoMs
      for (size_t i=0; i<dims.size(); i++)
	{
	  variable_cmd = "variable ud" + dims[i] + "_com_" + mother;
	  variable_cmd += " equal ";
	  variable_cmd += "v_d" + dims[i] + "_com_" + mother;
	  variable_cmd += "/v_d_com_" + mother;

	  std::cout << variable_cmd << std::endl;
	  command(variable_cmd);
	}

    }
}


// apply/remove forces to fork partitions
void LAMMPS_simulator::switch_fork_partition_force(bool s)
{

  std::vector<fork_partition> f_ps = lmp_sys->get_all_fork_partitions();

  command("include ${DNA_model_dir}/potentials/lmp.fork_partitioning");

  if (s == true)
    {

      std::string mother, ld, rd;
      std::string fix_cmd, temp_var, variable_cmd;

      std::array<std::string,3> dims = {"x","y","z"};

      for (fork_partition f_p : f_ps)
	{

	  mother = f_p.fork;
	  ld = mother + "l";
	  rd = mother + "r";

	  for (size_t i=0; i<dims.size(); i++)
	    {
	      temp_var = "lf" + dims[i] + "_" + mother;
	      variable_cmd = "variable " + temp_var + " equal ";
	      variable_cmd += "v_fork_force*v_ud" + dims[i]
		+ "_com_" + mother;
	      std::cout << variable_cmd << std::endl;
	      command(variable_cmd);

	      temp_var = "rf" + dims[i] + "_" + mother;
	      variable_cmd = "variable " + temp_var + " equal ";
	      variable_cmd += "-v_fork_force*v_ud" + dims[i]
		+ "_com_" + mother;
	      std::cout << variable_cmd << std::endl;
	      command(variable_cmd);
	    }


	  // apply force to left daughter and descendants
	  fix_cmd = "fix partition_" + ld + " " + ld + " addforce";

	  for (size_t i=0; i<dims.size(); i++)
	    {
	      temp_var = "${lf" + dims[i] + "_" + mother + "}";
	      fix_cmd += " " + temp_var;
	    }

	  fix_cmd += " region sphere_" + mother;

	  std::cout << fix_cmd << std::endl;
	  command(fix_cmd);

	  // apply force to right daughter and descendants
	  fix_cmd = "fix partition_" + rd + " " + rd + " addforce";

	  for (size_t i=0; i<dims.size(); i++)
	    {
	      temp_var = "${rf" + dims[i] + "_" + mother + "}";
	      fix_cmd += " " + temp_var;
	    }

	  fix_cmd += " region sphere_" + mother;

	  std::cout << fix_cmd << std::endl;
	  command(fix_cmd);
	}
    }
  else
    {
      std::string mother, ld, rd;
      std::string unfix_cmd;
      
      for (fork_partition f_p : f_ps)
	{

	  mother = f_p.fork;
	  ld = mother + "l";
	  rd = mother + "r";

	  unfix_cmd = "unfix partition_" + ld;
	  std::cout << unfix_cmd << std::endl;
	  command(unfix_cmd);

	  unfix_cmd = "unfix partition_" + rd;
	  std::cout << unfix_cmd << std::endl;
	  command(unfix_cmd);
	  
	}
    }
}

void LAMMPS_simulator::set_boundary_radius(double radius) {
  std::string region_cmd, fix_cmd;
  // Convert radius to a string
  std::string radius_str = std::to_string(radius+100);


  command("include ${DNA_model_dir}/lmp.DNA_physical_params");
  command("include ${DNA_model_dir}/potentials/lmp.DNA_pair_params");
  // region_cmd = "region boundary_sphere sphere/kk 0.0 0.0 0.0 " + radius_str;
  region_cmd = "region boundary_box block/kk -1050.0 1050.0 -1050.0 1050.0 -1050.0 1050.0";
  std::cout << region_cmd << std::endl;
  command(region_cmd);
  command ("info system");
  command ("info regions");
  // fix_cmd = "fix boundary_force all wall/region/kk boundary_sphere harmonic ${epsilon_soft} 0.0 100.0";
  fix_cmd = "fix boundary_force all wall/region/kk boundary_box harmonic 10 0.0 100.0";
  std::cout << fix_cmd << std::endl;
  command(fix_cmd);
  command("fix_modify boundary_force energy yes");
}
