#include <LAMMPS_simulator.hpp>

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
  
  // finalize MPI instance
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
    "-pk", "kokkos"          // Kokkos package
  };
  int lmpargc = sizeof(lmpargv)/sizeof(const char *);

  // Initialize LAMMPS with these arguments
  lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv, MPI_COMM_WORLD);

}


// destroy LAMMPS object
void LAMMPS_simulator::LAMMPS_destroy()
{
  if (lmp != nullptr)
    {
      delete lmp;
      lmp = nullptr;
    }
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

  lmp->input->one(("read_data " + data_file + " extra/bond/per/atom 4").c_str());


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


// clear the LAMMPS system state
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
void LAMMPS_simulator::minimize_topoDNA_harmonic(thermo_dump_parameters t_d_p)
{

  std::cout << "---[ minimizing topoDNA_HARMONIC ]---" << std::endl;

  // setup the minimization
  setup_minimize(t_d_p);

  // include minimization subroutine
  lmp->input->one("include ${DNA_model_dir}/minimize_subroutines/subroutine.min_topoDNA_harmonic");

  // cleanup the minimization
  cleanup_minimize();
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


// run with soft potentials and harmonic bonds
void LAMMPS_simulator::run_soft_harmonic(unsigned long N_steps, thermo_dump_parameters t_d_p)
{
  std::cout << "---[ running SOFT_HARMONIC ]---" << std::endl;
  // setup for run
  setup_run(N_steps,t_d_p);

  // include run subroutine
  lmp->input->one("include ${DNA_model_dir}/run_subroutines/subroutine.run_soft_harmonic");

  // set the timestep
  lmp->input->one("timestep ${delta_t}");

  // run for N_steps
  lmp->input->one(("run " + std::to_string(N_steps)).c_str());

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


// read the replication model
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


// form loops and/or extruded beads bonds
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
				bond_command += " special yes";
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

  // delete the existing loop bonds
  if (new_bonds == false)
    {
      lmp->input->one("delete_bonds DNA bond 2 remove");
    }

  // type array for scatter
  int N = lmp_sys->get_N_total();
  int *types = new int[N];

  // 1 for per-atom type, 1 for size of data (t)
  // lammps_gather_atoms(lmp, const_cast<char*>("type"), 1, 1, types);

  lmp_sys->get_types(types);

  // get the loop bonds
  // std::vector<bond> loop_bonds = lmp_sys->get_loop_bonds();
  std::vector<bond> loop_bonds = sim_ls.get_loop_bonds();
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

  delete[] types;
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
