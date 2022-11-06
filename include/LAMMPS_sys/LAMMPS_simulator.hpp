#ifndef INCLUDE_LAMMPS_SIMULATOR_HPP
#define INCLUDE_LAMMPS_SIMULATOR_HPP

// standard library include files
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <memory>

// OpenMPI include files
#include <mpi.h>

// LAMMPS include files
#include <lammps.h>
#include "domain.h"
#include "group.h"
#include "input.h"
#include "atom.h"
#include "library.h"

// btree_chromo include files
#include <LAMMPS_sys/LAMMPS_sys.hpp>

using namespace std;

struct thermo_dump_parameters
{
  // bool multi_rep; // multiple replicates requiring replicate number
  bool append, write_first; // append to dump file, write first timestep
  // int rep, rep_padding; // replicate number, padding for replicate label
  int dump_freq, thermo_freq; // dump frequency and thermo frequency
};

struct compute_tracker
{
  bool quats, ids, types, MSD;
};

struct dump_tracker
{
  bool lammpstrj;
};

class LAMMPS_simulator
{
public:

  // constructor and destructor
  LAMMPS_simulator();
  ~LAMMPS_simulator();

  void LAMMPS_initialize(string logfile);
  void LAMMPS_destroy();

  void set_lmp_sys(LAMMPS_sys *lmp_sys);

  // sync the simulator to the system
  void sim_to_sys();
  void sim_to_sys_atom_counts();

  // set nProc for simulator
  void set_nProc(int nProc);
  // set the DNA model
  void set_DNA_model_dir(string DNA_model_dir);
  // set output location
  void set_output_details(string output_dir, string output_file_label);
  // set PRNG seed
  void set_prng_seed(int s);
  // set delta_t
  void set_delta_t(double delta_t);

  // simulation protocol
  void reset_protocol_variables();
  void global_setup();
  void standard_computes();
  void reset_Nt(unsigned long Nt);
  void command(string command);
  void include_file(string filename);
  void clear();
  void store_Nt();
  void restore_Nt();

  // read_data
  void read_data(string data_file);

  // minimization routines
  void minimize_soft_harmonic(thermo_dump_parameters t_d_p);
  void minimize_hard_harmonic(thermo_dump_parameters t_d_p);
  void minimize_topoDNA_harmonic(thermo_dump_parameters t_d_p);
  void minimize_soft_FENE(thermo_dump_parameters t_d_p);
  void minimize_hard_FENE(thermo_dump_parameters t_d_p);
  void minimize_topoDNA_FENE(thermo_dump_parameters t_d_p);

  // run routines
  void run_soft_harmonic(unsigned long N_steps, thermo_dump_parameters t_d_p);
  void run_hard_harmonic(unsigned long N_steps, thermo_dump_parameters t_d_p);
  void run_topoDNA_harmonic(unsigned long N_steps, thermo_dump_parameters t_d_p);
  void run_soft_FENE(unsigned long N_steps, thermo_dump_parameters t_d_p);
  void run_hard_FENE(unsigned long N_steps, thermo_dump_parameters t_d_p);
  void run_topoDNA_FENE(unsigned long N_steps, thermo_dump_parameters t_d_p);

  // loop system
  int read_loop_params(string loop_param_filename);
  void set_loop_sim_params(loop_sim_params &l_sim_p);
  void update_loop_bonds(bool new_bonds);
  void run_loops(int N_loops, unsigned long N_steps, thermo_dump_parameters t_d_p);
  
private:

  void compute_trigger(string compute_label);
  void prepare_dump(thermo_dump_parameters &t_d_p);
  void setup_run(thermo_dump_parameters &t_d_p);
  void setup_minimize(thermo_dump_parameters &t_d_p);
  void reset_timestep_to_Nt();
  void set_T_freq(int T_freq);
  void set_D_freq(int D_freq);

  bool T_freq_specified, D_freq_specified;

  int sim_MPI_initialized, sim_MPI_finalized;
  int sim_MPI_size; // MPI size
  int sim_MPI_rank; // current MPI rank

  unsigned long stored_Nt;
  unsigned long Nt;
  int nProc; // number of processors for OpenMP
  int prng_seed; // seed for PRNG within LAMMPS object
  string DNA_model_dir, output_dir, output_file_label; // DNA model, output dir, and label for output files
  double delta_t; // timestep size

  // objects

  loop_sim_params l_sim_p;
  compute_tracker computes_active;
  dump_tracker dumps_active;
  
  LAMMPS_NS::LAMMPS *lmp;
  LAMMPS_sys *lmp_sys;
  
  

};

#endif
