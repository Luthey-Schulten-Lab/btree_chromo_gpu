#ifndef INCLUDE_BTREE_DRIVER_HPP
#define INCLUDE_BTREE_DRIVER_HPP

#include <unordered_map>

#include <btree/btree.hpp>
#include <btree/mapper.hpp>
#include <rep_kinetics/replicator.hpp>
#include <LAMMPS_sys/LAMMPS_sys.hpp>
#include <LAMMPS_sys/LAMMPS_simulator.hpp>

using namespace std;

struct lock
{
  string key;
  bool s;
};


class btree_driver //: protected btree
{
public:

  // constructor and destructor
  btree_driver();
  ~btree_driver();

  // read directives from a file
  void read_directives(string drctvs_filename);

  // print the directives
  void print_directives();

  // parse the directives
  void parse_directives();

  // validate the parameters in the command sequence
  int validate_command_sequence_parameters();

  // validate the command sequence
  int validate_command_sequence();

  // expand the metacommands
  int expand_metacommands();

  // print the command sequence
  void print_commands();
  
  // execute commands
  int execute_commands();

private:

  // parse a single directive
  void parse_single_directive(string drctv, string &command, vector<string> &params);

  // test if metacommands are properly paired
  int test_paired_metacommands(string paired_command);
  // compose any metacommands, loops etc.
  void compose_metacommands();
  // expand repeat metacommands
  void expand_repeat_metacommands();
  void expand_repeat_replicates_metacommands();
  void update_replicate_modified_params(string &rep_mod, string &command, vector<string> &params);
  void append_replicate_modifier(string &rep_mod, string &mod_param);
  void insert_replicate_modifier(string &rep_mod, string &mod_param);
  
  // execute a single command
  int execute_single_command(string &command, vector<string> &params);

  // reset the command locks and updates
  void reset_command_locks_and_updates();
  // set the command requirements
  void prepare_command_requirements();
  lock new_lock(string key, bool s);

  // test the command for parameter validity
  int test_command_parameter_validity(string &command, vector<string> &params);
  // test the command for lock state validity
  int test_command_lock_validity(string &command);
  // update the lock state following a command
  void update_lock_state_post_command(string &command);
  

  //////////////
  // commands // 
  //////////////
  

  // create a new chromosome
  int new_chromo(vector<string> &params);
  
  // input-ouput
  int input_state(vector<string> &params);
  int output_state(vector<string> &params);

  // transforms
  int transforms_file(vector<string> &params);
  int transform(vector<string> &params);
  int random_transforms(vector<string> &params);

  // regions
  int regions_file(vector<string> &params);
  int dump_regions(vector<string> &params);

  // topology
  int dump_topology(vector<string> &params);
  int update_topology();

  // coarse-graining
  int update_CG_map(vector<string> &params);
  int dump_CG_map(vector<string> &params);

  // miscellaneous
  int btree_prng_seed(vector<string> &params);
  int replicator_prng_seed(vector<string> &params);
  int print_state();

  // replication model
  int load_rep_model(vector<string> &params);
  int replicate(vector<string> &params);

  // LAMMPS system
  // loading coordinates and quaternions
  int load_mono_coords(vector<string> &params);
  int load_mono_quats(vector<string> &params);
  int load_ribo_coords(vector<string> &params);
  int load_ribo_quats(vector<string> &params);
  int load_bdry_coords(vector<string> &params);
  // writing coordinates and quaternions
  int write_mono_coords(vector<string> &params);
  int write_mono_quats(vector<string> &params);
  int write_ribo_coords(vector<string> &params);
  int write_ribo_quats(vector<string> &params);
  int write_bdry_coords(vector<string> &params);
  // manual boundary specification
  int spherical_bdry(vector<string> &params);
  // loading BD lengths
  int load_BD_lengths(vector<string> &params);
  // manipulate system interactions
  int switch_bonds(vector<string> &params);
  int switch_bending_angles(vector<string> &params);
  int switch_twisting_angles(vector<string> &params);
  // writing LAMMPS data file
  int write_LAMMPS_data(vector<string> &params);
  int write_mono_xyz(vector<string> &params);

  // mapper
  int set_initial_state();
  int set_final_state();
  int map_replication();

  // simulator
  int prepare_simulator(vector<string> &params);
  int simulator_include_file(vector<string> &params);
  int sync_simulator_and_system();
  int clear_simulator();
  int simulator_read_data(vector<string> &params);
  int simulator_set_nProc(vector<string> &params);
  int simulator_set_prng_seed(vector<string> &params);
  int simulator_set_DNA_model(vector<string> &params);
  int simulator_set_output_details(vector<string> &params);
  int simulator_set_delta_t(vector<string> &params);
  int simulator_store_timestep();
  int simulator_restore_timestep();
  int simulator_reset_prev_dump_timestep(vector<string> &params);
  int simulator_reset_timestep(vector<string> &params);
  
  // simulator minimization routines
  template <int SOFT_HARD, int HARMONIC_FENE>
  int simulator_minimize(vector<string> &params);
  
  // simulator run routines
  template <int SOFT_HARD, int HARMONIC_FENE>
  int simulator_run(vector<string> &params);

  // simulator looped DNA routines
  int simulator_load_loop_params(vector<string> &params);
  int simulator_run_loops(vector<string> &params);
  

  /////////////
  // objects //
  /////////////

  // internal classes for directive execution
  replicator driver_replicator;
  btree driver_bt;
  mapper driver_mapper;
  LAMMPS_sys driver_lmp_sys;
  LAMMPS_simulator driver_lmp_simulator;

  // internal variables for directive execution
  btree_state driver_st; // state structure
  btree_transforms driver_tr; // transform structure
  vector<chromo_region> driver_rg; // vector of chromo_regions
  CG_map driver_CG; // coarse-graining map

  // directives, commands, and parameters
  vector<string> drctvs; // set of directives
  vector<string> commands; // vector of commands as strings
  vector<vector<string>> command_params; // vector command parameters as vectors of strings

  // variables to hold requirements for directives
  unordered_map<string,bool> lock_state; // locks for the command sequence
  unordered_map<string,vector<lock>> lock_tests; // lock requirements for command execution
  unordered_map<string,vector<lock>> lock_updates; // lock updates given successful command execution
  unordered_map<string,size_t> N_param_reqs; // requirements for the number of parameters

};

#endif
