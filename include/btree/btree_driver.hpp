#ifndef INCLUDE_BTREE_DRIVER_HPP
#define INCLUDE_BTREE_DRIVER_HPP

#include <btree/btree.hpp>
#include <btree/mapper.hpp>
#include <rep_kinetics/replicator.hpp>
#include <LAMMPS_sys/LAMMPS_sys.hpp>
#include <LAMMPS_sys/LAMMPS_simulator.hpp>

using namespace std;

struct drctv_reqs
{
  bool btree_initialized;
  bool topo_update;
  bool CG_update;
  bool regions_present;
  bool rep_model_present;
  bool BD_lengths_present;
  bool map_initial_present;
  bool map_final_present;
  bool simulator_prepared;
  bool DNA_model;
  bool output_details;
  bool delta_t;
  bool lmp_data_present;
  bool loop_params_present;
};


class btree_driver //: protected btree
{
public:

  // constructor and destructor
  btree_driver();
  ~btree_driver();

  // read_directives
  void read_directives(string drctvs_filename);

  // print_directives
  void print_directives();

  // execute directives
  int execute_directives();

private:

  ////////////////
  // directives // 
  ////////////////

  // create a new chromosome
  int new_chromo(vector<string> &params, drctv_reqs &reqs);
  
  // input-ouput
  int input_state(vector<string> &params, drctv_reqs &reqs);
  int output_state(vector<string> &params);

  // transforms
  int transforms_file(vector<string> &params, drctv_reqs &reqs);
  int transform(vector<string> &params, drctv_reqs &reqs);
  int random_transforms(vector<string> &params, drctv_reqs &reqs);

  // regions
  int regions_file(vector<string> &params, drctv_reqs &reqs);
  int dump_regions(vector<string> &params, drctv_reqs &reqs);

  // topology
  int dump_topology(vector<string> &params, drctv_reqs &reqs);
  int update_topology(drctv_reqs &reqs);

  // coarse-graining
  int update_CG_map(vector<string> &params, drctv_reqs &reqs);
  int dump_CG_map(vector<string> &params, drctv_reqs &reqs);

  // miscellaneous
  int btree_prng_seed(vector<string> &params);
  int replicator_prng_seed(vector<string> &params);
  int print_state(drctv_reqs &reqs);

  // replication model
  int load_rep_model(vector<string> &params, drctv_reqs &reqs);
  int replicate(vector<string> &params, drctv_reqs &reqs);

  // LAMMPS system
  // loading coordinates and quaternions
  int load_mono_coords(vector<string> &params, drctv_reqs &reqs);
  int load_mono_quats(vector<string> &params, drctv_reqs &reqs);
  int load_ribo_coords(vector<string> &params, drctv_reqs &reqs);
  int load_ribo_quats(vector<string> &params, drctv_reqs &reqs);
  int load_bdry_coords(vector<string> &params, drctv_reqs &reqs);
  // writing coordinates and quaternions
  int write_mono_coords(vector<string> &params);
  int write_mono_quats(vector<string> &params);
  int write_ribo_coords(vector<string> &params);
  int write_ribo_quats(vector<string> &params);
  int write_bdry_coords(vector<string> &params);
  // loading BD lengths
  int load_BD_lengths(vector<string> &params, drctv_reqs &reqs);
  // writing LAMMPS data file
  int write_LAMMPS_data(vector<string> &params, drctv_reqs &reqs);
  int write_mono_xyz(vector<string> &params, drctv_reqs &reqs);

  // mapper
  int set_initial_state(drctv_reqs &reqs);
  int set_final_state(drctv_reqs &reqs);
  int map_replication(drctv_reqs &reqs);

  // simulator
  int prepare_simulator(vector<string> &params, drctv_reqs &reqs);
  int simulator_include_file(vector<string> &params, drctv_reqs &reqs);
  int sync_simulator_and_system(drctv_reqs &reqs);
  int clear_simulator(drctv_reqs &reqs);
  int simulator_read_data(vector<string> &params, drctv_reqs &reqs);
  int simulator_set_nProc(vector<string> &params, drctv_reqs &reqs);
  int simulator_set_prng_seed(vector<string> &params, drctv_reqs &reqs);
  int simulator_set_DNA_model(vector<string> &params, drctv_reqs &reqs);
  int simulator_set_output_details(vector<string> &params, drctv_reqs &reqs);
  int simulator_set_delta_t(vector<string> &params, drctv_reqs &reqs);
  int simulator_store_timestep(drctv_reqs &reqs);
  int simulator_restore_timestep(drctv_reqs &reqs);
  
  // simulator minimization routines
  template <int SOFT_HARD, int HARMONIC_FENE>
  int simulator_minimize(vector<string> &params, drctv_reqs &reqs);
  
  // simulator run routines
  template <int SOFT_HARD, int HARMONIC_FENE>
  int simulator_run(vector<string> &params, drctv_reqs &reqs);

  // simulator loop routines
  int simulator_load_loop_params(vector<string> &params, drctv_reqs &reqs);
  int simulator_run_loops(vector<string> &params, drctv_reqs &reqs);
  

  /////////////
  // objects //
  /////////////

  // classes
  replicator driver_replicator;
  btree driver_bt;
  mapper driver_mapper;
  LAMMPS_sys driver_lmp_sys;
  LAMMPS_simulator driver_lmp_simulator;

  // structs
  btree_state driver_st; // state structure
  btree_transforms driver_tr; // transform structure
  vector<chromo_region> driver_rg; // vector of chromo_regions
  vector<string> drctvs; // set of directives
  CG_map driver_CG; // coarse-graining map

};

#endif
