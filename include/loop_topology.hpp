#ifndef INCLUDE_LOOP_TOPOLOGY_HPP
#define INCLUDE_LOOP_TOPOLOGY_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <btree/btree.hpp>
#include <LAMMPS_sys/vec_quat_manipulator.hpp>
#include <LAMMPS_sys/loop.hpp>
#include <LAMMPS_sys/binding_region.hpp>
#include <LAMMPS_sys/intra_step_distribution.hpp>

using namespace std;

struct loop_sys_params
{
  int min_dist;
  string family;
  double ext_avg;
  int  ext_max;
  double p_unbinding, r_g;
};

struct loop_sim_params
{
  double r_0, k;
  unsigned long freq_loop, freq_topo, dNt_topo;
};

class loop_topology
{
public:

  // constructor and destructor
  loop_topology();
  ~loop_topology();

  // prng seeding
  void prng_seed(int s);

  // set the step distribution
  void set_step_dist(string family, double l, int k_max);

  // set the coords
  void set_coords(vector<vec> coords);

  // loop topology functions
  void prepare_binding_regions(vector<string> leaves, vector<theta_topo> leaf_topos, int *&t);
  void initialize_loops(int N_loops, int min_dist);
  void update_loops(int ext_max, int min_dist, double p_unbinding, double r_g);
  vector<loop> get_loops();
  vector<binding_region> get_regions();
  
private:

  // loop topology
  vector<loop> loops;
  vector<binding_region> regions;
  vector<vec> coords;

  // step distribution
  intra_step_distribution step_dist;

  // prng
  mt19937 rand_eng;

};

#endif
