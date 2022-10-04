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

using namespace std;


class loop_topology
{
public:

  // constructor and destructor
  loop_topology();
  ~loop_topology();

  // prng seeding
  void prng_seed(int s);

  // loop topology functions
  void prepare_binding_regions(vector<string> leaves, vector<theta_topo> leaf_topos, int *t);
  void initialize_loops(int N_loops);
  void update_loops();
  vector<loop> get_loops();
  
private:

  // loop topology
  vector<loop> loops;
  vector<binding_region> regions;

  vec_quat_manipulator vqm;
  mt19937 rand_eng;

};

#endif
