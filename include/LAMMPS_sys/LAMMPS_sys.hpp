#ifndef INCLUDE_LAMMPS_SYS_HPP
#define INCLUDE_LAMMPS_SYS_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <memory>

#include <btree/btree.hpp>
#include <LAMMPS_sys/vec_quat_manipulator.hpp>
#include <LAMMPS_sys/boundary_surface.hpp>
#include <LAMMPS_sys/atom_array.hpp>
#include <LAMMPS_sys/ellipsoid_array.hpp>
#include <LAMMPS_sys/bond_array.hpp>
#include <LAMMPS_sys/angle_array.hpp>

using namespace std;

// structure defining bounding box for the system
struct sys_bbox
{
  vec r_min, r_max;
};

// structure containing lengths for spatial system
struct BD_lengths
{
  double r_sphere, r_bdry;
  vec mono_shape, ribo_shape;
};

// structure defining a binding region for a loop
struct binding_region
{
  string leaf; // leaf that binding region belongs to
  int ll, ul, size; // lower limit and upper limit of indices
  bool completed;
  bool ter_crossing;
  int mid_ll, mid_ul;
};

// structure defining a loop
struct loop
{
  int a, h; // indices of anchor and hinge
  int d; // direction of loop extrusion
  binding_region a_region; // region containing anchor
  binding_region h_region; // region containing hinge
};

class LAMMPS_sys
{
public:

  // constructor and destructor
  LAMMPS_sys();
  ~LAMMPS_sys();

  // prng seeding
  void prng_seed(int s);

  // set the btree
  void set_btree(btree_state in_state);

  // read the Brownian dynamics lengths
  void read_BD_lengths(string lengths_filename);

  // read the monomer coordinates
  int read_mono_coords(string coords_filename, string order);
  // read the ribo coordinates
  int read_ribo_coords(string coords_filename, string order);
  // read the bdry coordinates
  int read_bdry_coords(string coords_filename, string order);

  // apply a mapping to the monomers
  void apply_mono_mapping(vector<vector<array<int,3>>> map);

  // read the monomer coordinates
  void write_mono_xyz(string data_filename);

  // write the data
  void write_data(string data_filename);

  // getters for the total system state
  int get_N_total();
  int get_N_mono();
  int get_N_ribo();
  int get_N_bdry();

  // setters for the total system state
  void set_coords_arr_total(double *&x, string order);
  void set_quats_arr_total(double *&q, string order);

  // sync the subarrays
  void sync_subarrays();

  // loop topology functions
  void prepare_binding_regions();
  void initialize_loop_topo(int N_loop);
  void update_loop_topo();
  vector<loop> get_loops();
  
private:

  void prepare_test_data();

  // concatenate arrays
  void cat_atom_array(atom_array &in_atoms, atom_array &cat_atoms);
  void cat_ellipsoid_array(ellipsoid_array &in_ellipsoids, ellipsoid_array &cat_ellipsoids);

  // determine the system's topology
  void prepare_topology();
  void set_bonds();
  void set_angles();
  void set_mono_types(int base_type);

  // prepare the system prior to output
  void finalize_system();
  void merge_system_components();
  void calc_bbox();

  int N_atom_types;
  int N_angle_types, N_bond_types;

  atom_array atoms, mono_atoms, ribo_atoms, bdry_atoms;
  ellipsoid_array ellipsoids, mono_ellipsoids, ribo_ellipsoids;
  bond_array bonds;
  angle_array angles;

  sys_bbox bbox;
  BD_lengths BD_l;
  btree internal_btree;
  boundary_surface b_surf;

  // loop topology
  vector<loop> loops;
  vector<binding_region> regions;

  vec_quat_manipulator vqm;
  mt19937 rand_eng;

};

#endif
