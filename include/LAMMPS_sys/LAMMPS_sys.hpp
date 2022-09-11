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

struct sys_bbox
{
  vec r_min, r_max;
};

struct BD_lengths
{
  double r_sphere, r_bdry;
  vec mono_shape, ribo_shape;
};

class LAMMPS_sys
{
public:

  // constructor and destructor
  LAMMPS_sys();
  ~LAMMPS_sys();

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

  vec_quat_manipulator vqm;

};

#endif
