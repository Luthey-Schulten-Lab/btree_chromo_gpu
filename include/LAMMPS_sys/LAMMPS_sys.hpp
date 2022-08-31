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

class LAMMPS_sys
{
public:

  // constructor and destructor
  LAMMPS_sys();
  ~LAMMPS_sys();

  // set the btree
  void set_btree(btree_state in_state);

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

  btree internal_btree;
  boundary_surface b_surf;

};

#endif
