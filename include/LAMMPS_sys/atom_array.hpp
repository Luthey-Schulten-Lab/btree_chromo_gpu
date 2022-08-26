#ifndef INCLUDE_ATOM_ARRAY_HPP
#define INCLUDE_ATOM_ARRAY_HPP

#include <memory>
#include <fstream>

#include <LAMMPS_sys/vec_quat_manipulator.hpp>

using namespace std;

struct atom
{
  int id, type, mol_id, ellipsoid_flag;
  double density;
  vec r;
};

class atom_array
{
public:

  // constructor and destructor
  atom_array();
  ~atom_array();

  // array sizing
  void set_N(int N);
  int get_N();

  // manipulate columns
  void reset_ids();
  void set_types(int *type);
  void set_type_all(int type);
  void set_ellipsoid_flag_all(int ellipsoid_flat);
  void set_mol_ids(int *mol_id);
  void set_mol_id_all(int mol_id);
  void set_densities(double *density);
  void set_density_all(double density);
  

  // setter and getter for elements
  void set_atom(int i, atom a);
  atom get_atom(int i);

  // write atom data
  void write(fstream &data_file);

private:

  void initialize_atoms();
  void destroy_atoms();

  int N;
  atom *atoms;
  
};
#endif
