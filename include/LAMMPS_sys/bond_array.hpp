#ifndef INCLUDE_BOND_ARRAY_HPP
#define INCLUDE_BOND_ARRAY_HPP

#include <memory>
#include <fstream>

#include <LAMMPS_sys/vec_quat_manipulator.hpp>

using namespace std;

struct bond
{
  int id, type, i, j;
};

class bond_array
{
public:

  // constructor and destructor
  bond_array();
  ~bond_array();

  // array sizing
  void set_N(int N);
  int get_N();

  // setter for elements
  void set_bond(int i, bond b);

  // write bond data
  void write(fstream &data_file);

private:

  void initialize_bonds();
  void destroy_bonds();

  int N;
  bond *bonds;
  
};
#endif
