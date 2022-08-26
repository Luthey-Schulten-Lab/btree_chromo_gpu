#ifndef INCLUDE_BOND_ARRAY_HPP
#define INCLUDE_BOND_ARRAY_HPP

#include <memory>

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

  void initialize_bonds();
  void destroy_bonds();

  void set_N(int N);
  int get_N();

private:

  int N;
  bond *bonds;
  
};
#endif
