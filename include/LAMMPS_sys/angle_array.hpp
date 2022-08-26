#ifndef INCLUDE_ANGLE_ARRAY_HPP
#define INCLUDE_ANGLE_ARRAY_HPP

#include <memory>

#include <LAMMPS_sys/vec_quat_manipulator.hpp>

using namespace std;

struct angle
{
  int id, type, i, j, k;
};

class angle_array
{
public:

  // constructor and destructor
  angle_array();
  ~angle_array();

  void initialize_angles();
  void destroy_angles();

  void set_N(int N);
  int get_N();

private:

  int N;
  angle *angles;
  
};
#endif
