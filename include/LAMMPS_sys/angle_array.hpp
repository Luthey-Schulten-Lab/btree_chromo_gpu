#ifndef INCLUDE_ANGLE_ARRAY_HPP
#define INCLUDE_ANGLE_ARRAY_HPP

#include <memory>
#include <fstream>

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

  // array sizing
  void set_N(int N);
  int get_N();

  // setter for elements
  void set_angle(int i, angle a);

  // write angle data
  void write(fstream &data_file);

private:

  void initialize_angles();
  void destroy_angles();

  int N;
  angle *angles;
  
};
#endif
