#ifndef INCLUDE_ELLIPSOID_ARRAY_HPP
#define INCLUDE_ELLIPSOID_ARRAY_HPP

#include <memory>

#include <LAMMPS_sys/vec_quat_manipulator.hpp>

using namespace std;

struct ellipsoid
{
  int id;
  vec s;
  quat q;
};

class ellipsoid_array
{
public:

  // constructor and destructor
  ellipsoid_array();
  ~ellipsoid_array();

  void initialize_ellipsoids();
  void destroy_ellipsoids();

  void set_N(int N);
  int get_N();

private:

  int N;
  ellipsoid *ellipsoids;
  
};
#endif
