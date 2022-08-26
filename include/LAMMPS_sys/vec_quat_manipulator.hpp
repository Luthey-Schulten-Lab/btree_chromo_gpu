#ifndef INCLUDE_VEC_QUAT_MANIPULATOR
#define INCLUDE_VEC_QUAT_MANIPULATOR

using namespace std;

struct vec
{
  double x, y, z;
};

struct quat
{
  double w, i, j, k;
};

class vec_quat_manipulator
{
public:

  // constructor and destructor
  vec_quat_manipulator();
  ~vec_quat_manipulator();
  
private:
  
};

#endif
