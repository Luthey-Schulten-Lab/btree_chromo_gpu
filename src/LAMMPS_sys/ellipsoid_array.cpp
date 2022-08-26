#include <LAMMPS_sys/ellipsoid_array.hpp>

// constructor
ellipsoid_array::ellipsoid_array()
{
  N = -1;
  ellipsoids = nullptr;
}

// destructor
ellipsoid_array::~ellipsoid_array()
{
  destroy_ellipsoids();
}

// set the number of elements
void ellipsoid_array::set_N(int N)
{
  this->N = N;
}

// get the number of elements
int ellipsoid_array::get_N()
{
  return N;
}

// initialize the array
void ellipsoid_array::initialize_ellipsoids()
{
  destroy_ellipsoids();

  if (N > 0)
    {

      ellipsoids = new ellipsoid[N];

      for (int i=0; i<N; i++)
	{
	  // atom id
	  ellipsoids[i].id = i + 1;
	  // shape
	  ellipsoids[i].s.x = -1.0;
	  ellipsoids[i].s.y = -1.0;
	  ellipsoids[i].s.z = -1.0;
	  // orientation
	  ellipsoids[i].q.w = -1.0;
	  ellipsoids[i].q.i = -1.0;
	  ellipsoids[i].q.j = -1.0;
	  ellipsoids[i].q.k = -1.0;
	}
      
    }
  
}

// destroy the array
void ellipsoid_array::destroy_ellipsoids()
{
  if (ellipsoids != nullptr)
    {
      delete[] ellipsoids;
      ellipsoids = nullptr;
    }
}
