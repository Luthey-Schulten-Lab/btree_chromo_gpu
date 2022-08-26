#include <LAMMPS_sys/angle_array.hpp>

// constructor
angle_array::angle_array()
{
  N = -1;
  angles = nullptr;
}

// destructor
angle_array::~angle_array()
{
  destroy_angles();
}

// set the number of elements
void angle_array::set_N(int N)
{
  this->N = N;
}

// get the number of elements
int angle_array::get_N()
{
  return N;
}

// initialize the array
void angle_array::initialize_angles()
{
  destroy_angles();

  if (N > 0)
    {

      angles = new angle[N];

      for (int i=0; i<N; i++)
	{
	  // angle id
	  angles[i].id = i + 1;
	  // type
	  angles[i].type = -1;
	  // angle atom ids
	  angles[i].i = -1;
	  angles[i].j = -1;
	  angles[i].k = -1;
	}
      
    }
  
}

// destroy the array
void angle_array::destroy_angles()
{
  if (angles != nullptr)
    {
      delete[] angles;
      angles = nullptr;
    }
}

