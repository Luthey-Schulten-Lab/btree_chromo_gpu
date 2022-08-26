#include <LAMMPS_sys/bond_array.hpp>

// constructor
bond_array::bond_array()
{
  N = -1;
  bonds = nullptr;
}

// destructor
bond_array::~bond_array()
{
  destroy_bonds();
}

// set the number of elements
void bond_array::set_N(int N)
{
  this->N = N;
}

// get the number of elements
int bond_array::get_N()
{
  return N;
}

// initialize the array
void bond_array::initialize_bonds()
{
  destroy_bonds();

  if (N > 0)
    {

      bonds = new bond[N];

      for (int i=0; i<N; i++)
	{
	  // bond id
	  bonds[i].id = i + 1;
	  // type
	  bonds[i].type = -1;
	  // bonded atom ids
	  bonds[i].i = -1;
	  bonds[i].j = -1;
	}
      
    }
  
}

// destroy the array
void bond_array::destroy_bonds()
{
  if (bonds != nullptr)
    {
      delete[] bonds;
      bonds = nullptr;
    }
}

