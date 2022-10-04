#include <LAMMPS_sys/binding_region.hpp>

// constructor
binding_region::binding_region(string leaf, int ll, int ul, int size,
		 bool completed,
		 bool ter_crossing, int mid_ll, int mid_ul)
{
  this->leaf = leaf;
  this->ll = ll;
  this->ul = ul;
  this->size = size;
  this->completed = completed;
  this->ter_crossing = ter_crossing;
  this->mid_ll = mid_ll;
  this->mid_ul = mid_ul;
}


// destructor
binding_region::~binding_region()
{
}


// getter for size
int binding_region::get_size()
{
  return size;
}


// getter for leaf
string binding_region::get_leaf()
{
  return leaf;
}


// test if queried index is within region
bool binding_region::within_region(int q)
{

  bool within = false;
  
  if (q > 0)
    {
      within = true;
    }

  return within;
}
