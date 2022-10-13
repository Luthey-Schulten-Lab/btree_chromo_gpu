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


// set the PRNG
void binding_region::set_rand_eng(mt19937 &rand_eng)
{
  this->rand_eng = rand_eng;
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


// select a random monomer within the region to serve as the anchor
int binding_region::select_random_monomer()
{
  // randomly sample a point within the entire region
  uniform_int_distribution<int> unif_dist(0,size-1);
  int r = unif_dist(rand_eng);
  
  if (completed == true)
    {
      r = r + ll;
    }
  else
    {
      if (ter_crossing == true)
	{
	  if (r <= (ll - mid_ll + 1))
	    {
	      r = r + mid_ll;
	    }
	  else
	    {
	      r = (r - (ll - mid_ll + 1)) + ul;
	    }
	}
      else
	{
	  r = r + ll;
	}
    }
  return r;
}


// based on the position of the anchor and a minimum distance between the anchor and hinge, select a direction for the hinge to travel
int binding_region::select_direction(int pos, int min_dist)
{
  // randomly sample a point within the entire region
  uniform_real_distribution<> unif_dist(0.0, 1.0);
  int d;
  
  if (completed == true)
    {
      d = 1;
      if (unif_dist(rand_eng) < 0.5) d = -1;
    }
  else
    {
      if (ter_crossing == true)
	{
	  if ((pos + min_dist) > ll) // no space for hinge at lower-limit
	    {
	      d = -1;
	    }
	  else if ((pos - min_dist) < ul) // no space for hinge at upper-limit
	    {
	      d = 1;
	    }
	  else // space for hinge in either direction
	    {
	      d = 1;
	      if (unif_dist(rand_eng) < 0.5) d = -1;
	    }
	}
      else
	{
	  if ((pos - min_dist) < ll) // no space for hinge at lower-limit
	    {
	      d = 1;
	    }
	  else if ((pos + min_dist) > ul) // no space for hinge at upper-limit
	    {
	      d = -1;
	    }
	  else // space for hinge in either direction
	    {
	      d = 1;
	      if (unif_dist(rand_eng) < 0.5) d = -1;
	    }
	}
    }
  return d;
}


// determine a new relative pos based on the current position, distance, and direction
int binding_region::get_relative_monomer_pos(int pos, int dist, int dir)
{
  int r = pos + dir*dist;
  
  if (completed == true)
    {
      if (r < ll)
	{
	  r = ul - (ll - r - 1);
	}
      else if (r > ul)
	{
	  r = ll + (r - ul - 1);
	}
    }
  else
    {
      if (ter_crossing == true)
	{
	  if ((pos < ll) && (r > ll))
	    {
	      r = ll;
	    }
	  else if ((pos > ul) && (r < ul))
	    {
	      r = ul;
	    }
	  else if (r < mid_ll)
	    {
	      r = mid_ul - (mid_ll - r - 1);
	    }
	  else if (r > mid_ul)
	    {
	      r = mid_ll + (r - mid_ul - 1);
	    }

	}
      else
	{
	  if (r < ll)
	    {
	      r = ll;
	    }
	  else if (r > ul)
	    {
	      r = ul;
	    }
	}
    }
  return r;
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
