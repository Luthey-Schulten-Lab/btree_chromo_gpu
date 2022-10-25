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

  proximities = new int[size];

  reset_proximities();
  
  mono_idx = new int[size];

  if (ter_crossing == true)
    {
      int i = 0;
      for (int j=ul; j<(mid_ul+1); j++)
	{
	  mono_idx[i] = j;
	  reg_idx[j] = i;
	  i += 1;
	}
      for (int j=mid_ll; j<(mid_ll+1); j++)
	{
	  mono_idx[i] = j;
	  reg_idx[j] = i;
	  i += 1;
	}
    }
  else
    {
      for (int i=0; i<size; i++)
	{
	  mono_idx[i] = ll + i;
	  reg_idx[ll+i] = i;
	}
    }
  
}


// destructor
binding_region::~binding_region()
{
  if (mono_idx != nullptr)
    {
      delete[] mono_idx;
      mono_idx = nullptr;
    }
  if (proximities != nullptr)
    {
      delete[] proximities;
      proximities = nullptr;
    }
  reg_idx.clear();
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


// getter for region position from monomer position
int binding_region::get_reg_pos(int mono_pos)
{
  return reg_idx[mono_pos];
}


// getter for monomer position from region position
int binding_region::get_mono_pos(int reg_pos)
{
  return mono_idx[reg_pos];
}


// based on the position of the anchor and a minimum distance between the anchor and hinge, select a direction for the hinge to travel
int binding_region::select_direction(int mono_pos, int min_dist, double r_d)
{
  // randomly sample a point within the entire region
  int d;
  
  if (completed == true)
    {
      d = 1;
      if (r_d < 0.5) d = -1;
    }
  else
    {
      int reg_pos = get_reg_pos(mono_pos);

      if (reg_pos < min_dist)
	{
	  d = 1;
	}
      else if (reg_pos >= (size - min_dist))
	{
	  d = -1;
	}
      else
	{
	  d = 1;
	  if (r_d < 0.5) d = -1;
	}
    }
  return d;
}


// determine a new relative pos based on the current position, distance, and direction
int binding_region::get_relative_monomer_pos(int mono_pos, int dist, int dir)
{
  int reg_pos = get_reg_pos(mono_pos) + dir*dist;
  
  if (completed == true)
    {
      if (reg_pos < 0)
	{
	  reg_pos = size + reg_pos;
	}
      else if (reg_pos >= size)
	{
	  reg_pos = reg_pos - size;
	}
    }
  else
    {
      if (reg_pos < 0)
	{
	  reg_pos = 0;
	}
      else if (reg_pos >= size)
	{
	  reg_pos = size - 1;
	}
    }
  return get_mono_pos(reg_pos);
}


// get the distance between two monomer positions in a region
int binding_region::get_dist(int mono_pos_i, int mono_pos_j)
{
  int reg_pos_i = get_reg_pos(mono_pos_i);
  int reg_pos_j = get_reg_pos(mono_pos_j);
  int d;

  if (completed == true)
    {
      if (reg_pos_i > reg_pos_j)
	{
	  d = min(reg_pos_i-reg_pos_j,size-reg_pos_i+reg_pos_j);
	}
      else if (reg_pos_j > reg_pos_i)
	{
	  d = min(reg_pos_j-reg_pos_i,size-reg_pos_j+reg_pos_i);
	}
      else
	{
	  d = 0;
	}
    }
  else
    {
      d = abs(reg_pos_i - reg_pos_j);
    }
  return d;
}


// reset the proximities
void binding_region::reset_proximities()
{
  for (int i=0; i<size; i++)
    {
      proximities[i] = 0;
    }
}


// update the proximities given a grab radius and coordinate
void binding_region::update_proximities(double r_g, vec a_coord, vector<vec> &coords)
{
  double r_g_2 = pow(r_g,2.0);
  double test_dist_L2;
  
  for (int i=0; i<size; i++)
    {
      cout << "i = " << i << endl;
      cout << "m_i = " << get_mono_pos(i) << endl;
      test_dist_L2 = vqm.v_L2(vqm.v_xpy(a_coord,vqm.v_inv(coords[get_mono_pos(i)])));
      if (test_dist_L2 < r_g_2)
	{
	  proximities[i] = 1;
	}
      else
	{
	  proximities[i] = 0;
	}
    }
}


// filter the proximites too close to the anchor
void binding_region::filter_intra_proximities_near_a(int min_dist, int a_mono_pos)
{
  int a_reg_pos = get_reg_pos(a_mono_pos);
  
  for (int i=max(a_reg_pos-min_dist,0); i<min(a_reg_pos+min_dist+1,size); i++)
    {
      proximities[i] = 0;
    }
}


// get the intra candidates and remove them from the set of proximities
vector<int> binding_region::get_and_filter_intra_candidates(int ext_max, int h_mono_pos, int dir)
{
  vector<int> intra_candidates;
  intra_candidates.push_back(h_mono_pos);

  int h_reg_pos = get_reg_pos(h_mono_pos);

  if (dir == 1)
    {
      for (int i=h_reg_pos+1; i<min(h_reg_pos+ext_max+1,size); i++)
	{
	  if (proximities[i] == 1)
	    {
	      intra_candidates.push_back(get_mono_pos(i));
	      proximities[i] = 0;
	    }
	  else
	    {
	      break;
	    }
	}
    }
  else if (dir == -1)
    {
      for (int i=max(h_reg_pos-ext_max,0); i<h_reg_pos; i++)
	{
	  if (proximities[i] == 1)
	    {
	      intra_candidates.push_back(get_mono_pos(i));
	      proximities[i] = 0;
	    }
	  else
	    {
	      break;
	    }
	}
    }
  return intra_candidates;
}

// get the inter candidates
vector<int> binding_region::get_inter_candidates()
{
  vector<int> inter_candidates;

  for (int i=0; i<size; i++)
    {
      if (proximities[i] == 1) inter_candidates.push_back(get_mono_pos(i));
    }
  
  return inter_candidates;
}
