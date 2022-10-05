#include <LAMMPS_sys/loop_topology.hpp>

// constructor
loop_topology::loop_topology()
{
  rand_eng.seed(0);
}


// destructor
loop_topology::~loop_topology()
{
}


void loop_topology::prng_seed(int s)
{
  rand_eng.seed(s);
}


// prepare the vector of binding regions
void loop_topology::prepare_binding_regions(vector<string> leaves, vector<theta_topo> leaf_topos, int *t)
{

  regions.clear();

  vector<int> partitions;

  string leaf; // leaf that binding region belongs to
  int ll, ul, size; // lower limit and upper limit of indices
  bool completed;
  bool ter_crossing;
  int mid_ll, mid_ul;

  theta_topo leaf_topo;

  for (size_t i_leaf=0; i_leaf<leaves.size(); i_leaf++)
    {

      leaf = leaves[i_leaf];
      leaf_topo = leaf_topos[i_leaf];

      // test if the leaf is a completed chromosome
      if ((leaf_topo.start == leaf_topo.end_link) &&
	  (leaf_topo.end == leaf_topo.start_link))
	{

	  // determine the partitions along the leaf
	  partitions.clear();

	  for (int i=(leaf_topo.start+1); i<leaf_topo.end; i++)
	    {
	      if (t[i] == 6)
		{
		  partitions.push_back(i);
		}
	    }

	  // prepare a parititioning if forks are present
	  if (partitions.size() > 0)
	    {
	    
	      // create regions from the partitions
	      for (size_t i=0; i<(partitions.size()-1); i++)
		{
		  if ((partitions[i+1] - partitions[i]) > 3)
		    {
		      ll = partitions[i] + 1;
		      ul = partitions[i+1] - 1;
		      size = ul - ll + 1;
		      completed = false;
		      ter_crossing = false;
		      mid_ll = -1;
		      mid_ul = -1;
		      binding_region region(leaf,ll,ul,size,
					    completed,
					    ter_crossing,mid_ll,mid_ul);
		      regions.push_back(region);
		    }
		}

	      // test for crossing back over Ter
	      int completed_gap = 0;

	      completed_gap += (leaf_topo.end - partitions[partitions.size()-1]);
	      completed_gap += (partitions[0] - leaf_topo.start);

	      if (completed_gap > 1)
		{
		  ter_crossing = true;
		  if (partitions[partitions.size()-1] < leaf_topo.end)
		    {
		      ul = partitions[partitions.size()-1] + 1;
		      mid_ul = leaf_topo.end;
		    }
		  else
		    {
		      ter_crossing = false;
		      ll = leaf_topo.start;
		      ul = partitions[0] - 1;
		    }
		  if (partitions[0] > leaf_topo.start)
		    {
		      ll = partitions[0] - 1;
		      mid_ll = leaf_topo.start;
		    }
		  else
		    {
		      ter_crossing = false;
		      ul = leaf_topo.end;
		      ll = partitions[partitions.size()-1] + 1;
		    }
		  completed = false;
		  size = completed_gap;
		  binding_region region(leaf,ll,ul,size,
					completed,
					ter_crossing,mid_ll,mid_ul);
		  regions.push_back(region);
		}
	      
	    }
	  else // no forks are present
	    {

	      ll = leaf_topo.start;
	      ul = leaf_topo.end;
	      size = ul - ll + 1;
	      completed = true;
	      ter_crossing = false;
	      mid_ll = -1;
	      mid_ul = -1;
	      binding_region region(leaf,ll,ul,size,
				    completed,
				    ter_crossing,mid_ll,mid_ul);
	      regions.push_back(region);
	      
	    }
	  
	}
      else
	{

	  // determine the partitions along the leaf
	  partitions.clear();

	  partitions.push_back(leaf_topo.start);
	  
	  for (int i=(leaf_topo.start+1); i<leaf_topo.end; i++)
	    {
	      if (t[i] == 6)
		{
		  partitions.push_back(i);
		}
	    }

	  partitions.push_back(leaf_topo.end);

	  // create regions from the partitions
	  for (size_t i=0; i<(partitions.size()-1); i++)
	    {
	      if ((partitions[i+1] - partitions[i]) > 3)
		{
		  ll = partitions[i] + 1;
		  ul = partitions[i+1] - 1;
		  size = ul - ll + 1;
		  completed = false;
		  ter_crossing = false;
		  mid_ll = -1;
		  mid_ul = -1;
		  binding_region region(leaf,ll,ul,size,
					completed,
					ter_crossing,mid_ll,mid_ul);
		  region.set_rand_eng(rand_eng);
		  regions.push_back(region);
		}
	    }
	  
	} // end conditinal for complete leaf test
      
    } // end loop over leaves

  // delete the type array
  delete[] t;
  
}


// prepare the vector of binding regions
void loop_topology::initialize_loops(int N_loops)
{
  // clear the loops and binding regions
  loops.clear();

  size_t N_regions = regions.size();

  int total_binding_region_size = 0;

  for (size_t i_region=0; i_region<N_regions; i_region++)
    {
      total_binding_region_size += regions[i_region].get_size();
    }

  uniform_int_distribution<int> unif_dist(1,total_binding_region_size);
  vector<int> a_dist;

  // loop over the number of loops
  for (int i_loop=0; i_loop<N_loops; i_loop++)
    {
      a_dist.push_back(unif_dist(rand_eng));
    }

  // select an anchor region
  for (int i_loop=0; i_loop<N_loops; i_loop++)
    {
      int accumulator = 0;
      for (size_t i_region=0; i_region<N_regions; i_region++)
	{
	  accumulator += regions[i_region].get_size();
	  if (a_dist[i_loop] <= accumulator)
	    {
	      loop l(i_region,i_region);
	      loops.push_back(l);
	      break;
	    }
	}
    }


  int a_region, h_region;
  
  // select an anchor for each loop
  for (int i_loop=0; i_loop<N_loops; i_loop++)
    {
      a_region = loops[i_loop].get_a_region();
      loops[i_loop].set_a(regions[a_region].select_random_monomer());
    }

  // based on the position of the anchor and a minimum distance between the anchor and hinge, select a direction for the hinge to travel
  int min_a_d_dist = 1;
  for (int i_loop=0; i_loop<N_loops; i_loop++)
    {
      h_region = loops[i_loop].get_h_region();
      loops[i_loop].set_d(regions[h_region].select_direction(loops[i_loop].get_a(),
							     min_a_d_dist));
    }

  // select a compatible hinge
  for (int i_loop=0; i_loop<N_loops; i_loop++)
    {
      h_region = loops[i_loop].get_h_region();
      loops[i_loop].set_h(regions[h_region].get_relative_monomer_pos(loops[i_loop].get_a(),
								     min_a_d_dist,
								     loops[i_loop].get_d()));
    }
}


// prepare the vector of binding regions
void loop_topology::update_loops()
{

  // loop over the loops

  // select an updated hinge
  
}


// prepare the vector of binding regions
vector<loop> loop_topology::get_loops()
{
  return loops;
}

