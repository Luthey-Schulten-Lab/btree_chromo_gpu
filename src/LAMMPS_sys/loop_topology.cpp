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


void loop_topology::set_coords(vector<vec> coords)
{
  this->coords = coords;
}


// prepare the vector of binding regions
void loop_topology::prepare_binding_regions(vector<string> leaves, vector<theta_topo> leaf_topos, int *&t)
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
  
}


// prepare the vector of binding regions
void loop_topology::initialize_loops(int N_loops, int min_dist)
{
  // clear the loops and binding regions
  loops.clear();
  
  size_t N_regions = regions.size();

  // only select binding regions that can contain a full loop (hinge and anchor) upon initialization

  int total_binding_region_size = 0;

  for (size_t i_region=0; i_region<N_regions; i_region++)
    {
      if (regions[i_region].get_size() > 2*min_dist)
	{
	  total_binding_region_size += regions[i_region].get_size();
	}
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
	  if (regions[i_region].get_size() > 2*min_dist)
	    {
	      accumulator += regions[i_region].get_size();
	    } 
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
      // get the anchor region
      a_region = loops[i_loop].get_a_region();

      // select a random monomer for the anchor
      loops[i_loop].set_a(regions[a_region].select_random_monomer());
      
      // bind the anchor
      loops[i_loop].set_a_bound(true);

      // get the hinge region
      h_region = loops[i_loop].get_h_region();

      // based on the position of the anchor and a minimum distance between the anchor and hinge, select a direction for the hinge to travel
      loops[i_loop].set_d(regions[h_region].select_direction(loops[i_loop].get_a(),
							     min_dist));

      // select a hinge based on the anchor position and direction
      loops[i_loop].set_h(regions[h_region].get_relative_monomer_pos(loops[i_loop].get_a(),
								     min_dist,
								     loops[i_loop].get_d()));
      // bind the hinge
      loops[i_loop].set_h_bound(true);
    }

}


// prepare the vector of binding regions
void loop_topology::update_loops(int ext_avg, int ext_max, int min_dist, double p_unbinding, double r_g)
{

  uniform_real_distribution<> unif_dist(0.0, 1.0);
  vector<int> intra_updates;
  vector<vector<int>> inter_updates;
  int a_mono, a_reg, h_mono, h_reg, d;
  vec a_coord;
  bool h_bound;

  cout << ext_avg << endl;
  
  // loop over the loops
  for (size_t i_loop=0; i_loop<loops.size(); i_loop++ )
    {

      // get the anchor, anchor region, and anchor coordinate
      a_mono = loops[i_loop].get_a();
      a_reg = loops[i_loop].get_a_region();
      a_coord = coords[a_mono];

      // get the hinge, hinge region, and hinge binding state
      h_mono = loops[i_loop].get_h();
      h_reg = loops[i_loop].get_h_region();
      h_bound = loops[i_loop].get_h_bound();

      // hinge is already bound
      if (h_bound == true)
	{

	  // get the direction
	  d = loops[i_loop].get_d();

	  // update the proximities in the hinge region
	  regions[h_reg].update_proximities(r_g,a_coord,coords);

	  // get the intra updates
	  intra_updates = regions[h_reg].get_and_filter_intra_candidates(ext_max,h_mono,d);

	  if (unif_dist(rand_eng) > p_unbinding)
	    {
	  
	      // select an updated hinge for 1D motion along strand
	      
	  
	    }
	  else
	    {
	  
	      // select an updated hinge for 3D motion between strands

	      int N_inter_total = 0;

	      // determine the candidate inter-strand hinge updates for all the regions
	      for (size_t i_reg=0; i_reg<regions.size(); i_reg++)
		{
		  // update proximities not in the hinge region
		  if (i_reg != h_reg)
		    {
		      regions[i_reg].update_proximities(r_g,a_coord,coords);
		    }

		  // filter proximities violating minimum distance within anchor region
		  if (i_reg == a_reg)
		    {
		      regions[i_reg].filter_intra_proximities_near_a(min_dist,a_mono);
		    }

		  // add the region's candidates to the total set of possible inter updates
		  inter_updates.push_back(regions[i_reg].get_inter_candidates());

		  N_inter_total += inter_updates[i_reg].size();
	      
		} // end loop to determine inter candidates

	      if (N_inter_total > 0)
		{
		  uniform_int_distribution<int> unif_dist(1,N_inter_total);
		  int accumulator;

		  int h_inter = unif_dist(rand_eng);

		  // determine the candidate inter-strand hinge updates for all the regions
		  for (size_t i_reg=0; i_reg<regions.size(); i_reg++)
		    {

		      if ((accumulator + inter_updates[i_reg].size()) >= h_inter)
			{
			  h_inter -= accumulator;
			  loops[i_loop].set_h(regions[i_reg].get_mono_pos(h_inter));
			  loops[i_loop].set_h_region(i_reg);
			}
		      else
			{
			  accumulator += inter_updates[i_reg].size();
			}
	      
		    }
		  
		  
		}
	      else
		{

		  // make the hinge unbound
		  loops[i_loop].set_h_bound(false);
		  
		}

	    }
	  
	}
      else // hinge is unbound
	{

	  // select an updated hinge for 3D motion between strands
	  int N_inter_total = 0;

	  // determine the candidate inter-strand hinge updates for all the regions
	  for (size_t i_reg=0; i_reg<regions.size(); i_reg++)
	    {
	      // update proximities for all regions
	      regions[i_reg].update_proximities(r_g,a_coord,coords);

	      // filter proximities violating minimum distance within anchor region
	      if (i_reg == a_reg)
		{
		  regions[i_reg].filter_intra_proximities_near_a(min_dist,a_mono);
		}

	      // add the region's candidates to the total set of possible inter updates
	      inter_updates.push_back(regions[i_reg].get_inter_candidates());

	      N_inter_total += inter_updates[i_reg].size();
	      
	    } // end loop over regions

	  if (N_inter_total > 0)
	    {
	      uniform_int_distribution<int> unif_dist(1,N_inter_total);
	      int accumulator;

	      int h_inter = unif_dist(rand_eng);

	      // determine the candidate inter-strand hinge updates for all the regions
	      for (size_t i_reg=0; i_reg<regions.size(); i_reg++)
		{

		  if ((accumulator + inter_updates[i_reg].size()) >= h_inter)
		    {
		      h_inter -= accumulator;
		      // set the hinge
		      loops[i_loop].set_h(regions[i_reg].get_mono_pos(h_inter));
		      // set the hinge region
		      loops[i_loop].set_h_region(i_reg);
		      // bind the hinge
		      loops[i_loop].set_h_bound(true);
		    }
		  else
		    {
		      accumulator += inter_updates[i_reg].size();
		    }
	      
		}
	      
	    } // end conditional for nonzero inter candidates
	  
	} // end conditional for bound/unbound hinges
      
    } // end loop over loops
  
}


// get the loops
vector<loop> loop_topology::get_loops()
{
  return loops;
}

// get the binding regions
vector<binding_region> loop_topology::get_regions()
{
  return regions;
}

