#include <loop_topology.hpp>

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


void loop_topology::set_step_dist(std::string family, double l, int k_max)
{
  if (family == "poisson")
    {
      step_dist.poisson_dist(l,k_max+1);
    }
  else if (family == "uniform")
    {
      step_dist.uniform_dist(k_max+1);
    }
}


void loop_topology::set_coords(std::vector<vec> coords)
{
  this->coords = coords;
}


// prepare the vector of binding regions
void loop_topology::prepare_binding_regions(std::vector<std::string> leaves, std::vector<theta_topo> leaf_topos, int *&t)
{
  std::cout << "Preparing binding regions..."<< std::endl;
  std::cout << "regions size, pre-clear: " << regions.size() << std::endl;
  regions.clear();
  std::cout << "regions size, post-clear: " << regions.size() << std::endl;
  std::cout << "Finished clearing binding regions..."<< std::endl;
  std::vector<int> partitions;

  std::string leaf; // leaf that binding region belongs to
  int ll, ul, size; // lower limit and upper limit of indices
  bool completed;
  bool ter_crossing;
  int mid_ll = -1;
  int mid_ul = -1;

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
		  regions.push_back(region);
		}
	    }

	} // end conditinal for complete leaf test

    } // end loop over leaves

  // std::cout << "at end of prepare_binding_regions" << std::endl;
  for (size_t i_reg=0; i_reg<regions.size(); i_reg++)
    {
      std::cout << "i_reg = " << i_reg << std::endl;
      regions[i_reg].prepare_idx();
      // regions[i_reg].print_region_map();
    }

}


// initialize the loops based on the set of binding regions
void loop_topology::initialize_loops(int N_loops, int min_dist)
{
  // clear the loops and binding regions
  loops.clear();

  size_t N_regions = regions.size();

  // std::cout << "at start of initialize_loops" << std::endl;
  // for (size_t i_reg=0; i_reg<regions.size(); i_reg++)
  //   {
  //     std::cout << "i_reg = " << i_reg << std::endl;
  //     regions[i_reg].print_region_map();
  //   }

  // only select binding regions that can contain a full loop (hinge and anchor) upon initialization

  int total_binding_region_size = 0;

  for (size_t i_region=0; i_region<N_regions; i_region++)
    {
      if (regions[i_region].get_size() > 2*min_dist)
	{
	  total_binding_region_size += regions[i_region].get_size();
	}
    }

  std::uniform_int_distribution<int> unif_dist(1,total_binding_region_size);
  std::vector<int> a_dist;

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


  int a_reg, h_reg;

  std::uniform_real_distribution<double> dir_dist(0.0, 1.0);
  double r_d;
  int d, r_mono, a_mono, h_mono;

  // select an anchor for each loop
  for (int i_loop=0; i_loop<N_loops; i_loop++)
    {
      // get the anchor region
      a_reg = loops[i_loop].get_a_region();

      // select a random monomer for the anchor
      std::uniform_int_distribution<int> unif_dist(0,regions[a_reg].get_size()-1);
      r_mono = unif_dist(rand_eng);
      a_mono = regions[a_reg].get_mono_pos(r_mono);
      loops[i_loop].set_a(a_mono);

      // bind the anchor
      loops[i_loop].set_a_bound(true);

      // get the hinge region
      h_reg = loops[i_loop].get_h_region();

      // based on the position of the anchor and a minimum distance between the anchor and hinge, select a direction for the hinge to travel
      r_d = dir_dist(rand_eng);
      d = regions[h_reg].select_direction(a_mono,
					  min_dist,
					  r_d);
      loops[i_loop].set_d(d);

      // select a hinge based on the anchor position and direction
      h_mono = regions[h_reg].get_relative_monomer_pos(a_mono,
						       min_dist,
						       d);
      loops[i_loop].set_h(h_mono);

      // std::cout << "a=" << a_mono << ", h=" << h_mono << ", d=" << d << std::endl;
      // regions[a_reg].print_region_map();

      // bind the hinge
      loops[i_loop].set_h_bound(true);
    }

}


void loop_topology::initialize_loops_on_periphery(int N_loops, int min_dist, double radius)
{
  std::cout << "Setting the anchors on periphery..." << std::endl;
  // Clear existing loops
  loops.clear();

  std::vector<std::pair<int, int>> candidate_monomers; // Pair of (region_index, monomer_position)

  // Collect all monomers from all regions that lie outside the specified radius
  for (size_t i_region = 0; i_region < regions.size(); i_region++)
  {
    binding_region region = regions[i_region];
    for (int i = 0; i < region.get_size(); i++)
    {
      int mono_pos = region.get_mono_pos(i);
	  double mono_dist = sqrt(vqm.v_L2(coords[mono_pos]));
      if (mono_dist > radius) // Replace with actual distance calculation
      {
      	// std::cout << "Accepted a monomer with distance from center = " << mono_dist << std::endl;
        candidate_monomers.emplace_back(i_region, mono_pos);
      }
    }
  }

  std::cout << "N_candidate_monomers = " << candidate_monomers.size() << std::endl;

  // Check if there are enough candidates
  if (candidate_monomers.size() < static_cast<std::size_t>(N_loops))
  {
    throw std::runtime_error("Not enough monomers outside the specified radius to initialize loops.");
  }

  // Shuffle the candidate pool to ensure randomness
  std::shuffle(candidate_monomers.begin(), candidate_monomers.end(), rand_eng);

  // Select anchors for the loops
  for (int i_loop = 0; i_loop < N_loops; i_loop++)
  {
    int a_region = candidate_monomers[i_loop].first;
    int a_mono = candidate_monomers[i_loop].second;

    // Create a new loop with the same region for anchor and hinge (for now)
    loop l(a_region, a_region);
    l.set_a(a_mono);       // Set the anchor monomer
    l.set_a_bound(true);   // Bind the anchor
    loops.push_back(l);    // Add loop to the list
  }

  // Now assign hinges for the loops
  std::uniform_real_distribution<double> dir_dist(0.0, 1.0);
  for (int i_loop = 0; i_loop < N_loops; i_loop++)
  {
    int a_mono = loops[i_loop].get_a();
    int h_region = loops[i_loop].get_h_region();

    // Select direction for the hinge to travel
    double r_d = dir_dist(rand_eng);
    int d = regions[h_region].select_direction(a_mono, min_dist, r_d);
    loops[i_loop].set_d(d);

    // Select hinge monomer based on anchor and direction
    int h_mono = regions[h_region].get_relative_monomer_pos(a_mono, min_dist, d);
    loops[i_loop].set_h(h_mono);
    loops[i_loop].set_h_bound(true); // Bind the hinge
  }
}

void loop_topology::initialize_loops_interior(int N_bound, int N_unbound, int min_dist, double radius)
{
  std::cout << "Setting the anchors within the periphery (interior of cell)..." << std::endl;
  // Clear existing loops
  // loops.clear();

  int N_original_loops = loops.size();

  std::vector<std::pair<int, int>> candidate_monomers; // Pair of (region_index, monomer_position)

  // Collect all monomers from all regions that lie outside the specified radius
  std::cout << "Looping through regions looking for suitable monomers..." << std::endl;
  for (size_t i_region = 0; i_region < regions.size(); i_region++)
  {
  	std::cout << "Obtaining the region..." << std::endl;
    binding_region region = regions[i_region];
  	std::cout << "Region obtained, checking monomers..." << std::endl;
    for (int i = 0; i < region.get_size(); i++)
    {
      int mono_pos = region.get_mono_pos(i);
      // std::cout << "mono pos:" << mono_pos << std::endl;
	  double mono_dist = sqrt(vqm.v_L2(coords[mono_pos]));
      if (mono_dist < radius)
      {
      	 candidate_monomers.emplace_back(i_region, mono_pos);

      }
    }
  }

  std::cout << "N_candidate_monomers = " << candidate_monomers.size() << std::endl;

  // Check if there are enough candidates
  if (candidate_monomers.size() < static_cast<std::size_t>(N_bound+N_unbound))
  {
    throw std::runtime_error("Not enough monomers inside the specified radius and in anchor candidates to initialize loops.");
  }

  // Shuffle the candidate pool to ensure randomness
  std::shuffle(candidate_monomers.begin(), candidate_monomers.end(), rand_eng);

  // Select anchors for the loops
  for (int i_loop = 0; i_loop < N_bound; i_loop++)
  {
    int a_region = candidate_monomers[i_loop].first;
    int a_mono = candidate_monomers[i_loop].second;

    // Create a new loop with the same region for anchor and hinge (for now)
    loop l(a_region, a_region);
    l.set_a(a_mono);       // Set the anchor monomer
    l.set_a_bound(true);   // Bind the anchor
    loops.push_back(l);    // Add loop to the list
  }

  // Now assign hinges for the loops
  std::uniform_real_distribution<double> dir_dist(0.0, 1.0);
  for (int i_loop = N_original_loops; i_loop < N_original_loops + N_bound; i_loop++)
  {
    int a_mono = loops[i_loop].get_a();
    int h_region = loops[i_loop].get_h_region();

    // Select direction for the hinge to travel
    double r_d = dir_dist(rand_eng);
    int d = regions[h_region].select_direction(a_mono, min_dist, r_d);
    loops[i_loop].set_d(d);

    // Select hinge monomer based on anchor and direction
    int h_mono = regions[h_region].get_relative_monomer_pos(a_mono, min_dist, d);
    loops[i_loop].set_h(h_mono);
    loops[i_loop].set_h_bound(true); // Bind the hinge
  }

  // Add unbound loops
  for (int i_loop = N_bound; i_loop < N_bound+N_unbound; i_loop++)
  {
    int a_region = candidate_monomers[i_loop].first; // dummy region

    // Create a new loop with the same region for anchor and hinge (for now)
    loop l(a_region, a_region);
    l.set_a_bound(false);
    loops.push_back(l);    // Add loop to the list
  }
}

void loop_topology::initialize_loops_conc(int N_loops, int min_dist)
{
  std::cout << "Setting the anchors bound to linear strands of DNA (low concentration voxels)..." << std::endl;
  // Clear existing loops
  // loops.clear();
  int N_original_loops = loops.size();

  std::vector<std::pair<int, int>> candidate_monomers; // Pair of (region_index, monomer_position)

  // Only add to candidates if monomer is part of predetermined anchor_candidates
  for (size_t i_region = 0; i_region < regions.size(); i_region++)
  {
    binding_region region = regions[i_region];
    for (int i = 0; i < region.get_size(); i++)
    {
      int mono_pos = region.get_mono_pos(i);

      // Check if mono_pos is in anchor_candidates
      if (std::find(anchor_candidates.begin(), anchor_candidates.end(), mono_pos) != anchor_candidates.end())
      {
        candidate_monomers.emplace_back(i_region, mono_pos);
      }
    }
  }

  std::cout << "N_candidate_monomers = " << candidate_monomers.size() << std::endl;

  // Shuffle the candidate pool to ensure randomness
  std::shuffle(candidate_monomers.begin(), candidate_monomers.end(), rand_eng);

  // Select anchors for the loops
  for (int i_loop = 0; i_loop < N_loops; i_loop++)
  {
    int a_region = candidate_monomers[i_loop].first;
    int a_mono = candidate_monomers[i_loop].second;

    // Create a new loop with the same region for anchor and hinge (for now)
    loop l(a_region, a_region);
    l.set_a(a_mono);       // Set the anchor monomer
    l.set_a_bound(true);   // Bind the anchor
    loops.push_back(l);    // Add loop to the list
  }

  // Now assign hinges for the loops
  std::uniform_real_distribution<double> dir_dist(0.0, 1.0);
  for (int i_loop = N_original_loops; i_loop < N_original_loops + N_loops; i_loop++)
  {
    int a_mono = loops[i_loop].get_a();
    int h_region = loops[i_loop].get_h_region();

    // Select direction for the hinge to travel
    double r_d = dir_dist(rand_eng);
    int d = regions[h_region].select_direction(a_mono, min_dist, r_d);
    loops[i_loop].set_d(d);

    // Select hinge monomer based on anchor and direction
    int h_mono = regions[h_region].get_relative_monomer_pos(a_mono, min_dist, d);
  	std::cout << "Hinge set to " << h_mono << " for anchor " << a_mono << std::endl;
    loops[i_loop].set_h(h_mono);
    loops[i_loop].set_h_bound(true); // Bind the hinge
  }
}


// update the loops
void loop_topology::update_loops(int ext_max, int min_dist, double p_unbinding, double p_off, double p_on, double r_g, double threshold)
{
  bool anchor_unbind_rebind = false;
  std::cout << "---[ updating loops ]---" << std::endl;
  std::uniform_real_distribution<double> unif_dist(0.0, 1.0);
  size_t N_regions = regions.size();
  int total_binding_region_size = 0;

  for (size_t i_region=0; i_region<N_regions; i_region++)
  {
      if (regions[i_region].get_size() > 2*min_dist)
	{
	  total_binding_region_size += regions[i_region].get_size();
	}
  }

  std::uniform_int_distribution<int> anchor_dist(1,total_binding_region_size);
  std::vector<int> intra_updates;
  std::vector<int> new_bound;
  std::vector<std::vector<int>> inter_updates;
  int a_mono, h_mono, d;
  size_t a_reg, h_reg;
  vec a_coord;
  bool a_bound, h_bound;
  double r_intra, r_d;
  size_t r_inter, accumulator;
  int h_intra, h_inter;

  // loop over the loops, deciding which bound loops to unbind and which unbound loops to bind
  for (size_t i_loop=0; i_loop<loops.size(); i_loop++)
    {
  	  anchor_unbind_rebind = false;
  	  a_bound = loops[i_loop].get_a_bound();

      if (a_bound == true)
        {
      	  vec h_pos = coords[loops[i_loop].get_h()];
      	  vec a_pos = coords[loops[i_loop].get_a()];
      	  double a_h_dist = sqrt(vqm.v_L2(vqm.v_axpy(-1, h_pos, a_pos)));
      	  double a_origin_dist = sqrt(vqm.v_L2(a_pos));
          if (unif_dist(rand_eng) < p_off) // anchor is bound, but it becomes unbound due to random unbinding
            {
			  loops[i_loop].set_a_bound(false);
          	  int a_mono = loops[i_loop].get_a();
          	  std::cout << "Anchor came unbound from monomer " << a_mono << std::endl;
            }
      	  else if (a_h_dist > 200) // anchor is bound, but it becaoems unbound because anchor-hinge pair overstreched
      	  {
      		  loops[i_loop].set_a_bound(false);
      		  int a_mono = loops[i_loop].get_a();
      		  std::cout << "Anchor came unbound from monomer " << a_mono << " because anchor-hinge pair was stretched > 20 nm. Rebinding elsewhere..." << std::endl;
      	  	  anchor_unbind_rebind = true;
      	  }
      	  else if (a_origin_dist > threshold) // anchor is bound, but it becaoems unbound because anchor-hinge pair overstreched
      	  {
      		loops[i_loop].set_a_bound(false);
      		int a_mono = loops[i_loop].get_a();
      		std::cout << "Anchor came unbound from monomer " << a_mono << " because anchor-hinge pair was too close to boundary (" << a_origin_dist << "). Rebinding elsewhere..."  << std::endl;
      	  	anchor_unbind_rebind = true;
      	  }
      	  else if (loops[i_loop].get_a_region() != loops[i_loop].get_h_region()) // replisome disrupts SMC, causing it to unbind
      	  {
      		loops[i_loop].set_a_bound(false);
      		int a_mono = loops[i_loop].get_a();
      		std::cout << "Anchor came unbound from monomer " << a_mono << " because a replication fork forced the SMC off the DNA. Rebinding elsewhere..."  << std::endl;
      		anchor_unbind_rebind = true;
      	  }
          for (size_t j_loop = 0; j_loop < i_loop; j_loop++) {
    		// If we are comparing two bound SMC, get anchor and hinge positions for both loops
          if (loops[j_loop].get_a_bound()) {
    		int a_mono = loops[i_loop].get_a();
    		int h_mono = loops[i_loop].get_h();
   		    int a_mono2 = loops[j_loop].get_a();
    		int h_mono2 = loops[j_loop].get_h();

    		// Check for overlap in a circular chromosome
    		bool overlap = false;
            int N_extruded_mono1 = abs(a_mono - h_mono);
            int N_extruded_mono2 = abs(a_mono2 - h_mono2);
  		    // int dir1 = (h_mono - a_mono) / N_extruded_mono1;
            // int dir2 = (h_mono2 - a_mono2) / N_extruded_mono2;
  	        bool wrap1 = N_extruded_mono1 > 40000; // we can use the ter_crossing logic
            bool wrap2 = N_extruded_mono2 > 40000; // we can use the ter_crossing logic

            int start1 = std::min(a_mono, h_mono);
        	int end1   = std::max(a_mono, h_mono);
       		int start2 = std::min(a_mono2, h_mono2);
        	int end2   = std::max(a_mono2, h_mono2);

    		if (!wrap1 && !wrap2) {
        		// Both intervals are "normal" (non-wrapping), so sort the endpoints.
        		overlap = (start1 <= end2 && end1 >= start2);
    		}
    		else if (wrap1 && wrap2) {
        		// If both intervals wrap, they both cross the terminus.
          		overlap = true;
    		}
    		else if (wrap1) {
        		// Interval 1 wraps (i.e. from start1 to chromo_length-1 and from 0 to end1)
        		// Overlap occurs if no overlap bewteen start1 to 54338, or 0 to end1
                int start1 = std::max(a_mono, h_mono);
        	    int end1   = std::min(a_mono, h_mono);
        		overlap = (start2 <= start1 && end2 >= start1)
                         || (start2 <= end1 && end2 >= end1)
                         || (start2 >= start1 && end2 <= 54338)
                         || (start2 >= 0 && end2 <= end1);
   			 }
    		else { // wrap2 is true (and wrap1 is false)
                int start2 = std::max(a_mono2, h_mono2);
        	    int end2   = std::min(a_mono2, h_mono2);
        		overlap = (start1 <= start2 && end1 >= start2)
                         || (start1 <= end2 && end1 >= end2)
                         || (start1 >= start2 && end1 <= 54338)
                         || (start1 >= 0 && end1 <= end2);
    		}

    		if (overlap) {
        		loops[i_loop].set_a_bound(false);
      			int a_mono = loops[i_loop].get_a();
      			std::cout << "Anchor came unbound from monomer " << a_mono << " because it overlapped with another SMC. Rebinding elsewhere..."  << std::endl;
      			anchor_unbind_rebind = true;
    		}
		  }
		  }

        }
      if (a_bound == false || anchor_unbind_rebind) {
          if (unif_dist(rand_eng) < p_on || anchor_unbind_rebind) // anchor is unbound, but it becomes bound
            {
              int a_mono = anchor_dist(rand_eng);

			  loops[i_loop].set_a_bound(true);
			  if (!anchor_unbind_rebind) {
          		  int accumulator = 0;
      			  for (size_t i_region=0; i_region<N_regions; i_region++)
	 			  {
	  			  if (regions[i_region].get_size() > 2*min_dist)
	  			  {
	  			   accumulator += regions[i_region].get_size();
	  			  }
	  			  if (a_mono <= accumulator)
	  			  {
	  			    loops[i_loop].set_a_region(i_region);
	  		  		std::cout << "Anchor/Hinge binding to region " << i_region << std::endl;
	 			    break;
				  }
				  }
			  } else {
			  	  std::cout << "Anchor/Hinge rebinds to the same region, " << loops[i_loop].get_a_region() << std::endl;
          	  }

          	  // get the anchor region
          	  a_reg = loops[i_loop].get_a_region();

          	  // select a random monomer for the anchor
          	  std::uniform_int_distribution<int> unif_dist(0,regions[a_reg].get_size()-1);
          	  int r_mono = unif_dist(rand_eng);
          	  a_mono = regions[a_reg].get_mono_pos(r_mono);
          	  loops[i_loop].set_a(a_mono);

              // Now assign hinges for the loops
          	  loops[i_loop].set_h_region(a_reg);
          	  int h_region = loops[i_loop].get_a_region();

              // Select direction for the hinge to travel
              double r_d = unif_dist(rand_eng);
              int d = regions[h_region].select_direction(a_mono, min_dist, r_d);
              loops[i_loop].set_d(d);

              // Select hinge monomer based on anchor and direction
              int h_mono = regions[h_region].get_relative_monomer_pos(a_mono, min_dist, d);
  	          std::cout << "Anchor bound. Hinge set to " << h_mono << " for anchor " << a_mono << std::endl;
              loops[i_loop].set_h(h_mono);
              loops[i_loop].set_h_bound(true); // Bind the hinge
			  new_bound.push_back(i_loop); // keep track that this one just bound
            }
      }
    }

  // loop over the loops to update hinge positions
  for (size_t i_loop=0; i_loop<loops.size(); i_loop++)
    {

      // get the anchor, anchor region, and anchor coordinate
      a_mono = loops[i_loop].get_a();
      a_reg = static_cast<size_t>(loops[i_loop].get_a_region());
      a_bound = loops[i_loop].get_a_bound();
      a_coord = coords[a_mono];

      // get the hinge, hinge region, and hinge binding state
      h_mono = loops[i_loop].get_h();
      h_reg = static_cast<size_t>(loops[i_loop].get_h_region());
      h_bound = loops[i_loop].get_h_bound();
      // Check if the anchor has bound this timestep
      int new_bound_anchor = (std::find(new_bound.begin(), new_bound.end(), i_loop) != new_bound.end());

	  if (a_bound && !new_bound_anchor)
    {

      // hinge is already bound
      if (h_bound == true)
	{

	  // get the direction
	  d = loops[i_loop].get_d();

	  // update the proximities in the hinge region
	  // std::cout << "h_reg = " << h_reg << std::endl;
	  // regions[h_reg].print_region_map();
	  regions[h_reg].update_proximities(r_g,a_coord,coords);

	  if (a_reg == h_reg) regions[h_reg].filter_proximities_near_a(min_dist,a_mono);

	  // get the intra updates
	  intra_updates = regions[h_reg].get_and_filter_intra_candidates(ext_max,h_mono,d);

	  // std::cout << "intra_updates, " << h_mono << ", " << intra_updates.size() << " candidates" << std::endl;
	  // for (size_t i_intra=0; i_intra<intra_updates.size(); i_intra++)
	  //   {
	  //     std::cout << intra_updates[i_intra] << std::endl;
	  //   }

	  if (unif_dist(rand_eng) > p_unbinding)
	    {

	      // select an updated hinge for 1D motion along strand
	      r_intra = unif_dist(rand_eng);
	      //std::cout << "r_intra = " << r_intra << std::endl;
	      h_intra = intra_updates[step_dist.get_k(r_intra,static_cast<int>(intra_updates.size()))];
	      //std::cout << "h_intra_final = " << h_intra << std::endl;

	      loops[i_loop].set_h(h_intra);

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
		      regions[i_reg].filter_proximities_near_a(min_dist,a_mono);
		    }

		  // add the region's candidates to the total set of possible inter updates
		  inter_updates.push_back(regions[i_reg].get_inter_candidates());

		  N_inter_total += inter_updates[i_reg].size();

		} // end loop to determine inter candidates

	      if (N_inter_total > 0)
		{
		  std::uniform_int_distribution<size_t> unif_dist(0,N_inter_total);

		  accumulator = 0;
		  r_inter = unif_dist(rand_eng);
		  std::cout << "r_inter = " << r_inter  << ", N_inter_total = " << N_inter_total << std::endl;

		  // determine the candidate inter-strand hinge updates for all the regions
		  for (size_t i_reg=0; i_reg<regions.size(); i_reg++)
		    {

		      if ((accumulator + inter_updates[i_reg].size()) > r_inter)
			{
			  h_inter = static_cast<int>(r_inter-accumulator);
			  // set the hinge
			  std::cout << "accumulator = " << accumulator  << ", i_reg = " << i_reg << std::endl;
			  std::cout << "size = " << inter_updates[i_reg].size()  << ", h_inter = " << h_inter << std::endl;
			  loops[i_loop].set_h(inter_updates[i_reg][h_inter]);
			  // set the hinge region
			  loops[i_loop].set_h_region(i_reg);
			  // bind the hinge
			  loops[i_loop].set_h_bound(true);
			  // set the direction
			  d = 1;
			  r_d = unif_dist(rand_eng);
			  if (r_d < 0.5) d = -1;
			  loops[i_loop].set_d(d);
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
		  regions[i_reg].filter_proximities_near_a(min_dist,a_mono);
		}

	      // add the region's candidates to the total set of possible inter updates
	      inter_updates.push_back(regions[i_reg].get_inter_candidates());

	      N_inter_total += inter_updates[i_reg].size();

	    } // end loop over regions

	  if (N_inter_total > 0)
	    {
	      std::uniform_int_distribution<int> unif_dist(0,N_inter_total);

	      accumulator = 0;
	      r_inter = unif_dist(rand_eng);

	      // determine the candidate inter-strand hinge updates for all the regions
	      for (size_t i_reg=0; i_reg<regions.size(); i_reg++)
		{

		  if ((accumulator + inter_updates[i_reg].size()) >= r_inter)
		    {
		      h_inter = static_cast<int>(r_inter-accumulator);
		      // set the hinge
		      loops[i_loop].set_h(inter_updates[i_reg][h_inter]);
		      // set the hinge region
		      loops[i_loop].set_h_region(i_reg);
		      // bind the hinge
		      loops[i_loop].set_h_bound(true);
		      // set the direction
		      d = 1;
		      r_d = unif_dist(rand_eng);
		      if (r_d < 0.5) d = -1;
		      loops[i_loop].set_d(d);
		    }
		  else
		    {
		      accumulator += inter_updates[i_reg].size();
		    }

		}

	    } // end conditional for nonzero inter candidates

	} // end conditional for bound/unbound hinges

    } // end conditional for bound anchors

    } // end loop over loops
}


// get the loops
std::vector<loop> loop_topology::get_loops()
{
  return loops;
}

// get the binding regions
std::vector<binding_region> loop_topology::get_regions()
{
  return regions;
}

void loop_topology::set_anchor_candidates(std::vector<int> anchor_candidates)
{
  this->anchor_candidates = anchor_candidates;
}

void loop_topology::update_anchor_hinge_regions()
{

  std::cout << "Updating anchor hinge pairs..." <<std::endl;
  for (size_t i_loop=0; i_loop<loops.size(); i_loop++)
  {
  	std::cout << "Updating loop " << i_loop << std::endl;
    if (loops[i_loop].get_a_bound()) {
	    int a_mono = loops[i_loop].get_a();
  		int h_mono = loops[i_loop].get_h();
    	std::cout << "a_mono: " << a_mono << ", h_mono: " << h_mono << std::endl;
  		int old_a_region = loops[i_loop].get_a_region();
  		int old_h_region = loops[i_loop].get_h_region();
  		std::cout << "old_a_region: " << old_a_region << ", old_h_region: " << old_h_region << std::endl;
  		int new_a_region = -1;
  		int new_h_region = -1;

  		for (size_t i_reg = 0; i_reg < regions.size(); i_reg++) {
        if (regions[i_reg].is_ter_crossing()) {  // Check if this region crosses the terminus
        // First part of the split region
        if ((a_mono >= regions[i_reg].get_ul()) && (a_mono <= regions[i_reg].get_mid_ul())) {
            new_a_region = i_reg;
        }
        if ((h_mono >= regions[i_reg].get_ul()) && (h_mono <= regions[i_reg].get_mid_ul())) {
            new_h_region = i_reg;
        }
        // Second part of the split region
        if ((a_mono >= regions[i_reg].get_mid_ll()) && (a_mono <= regions[i_reg].get_ll())) {
            new_a_region = i_reg;
        }
        if ((h_mono >= regions[i_reg].get_mid_ll()) && (h_mono <= regions[i_reg].get_ll())) {
            new_h_region = i_reg;
        }
        } else {
        // Non-split region logic (ter_crossing == false)
        if ((a_mono >= regions[i_reg].get_ll()) && (a_mono <= regions[i_reg].get_ul())) {
            new_a_region = i_reg;
        }
        if ((h_mono >= regions[i_reg].get_ll()) && (h_mono <= regions[i_reg].get_ul())) {
            new_h_region = i_reg;
        }
        }
        }
	    if (new_a_region == -1) {
          std::cerr << "Error: No region found for anchor at index " << a_mono << ", setting to 1" << std::endl;
          new_a_region = 1;
        }
        if (new_h_region == -1) {
          std::cerr << "Error: No region found for hinge at index " << h_mono << ", setting to 1" << std::endl;
          new_h_region = 1;
        }
    	std::cout << "new_a_region: " << new_a_region << ", new_h_region: " << new_h_region << std::endl;
  		loops[i_loop].set_a_region(new_a_region);
  		loops[i_loop].set_h_region(new_h_region);

  		if (new_a_region != old_a_region) {
  		  std::cout << "The region for anchor at index " << a_mono << " changed from " << old_a_region << " to " << new_a_region <<std::endl;
  		}
  		if (new_h_region != old_h_region) {
  			std::cout << "The region for hinge at index " << h_mono << " changed from " << old_h_region << " to " << new_h_region <<std::endl;
  		}

    	if (new_a_region == 2) {
    		loops[i_loop].set_a(loops[i_loop].get_a()+20);
    		std::cout << "Due to replication, index for anchor changed from " << a_mono << " to " << loops[i_loop].get_a() <<std::endl;
    	}
    	if (new_h_region == 2) {
    		loops[i_loop].set_h(loops[i_loop].get_h()+20);
    		std::cout << "Due to replication, index for hinge changed from " << h_mono << " to " << loops[i_loop].get_h() <<std::endl;
    	}
    }
  }

}
