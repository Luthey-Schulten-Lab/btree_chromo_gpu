#include <btree/mapper.hpp>

// constructor
mapper::mapper()
{
  m = nullptr;
  N_new = nullptr;
  N_initial = -1;
  N_final = -1;
  N_transforms = -1;
}


// destructor
mapper::~mapper()
{
  destroy_map();
}


// setter for initial state
void mapper::set_initial_state(btree_state st)
{
  initial_st = st;
}


// setter for final state
void mapper::set_final_state(btree_state st)
{
  final_st = st;
}


// determine the difference in states
btree_transforms mapper::state_diff(btree_state initial_state, btree_state final_state)
{

  btree_transforms diff_tr;
  fork_rho f_r_temp;
  bool fork_exists, fork_growth;
  
  for (fork_rho f_r_f : final_state.transforms)
    {
      
      fork_exists = false;
      fork_growth = false;
      // cout << "f_r_f: " << f_r_f.fork << "_cw" << f_r_f.rho_cw << "_ccw" << f_r_f.rho_ccw << endl;
      
      for (fork_rho f_r_i : initial_state.transforms)
	{
	  // cout << "f_r_i: " << f_r_i.fork << "_cw" << f_r_i.rho_cw << "_ccw" << f_r_i.rho_ccw << endl;
	  
	  if (f_r_f.fork == f_r_i.fork)
	    {
	      fork_exists = true;
	      f_r_temp = f_r_f;
	      f_r_temp.rho_cw = max(0,f_r_temp.rho_cw-f_r_i.rho_cw);
	      f_r_temp.rho_ccw = max(0,f_r_temp.rho_ccw-f_r_i.rho_ccw);
	      // cout << f_r_temp.fork << "_cw" << f_r_temp.rho_cw << "_ccw" << f_r_temp.rho_ccw << endl;
	      if ((f_r_temp.rho_cw > 0) || (f_r_temp.rho_ccw > 0)) fork_growth = true;
	      break;
	    }
	  
	}

      if (fork_exists == false) f_r_temp = f_r_f;

      // cout << f_r_temp.fork << "_cw" << f_r_temp.rho_cw << "_ccw" << f_r_temp.rho_ccw << endl;
      if ((fork_exists == false) || (fork_growth == true)) diff_tr.push_back(f_r_temp);
      
    }

  return diff_tr;
}


// initialize the map
void mapper::initialize_map()
{
  destroy_map();
  if ((N_transforms > 0) &&
      (N_final > 0))
    {
	  
      m = new int**[N_transforms];

      for (int i=0; i<N_transforms; i++)
	{
	  m[i] = new int*[N_final];
	  for (int j=0; j<N_final; j++)
	    {
	      m[i][j] = new int[3];
	    }
	}

      // for (int i=0; i<N_transforms; i++)
      // 	{
      // 	  for (int j=0; j<N_final; j++)
      // 	    {
      // 	      m[i][j] = new int[3];
      // 	    }
      // 	}

      // initialize the map to zero
      for (int i=0; i<N_transforms; i++)
	{
	  for (int j=0; j<N_final; j++)
	    {
	      for (int k=0; k<3; k++)
		{
		  m[i][j][k] = 0;
		}
	    }
	}	  
    }

  if (N_transforms > 0)
    {
      N_new = new int[N_transforms];
      for (int i=0; i<N_transforms; i++)
	{
	  N_new[i] = 0;
	}
    }
  
}


// destroy the map
void mapper::destroy_map()
{
  if (m != nullptr)
    {

      for (int i=0; i<N_transforms; i++)
	{
	  for (int j=0; j<N_final; j++)
	    {
	      delete[] m[i][j];
	    }
	  delete[] m[i];
	}
      delete[] m;
      m = nullptr;
    }

  if (N_new != nullptr)
    {
      delete[] N_new;
      N_new = nullptr;
    }
}


// prepare the mapping
int mapper::prepare_mapping()
{

  btree prev_bt, next_bt;

  // test if states are compatible
  if (initial_st.size != final_st.size)
    {
      return 1;
    }

  // set the btrees to the initial and final states
  prev_bt.prepare_state(initial_st);
  next_bt.prepare_state(final_st);
  
  // solve the topologies for the two states
  prev_bt.solve_topology();
  next_bt.solve_topology();

  // print the initial state
  cout << "initial tree\n" << endl;
  prev_bt.print_tree();

  // print the final state
  cout << "final tree\n" << endl;
  next_bt.print_tree();

  btree_transforms diff_tr = state_diff(initial_st,final_st);

  // destroy the map if it already exists
  destroy_map();
  
  N_initial = prev_bt.total_size();
  N_final = next_bt.total_size();
  N_transforms = static_cast<int>(diff_tr.size());

  initialize_map();

  int N_prev, N_next;
  string map_fork;
  string l_d_lmax_next, r_d_lmax_next;
  string l_d_lmax_prev, r_d_lmax_prev;
  bool new_fork, daughter_found, leaf_match;
  vector<string> prev_leaves;
  vector<string> next_leaves;
  theta_topo topo_prev, topo_next, topo_special;

  for (int i_trans=0; i_trans<N_transforms; i_trans++)
    {
      // set the next binary tree state to the previous state
      next_bt.prepare_state(prev_bt.dump_state());

      // store the fork location
      map_fork = diff_tr[i_trans].fork;
      
      // apply the transformation
      next_bt.single_transform(diff_tr[i_trans]);

      // solve the new topology
      next_bt.solve_topology();

      // determine the mapping
      N_prev = prev_bt.total_size();
      N_next = next_bt.total_size();
      N_new[i_trans] = N_next;

      cout << "i_trans = " << i_trans << endl;
      cout << diff_tr[i_trans].fork << "_cw" << diff_tr[i_trans].rho_cw << "_ccw" << diff_tr[i_trans].rho_ccw << endl;
      cout << "N_prev = " << N_prev << endl;
      cout << "N_new = " << N_new[i_trans] << endl;

      // determine the leaves in the two states
      prev_leaves = prev_bt.get_leaves();
      next_leaves = next_bt.get_leaves();
      
      new_fork = false;
      if (prev_leaves.size() < next_leaves.size()) new_fork = true;

      // following the left branch of the fork, determine the leftmost daughter in the next state
      l_d_lmax_next = map_fork + "l";

      daughter_found = false;
      while (daughter_found == false)
	{
	  leaf_match = false;
	  for (string leaf : next_leaves)
	    {
	      if (l_d_lmax_next == leaf.substr(0,l_d_lmax_next.length()))
		{
		  leaf_match = true;
		  break;
		}
	    }
	  if (leaf_match == false)
	    {
	      daughter_found = true;
	    }
	  else
	    {
	      l_d_lmax_next += "l";
	    }
	}

      l_d_lmax_next = l_d_lmax_next.substr(0,l_d_lmax_next.length()-1);

      // following the right branch of the fork, determine the leftmost daughter in the next state
      r_d_lmax_next = map_fork + "r";

      daughter_found = false;
      while (daughter_found == false)
	{
	  leaf_match = false;
	  for (string leaf : next_leaves)
	    {
	      if (r_d_lmax_next == leaf.substr(0,r_d_lmax_next.length()))
		{
		  leaf_match = true;
		  break;
		}
	    }
	  if (leaf_match == false)
	    {
	      daughter_found = true;
	    }
	  else
	    {
	      r_d_lmax_next += "l";
	    }
	}

      r_d_lmax_next = r_d_lmax_next.substr(0,r_d_lmax_next.length()-1);

      // determine leftmost daughters from the fork in the previous state
      if (new_fork == true)
	{
	  l_d_lmax_prev = map_fork;
	  r_d_lmax_prev = "none";
	}
      else
	{
	  l_d_lmax_prev = l_d_lmax_next;
	  r_d_lmax_prev = r_d_lmax_next;
	}

      cout << map_fork << endl;
      cout << "l_d_lmax_prev = " << l_d_lmax_prev << endl;
      cout << "r_d_lmax_prev = " << r_d_lmax_prev << endl;
      cout << "l_d_lmax_next = " << l_d_lmax_next << endl;
      cout << "r_d_lmax_next = " << r_d_lmax_next << endl;

      // print the next btree
      next_bt.print_tree();

      int prev_mono, next_mono;
      int start_link_offset, end_link_offset; //, prev_size;

      if (new_fork == true)
	{
	  cout << "new leaf will be created" << endl;
	  
	  // loop over the next leaves
	  for (string leaf : next_leaves)
	    {
	      // cout << "leaf = " << leaf << endl;
	      
	      if (leaf == l_d_lmax_next) // leaf is leftmost descendant of left daughter at growth fork
		
		{

		  topo_next = next_bt.get_leaf_topo(leaf);
		  topo_prev = prev_bt.get_leaf_topo(l_d_lmax_prev);

		  for (int j=0; j<(topo_next.end-topo_next.start+1); j++)
		    {
		      prev_mono = topo_prev.start + j;
		      next_mono = topo_next.start + j;
		      m[i_trans][next_mono][0] = next_mono; // mono in next state
		      m[i_trans][next_mono][1] = prev_mono; // mono in prev state
		      m[i_trans][next_mono][2] = 0; // direction
		    }

		}
	      
	      else if (leaf == r_d_lmax_next) // leaf is leftmost descendant of right daughter at growth fork
		
		{

		  topo_next = next_bt.get_leaf_topo(leaf);
		  topo_prev = prev_bt.get_leaf_topo(l_d_lmax_prev);

		  // replication was completed by the transform
		  if ((topo_next.start_link == topo_next.end) &&
		      (topo_next.end_link == topo_next.start))
		    {

		      for (int j=0; j<(topo_next.end-topo_next.start+1); j++)
			{
			  prev_mono = topo_prev.start + j;
			  next_mono = topo_next.start + j;
			  m[i_trans][next_mono][0] = next_mono; // mono in next state
			  m[i_trans][next_mono][1] = prev_mono; // mono in prev state
			  m[i_trans][next_mono][2] = -1; // direction
			  m[i_trans][prev_mono][2] = 1; // direction along opposite strand
			}
		      
		    }
		  
		  else // replication was not completed
		    
		    {
		      
		      for (int j=0; j<(topo_next.end-topo_next.start+1); j++)
			{
			  if ((topo_next.start_link <= topo_prev.mid) && (topo_next.end_link > topo_prev.mid))
			    {
			      prev_mono = (topo_next.start_link + 1) + j;
			    }
			  else if (topo_next.start_link > topo_prev.mid) // ccw past Ter
			    {
			      if (j < (topo_prev.end - topo_next.start_link))
				{
				  prev_mono = (topo_next.start_link + 1) + j;
				}
			      else
				{
				  prev_mono = topo_prev.start + j - (topo_prev.end - topo_next.start_link);
				}
			    }
			  else if (topo_next.end_link <= topo_prev.mid) // cw past ter
			    {
			      if (j < (topo_prev.end - topo_next.start_link))
				{
				  prev_mono = topo_next.start_link + j + 1;
				}
			      else
				{
				  prev_mono = topo_prev.start + j - (topo_prev.end - topo_next.start_link);
				}
			    }
			  next_mono = topo_next.start + j;
			  m[i_trans][next_mono][0] = next_mono; // mono in next state
			  m[i_trans][next_mono][1] = prev_mono; // mono in prev state
			  m[i_trans][next_mono][2] = -1; // direction
			  m[i_trans][prev_mono][2] = 1; // direction along opposite strand
			}
		      
		    } // end conditional for replication completion
		  
		}
	      
	      else // leaf is not subject to monomer additions
		
		{

		  topo_next = next_bt.get_leaf_topo(leaf);
		  topo_prev = prev_bt.get_leaf_topo(leaf);

		  for (int j=0; j<(topo_next.end-topo_next.start+1); j++)
		    {
		      prev_mono = topo_prev.start + j;
		      next_mono = topo_next.start + j;
		      m[i_trans][next_mono][0] = next_mono; // mono in next state
		      m[i_trans][next_mono][1] = prev_mono; // mono in prev state
		      m[i_trans][next_mono][2] = 0; // direction
		    }

		} // end conditional leaf type

	    } // end loop over new leaves

	} 
      
      else // no new leaf was created
	
	{
	  cout << "no new leaf will be created" << endl;
	  // loop over the next leaves
	  for (string leaf : next_leaves)
	    {
	      // cout << "leaf = " << leaf << endl;
	      
	      if (leaf == l_d_lmax_next) // leaf is leftmost descendant of left daughter at growth fork
		
		{

		  topo_next = next_bt.get_leaf_topo(leaf);
		  topo_prev = prev_bt.get_leaf_topo(l_d_lmax_prev);

		  for (int j=0; j<(topo_next.end-topo_next.start+1); j++)
		    {
		      prev_mono = topo_prev.start + j;
		      next_mono = topo_next.start + j;
		      m[i_trans][next_mono][0] = next_mono; // mono in next state
		      m[i_trans][next_mono][1] = prev_mono; // mono in prev state
		      m[i_trans][next_mono][2] = 0; // direction
		    }

		}
	      
	      else if (leaf == r_d_lmax_next) // leaf is leftmost descendant of right daughter at growth fork
		
		{

		  topo_next = next_bt.get_leaf_topo(leaf);
		  topo_prev = prev_bt.get_leaf_topo(r_d_lmax_prev);

		  // replication was completed by the transform
		  if ((topo_next.start_link == topo_next.end) &&
		      (topo_next.end_link == topo_next.start))
		    {

		      topo_special = next_bt.get_leaf_topo(l_d_lmax_next);

		      // ccw not past ter and cw not past ter in prev leaf
		      if ((topo_prev.start_link <= topo_special.mid) && (topo_prev.end_link > topo_special.mid))
			{

			  start_link_offset = topo_prev.start_link - topo_special.start;
			  end_link_offset = topo_special.end - topo_prev.end_link + 1;
			  cout << "start_link_offset = " << start_link_offset << endl;
			  cout << "end_link_offset = " << end_link_offset << endl;

			  for (int j=0; j<(topo_next.end-end_link_offset-(topo_next.start+start_link_offset)); j++)
			    {
			      prev_mono = topo_prev.start + j;
			      next_mono = (topo_next.start+start_link_offset) + j;
			      m[i_trans][next_mono][0] = next_mono; // mono in next state
			      m[i_trans][next_mono][1] = prev_mono; // mono in prev state
			      m[i_trans][next_mono][2] = 0; // direction
			    }

			  // replication fork traveling negtive direction along monomers towards ter
			  for (int j=0; j<start_link_offset; j++)
			    {
			      prev_mono = topo_special.start + j;
			      next_mono = topo_next.start + j;
			      m[i_trans][next_mono][0] = next_mono; // mono in next state
			      m[i_trans][next_mono][1] = prev_mono; // mono in prev state
			      m[i_trans][next_mono][2] = -1; // direction
			      m[i_trans][prev_mono][2] = 1; // direction along opposite strand
			    }

			  // replication fork traveling positive direction along monomers towards ter
			  for (int j=0; j<(end_link_offset+1); j++)
			    {
			      prev_mono = topo_prev.end_link + j - 1;
			      next_mono = topo_next.end-end_link_offset + j;
			      m[i_trans][next_mono][0] = next_mono; // mono in next state
			      m[i_trans][next_mono][1] = prev_mono; // mono in prev state
			      m[i_trans][next_mono][2] = -1; // direction
			      m[i_trans][prev_mono][2] = 1; // direction along opposite strand
			    }
		      
			}
		      else if (topo_prev.start_link > topo_special.mid) // ccw past Ter
			{
			      
			}
		      else if (topo_prev.end_link <= topo_special.mid) // cw past ter
			{
		      
			}

		    }
		  
		  else // replication was not completed
		    
		    {

		      topo_special = next_bt.get_leaf_topo(l_d_lmax_next);

		      // ccw not past ter and cw not past ter in prev leaf
		      if ((topo_prev.start_link <= topo_special.mid) && (topo_prev.end_link > topo_special.mid))
			{

			  start_link_offset = topo_prev.start_link - topo_next.start_link;
			  end_link_offset = topo_next.end_link - topo_prev.end_link;
			  cout << "start_link_offset = " << start_link_offset << endl;
			  cout << "end_link_offset = " << end_link_offset << endl;

			  for (int j=0; j<(topo_next.end-end_link_offset+1-(topo_next.start+start_link_offset)); j++)
			    {
			      prev_mono = topo_prev.start + j;
			      next_mono = (topo_next.start+start_link_offset) + j;
			      m[i_trans][next_mono][0] = next_mono; // mono in next state
			      m[i_trans][next_mono][1] = prev_mono; // mono in prev state
			      m[i_trans][next_mono][2] = 0; // direction
			    }

			  // replication fork traveling negtive direction along monomers towards ter
			  for (int j=0; j<start_link_offset; j++)
			    {
			      prev_mono = topo_next.start_link + j;
			      next_mono = topo_next.start + j;
			      m[i_trans][next_mono][0] = next_mono; // mono in next state
			      m[i_trans][next_mono][1] = prev_mono; // mono in prev state
			      m[i_trans][next_mono][2] = -1; // direction
			      m[i_trans][prev_mono][2] = 1; // direction along opposite strand
			    }

			  // replication fork traveling positive direction along monomers towards ter
			  for (int j=0; j<end_link_offset; j++)
			    {
			      prev_mono = topo_prev.end_link + j - 1;
			      next_mono = topo_next.end-end_link_offset + j + 1;
			      m[i_trans][next_mono][0] = next_mono; // mono in next state
			      m[i_trans][next_mono][1] = prev_mono; // mono in prev state
			      m[i_trans][next_mono][2] = -1; // direction
			      m[i_trans][prev_mono][2] = 1; // direction along opposite strand
			    }
		      
			}
		      else if (topo_prev.start_link > topo_special.mid) // ccw past Ter
			{
			      
			}
		      else if (topo_prev.end_link <= topo_special.mid) // cw past ter
			{
		      
			}
		      
		    } // end conditionals for replication completion
		  
		}
	      
	      else // leaf is not subject to monomer additions
		
		{

		  topo_next = next_bt.get_leaf_topo(leaf);
		  topo_prev = prev_bt.get_leaf_topo(leaf);

		  for (int j=0; j<(topo_next.end-topo_next.start+1); j++)
		    {
		      prev_mono = topo_prev.start + j;
		      next_mono = topo_next.start + j;
		      m[i_trans][next_mono][0] = next_mono; // mono in next state
		      m[i_trans][next_mono][1] = prev_mono; // mono in prev state
		      m[i_trans][next_mono][2] = 0; // direction
		    }

		} // end conditional for leaf type

	    } // end loop over leaves

	} // end conditional for a new leaf

      for (int k=0; k<N_new[i_trans]; k++)
	{
	  cout << k << " : "
	       << m[i_trans][k][0] << ","
	       << m[i_trans][k][1] << ","
	       << m[i_trans][k][2] << endl;
	}

      cout << "\n" << endl;

      // set the previous binary tree state to the next state
      prev_bt.prepare_state(next_bt.dump_state());

      // solve the new topology
      prev_bt.solve_topology();
      
    }
  
  return 0;
  
}


// get the map
vector<vector<array<int,3>>> mapper::get_map()
{
  vector<vector<array<int,3>>> map;
  vector<array<int,3>> map_per_transform;
  array<int,3> map_element;

  for (int i_trans=0; i_trans<N_transforms; i_trans++)
    {

      map_per_transform.clear();

      for (int i=0; i<N_new[i_trans]; i++)
	{

	  map_element[0] = m[i_trans][i][0];
	  map_element[1] = m[i_trans][i][1];
	  map_element[2] = m[i_trans][i][2];

	  map_per_transform.push_back(map_element);
	  
	}

      map.push_back(map_per_transform);
      
    }

  return map;
  
}
