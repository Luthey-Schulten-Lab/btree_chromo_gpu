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
  bool fork_exists;
  
  for (fork_rho f_r_f : final_state.transforms)
    {
      
      fork_exists = false;
      
      for (fork_rho f_r_i : initial_state.transforms)
	{
	  if (f_r_f.fork == f_r_i.fork)
	    {
	      f_r_temp = f_r_f;
	      f_r_temp.rho_cw -= f_r_i.rho_cw;
	      f_r_temp.rho_ccw -= f_r_i.rho_ccw;
	      fork_exists = true;
	      break;
	    }
	}

      if (fork_exists == false)
	{
	  f_r_temp = f_r_f;
	}

      diff_tr.push_back(f_r_temp);
      
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
	}

      for (int i=0; i<N_transforms; i++)
	{
	  for (int j=0; j<N_final; j++)
	    {
	      m[i][j] = new int[3];
	    }
	}

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
	  delete [] m[i];
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
  int N_prev, N_next;

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

  N_initial = prev_bt.total_size();
  N_final = next_bt.total_size();
  N_transforms = static_cast<int>(diff_tr.size());

  initialize_map();

  for (int i=0; i<N_transforms; i++)
    {
      // set the next binary tree state to the previous state
      next_bt.prepare_state(prev_bt.dump_state());

      // apply the transformation
      next_bt.single_transform(diff_tr[i]);

      // solve the new topology
      next_bt.solve_topology();

      // print the new btree
      next_bt.print_tree();

      // determine the mapping
      N_prev = prev_bt.total_size();
      N_next = next_bt.total_size();
      N_new[i] = N_next;

      cout << "N_prev = " << N_prev << endl;
      cout << "N_new = " << N_new[i] << endl;

      // set the previous binary tree state to the next state
      prev_bt.prepare_state(next_bt.dump_state());

      // solve the new topology
      prev_bt.solve_topology();
      
    }
  
  return 0;
  
}
