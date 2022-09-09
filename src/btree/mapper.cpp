#include <btree/mapper.hpp>

// constructor
mapper::mapper()
{
  m = nullptr;
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
  initial_bt.prepare_state(st);
}


// setter for final state
void mapper::set_final_state(btree_state st)
{
  final_bt.prepare_state(st);
}


// determine the difference in states
btree_transforms mapper::state_diff(btree_state initial_state, btree_state final_state)
{

  btree_transforms diff;
  fork_rho f_r_temp;
  bool fork_exists;
  
  for (fork_rho f_r_f : final_state.fork_rhos)
    {
      
      fork_exists = false;
      
      for (fork_rho f_r_i : initial_state.fork_rhos)
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

      diff.fork_rhos.push_back(f_r_temp);
      
    }

  return diff;
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
}


// prepare the mapping
void mapper::prepare_mapping()
{

  btree temp_bt;

  temp_bt.prepare_state(initial_bt.dump_state());
  
  // solve the topologies for the two states
  initial_bt.solve_topology();
  final_bt.solve_topology();

  // print the initial state
  cout << "initial tree\n" << endl;
  initial_bt.print_tree();

  // print the final state
  cout << "final tree\n" << endl;
  final_bt.print_tree();

  btree_transforms diff_tr = state_diff(initial_bt.dump_state(),final_bt.dump_state());

  N_initial = initial_bt.total_size();
  N_final = final_bt.total_size();
  N_transforms = static_cast<int>(diff_tr.fork_rhos.size());

  initialize_map();

  for (int i=0; i<N_transforms; i++)
    {
      temp_bt.single_transform(diff_tr.fork_rhos[i]);
      temp_bt.solve_topology();
      temp_bt.print_tree();
    }
  

  
}
