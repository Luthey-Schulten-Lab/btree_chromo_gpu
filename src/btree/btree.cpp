#include <iostream>

#include <btree/btree.hpp>

// constructor
btree::btree()
{
  root = NULL;
}


// destructor
btree::~btree()
{
  destroy_tree(root);
}


// function to initialize binary tree
void btree::initialize_tree(int s)
{
  if (root == NULL)
    {
      root = new node;
      initialize_branch(root,0,s);
    }
}


// function to prepare a state of the binary tree
void btree::prepare_state(btree_state st)
{

  int error_code;
  
  initialize_tree(st.size);
  for (fork_rho f_r: st.fork_rhos)
    {
      
      error_code = grow_at_branch(f_r.fork,f_r.rho);

      if (error_code == -1)
	{
	  cout << "ERROR: impossible state, destroying tree" << endl;
	  destroy_tree();
	  return;
	}
      
    }
}


// function to prepare a state of the binary tree
btree_state btree::dump_state()
{

  btree_state st;
  fork_rho f_r;
  
  if (root != NULL)
    {
      st.size = root->size;

      // add completed forks to state
      for (string s: get_completed_forks())
	{
	  
	  f_r.fork = s;
	  f_r.rho = get_branch(s)->rho;

	  st.fork_rhos.push_back(f_r);
	  
	}

      // add active forks to state
      for (string s: get_active_forks())
	{
	  
	  f_r.fork = s;
	  f_r.rho = get_branch(s)->rho;

	  st.fork_rhos.push_back(f_r);
	  
	}
      
    }
  else
    {
      st.size = -1;
    }
  
  return st;
}


// function to apply transforms to the binary tree
void btree::apply_transforms(btree_transforms tr)
{

  int error_code;
  btree_state st = dump_state();
  
  for (fork_rho f_r: tr.fork_rhos)
    {
      
      error_code = grow_at_branch(f_r.fork,f_r.rho);

      if (error_code == -1)
	{
	  cout << "ERROR: impossible state, destroying tree and returning to initial state" << endl;

	  // destroy the impossible tree and return to the initial state
	  destroy_tree();
	  prepare_state(st);
	  return;
	}
      
    }
}


// function to initialize a branch
void btree::initialize_branch(node *branch, int g, int s)
{
  branch->gen = g;
  branch->size = s;
  branch->rho = 0;
  branch->leaf = true;
  branch->parent = NULL;
  branch->left = NULL;
  branch->right = NULL;
}


// function to begin a branching procedure if it is valid
int btree::branch(string loc)
{
  node *branch;

  branch = get_branch(loc);

  if (branch == NULL)
    {
      cout << "invalid branch ("
	   << loc
	   << ")"
	   << endl;
      
      return -1;
    }
  else
    {
      
      cout << "branching at ("
	   << loc
	   << ")"
	   << endl;
      
      split_branch(branch);
      
      return 0;
    }
}


// function to create children at branch
void btree::split_branch(node *branch)
{
  if (branch->leaf == true)
    {
      branch->leaf = false;
  
      branch->left = new node;
      initialize_branch(branch->left,branch->gen+1,branch->size);
      branch->left->parent = branch;
  
      branch->right = new node;
      initialize_branch(branch->right,branch->gen+1,branch->size);
      branch->right->parent = branch;
    }
}


// function to increase rho in a manner dependent on the possible growth
int btree::grow_at_branch(string loc, int r)
{
  int max_growth, growth, rem, error_code;
  node *g_branch;

  error_code = branch(loc);
  if (error_code < 0)
    {
      return error_code;
    }

  // grow the branch if r > 0
  if (r > 0)
    {
      g_branch = get_branch(loc);

      if (g_branch->parent == NULL)
	{
	  max_growth = g_branch->size;
	}
      else
	{
	  max_growth = g_branch->parent->rho;
	}
  
      growth = min(max_growth-g_branch->rho,r);
      g_branch->rho += growth;
      rem = r - growth;
    }
  else // do nothing if r == 0
    {
      rem = 0;
    }
  
  return rem; // return the remainder of the growth
}


// function get pointer for a branch given the code
node *btree::get_branch(string loc)
{

  int n = loc.length();
  int i;

  node *tar;

  tar = root;

  if (n == 1)
    {
      return tar;
    }
  
  i = 1;
  while ((i < n) && (tar != NULL))
    {
      tar = parse_dir(tar,loc[i]);
      i += 1;
    }
  return tar;
  
}

node *btree::parse_dir(node *branch, char d)
{
  if (branch->leaf == false)
    {
      if (d == 'l')
	{
	  return branch->left;
	}
      else if (d == 'r')
	{
	  return branch->right;
	}
    }
  return NULL;
}


// functions to count the total number of leaves
int btree::count_total_leaves()
{
  return count_leaves("m");
}

int btree::count_leaves(string loc)
{
  return leaf_counter(get_branch(loc));
}

int btree::leaf_counter(node *branch)
{
  int leaf_count = 0;
  if (branch->leaf == true)
    {
      leaf_count += 1;
    }
  else
    {
      leaf_count += leaf_counter(branch->left);
      leaf_count += leaf_counter(branch->right);
    }
  return leaf_count;
}


// functions to calculate the total number of forks
int btree::count_total_forks()
{
  return count_total_forks("m");
}

int btree::count_total_forks(string loc)
{
  return total_fork_counter(get_branch(loc));
}

int btree::total_fork_counter(node *branch)
{
  int fork_count = 0;
  if (branch->leaf == false)
    {
      fork_count += 1;
      fork_count += total_fork_counter(branch->left);
      fork_count += total_fork_counter(branch->right);
    }
  return fork_count;
}


// counting functions for completed forks
int btree::count_completed_forks()
{
  return count_completed_forks("m");
}

int btree::count_completed_forks(string loc)
{
  return completed_fork_counter(get_branch(loc));
}

int btree::completed_fork_counter(node *branch)
{
  int fork_count = 0;
  if (branch->leaf == false)
    {
      if (branch->rho == branch->size) fork_count += 1;
      fork_count += completed_fork_counter(branch->left);
      fork_count += completed_fork_counter(branch->right);
    }
  return fork_count;
}


// counting functions for active forks
int btree::count_active_forks()
{
  return count_active_forks("m");
}

int btree::count_active_forks(string loc)
{
  return active_fork_counter(get_branch(loc));
}

int btree::active_fork_counter(node *branch)
{
  int fork_count = 0;
  if (branch->leaf == false)
    {
      if (branch->rho < branch->size) fork_count += 1;
      fork_count += active_fork_counter(branch->left);
      fork_count += active_fork_counter(branch->right);
    }
  return fork_count;
}


// functions to get labels for leaves
vector<string> btree::get_leaves()
{
  vector<string> leaves;
  if (count_total_leaves() > 0)
    {
      leaves.push_back("m");
      leaves = traverse_leaves(leaves,get_branch("m"));
    }
  return leaves;
}

vector<string> btree::traverse_leaves(vector<string> leaves, node *branch)
{
  if (branch->leaf == false)
    {
      string temp_parent = leaves.back();
      leaves.back() += "l";
      leaves = traverse_leaves(leaves,branch->left);
      leaves.push_back(temp_parent+"r");
      leaves = traverse_leaves(leaves,branch->right);
    }
  return leaves;
}


// functions to get labels for completed forks
vector<string> btree::get_completed_forks()
{
  vector<string> forks;
  if (count_completed_forks() > 0)
    {
      forks = traverse_completed_forks(forks,0,get_branch("m"));
    }
  return forks;
}

vector<string> btree::traverse_completed_forks(vector<string> forks, int i, node *branch)
{
  
  // case for root node
  if (forks.size() == 0)
    {
      if ((branch->leaf == false) &&
	  (branch->rho == branch->size))
	{ 
	  forks.push_back("m");
	}
    }

  // case for other nodes
  if ((branch->leaf == false) &&
      (branch->rho == branch->size))
    {
      forks.push_back(forks[i]+"l");
      forks = traverse_completed_forks(forks,forks.size()-1,branch->left);
      forks.push_back(forks[i]+"r");
      forks = traverse_completed_forks(forks,forks.size()-1,branch->right);
    }
  else
    {
      forks.pop_back();
    }
  return forks;
}


// functions to get labels for active forks
vector<string> btree::get_active_forks()
{
  vector<string> forks;
  if (count_active_forks() > 0)
    {
      forks = traverse_active_forks(forks,"m",get_branch("m"));
    }
  return forks;
}

vector<string> btree::traverse_active_forks(vector<string> forks, string f, node *branch)
{
  
  if (branch->leaf == false)
    { 
      if (branch->rho < branch->size)
	{ 
	  forks.push_back(f);
	}
      forks = traverse_active_forks(forks,f+"l",branch->left);
      forks = traverse_active_forks(forks,f+"r",branch->right);
    }
  return forks;
}


// function used by destructor
void btree::destroy_tree(node *branch)
{
  if (branch->leaf == false)
    {
      destroy_tree(branch->left);
      destroy_tree(branch->right);
      delete branch;
    }
}

void btree::destroy_tree()
{
  if (root != NULL)
    {
      destroy_tree(root);
      root = NULL;
    }
}


// function to print branches of the tree
void btree::print_branch(node *branch)
{

  string gen_offset;

  gen_offset = "";
  for (int i=0; i<branch->gen; i++)
    {
      gen_offset += "  ";
    }
  gen_offset += "| ";
  
  cout << gen_offset
	    << "generation = "
	    << branch->gen
	    << endl;
  
  if (branch->leaf == false)
    {
      int max_size;
      if (branch->parent != NULL)
	{
	  max_size = branch->parent->rho;
	}
      else
	{
	  max_size = branch->size;
	}
      cout << gen_offset
		<< "rho = "
		<< branch->rho
		<< "/" << max_size
		<< endl;
    }

  if (branch->leaf == false)
    {
      cout << gen_offset
		<< "left branch"
		<< endl;
      print_branch(branch->left);
      cout << gen_offset
		<< "right branch"
		<< endl;
      print_branch(branch->right);
    }

}


// function to print the entire tree
void btree::print_tree()
{
  
  if (root != NULL)
    {

      int N_leaves = count_total_leaves();
      int N_forks = count_total_forks();
      int N_completed_forks = count_completed_forks();
      int N_active_forks = count_active_forks();
      
      cout << "printing tree with "
	   << N_leaves
	   << " leaves and "
	   << N_forks
	   << " forks"
	   << endl;
      
      cout << "fork breakdown: "
	   << N_completed_forks
	   << " completed, "
	   << N_active_forks
	   << " active"
	   << endl;
		     
      cout << "leaves: "
	   << endl;
      
      for (string s: get_leaves())
	{
	  cout << s << " ";
	}
      cout << endl;

      if (N_completed_forks > 0)
	{
	  cout << "completed forks: "
	       << endl;
      
	  for (string s: get_completed_forks())
	    {
	      cout << s << " ";
	    }
	  cout << endl;
	}

      if (N_active_forks > 0)
	{
	  cout << "active forks: "
	       << endl;
      
	  for (string s: get_active_forks())
	    {
	      cout << s << " ";
	    }
	  cout << endl;
	}
      
      print_branch(root);
    }
  else
    {
      cout << "tree does not exist" << endl;
    }
}

void btree::foo()
{
  cout << "btree\n";
}
