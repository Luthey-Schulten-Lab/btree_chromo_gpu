#include <btree/btree.hpp>

// constructor
btree::btree()
{
  reset_root();
}

void btree::reset_root()
{
  root = nullptr;
}

// destructor
btree::~btree()
{
  // cout << "btree destructor before destroy" << endl;
  destroy_tree();
  // cout << "btree destructor after destroy" << endl;
}


// function to initialize binary tree
void btree::initialize_tree(int s)
{
  if (root != nullptr) destroy_tree();
  root = new node;
  initialize_branch(root,0,s);
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


// function to dump a state of the binary tree
btree_state btree::dump_state()
{

  btree_state st;
  fork_rho f_r;
  
  if (root != nullptr)
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


// function to write a binary tree state to a file
void btree::write_state(string st_filename, btree_state st)
{
  fstream st_file;

  st_file.open(st_filename, ios::out);

  if (!st_file)
    {
      cout << "ERROR: file not opened in write_state" << endl;
    }
  else
    {

      st_file << "size=" << st.size << endl;

      for (fork_rho f_r: st.fork_rhos)
	{
      
	  st_file << f_r.fork << "_" << f_r.rho << endl;
      
	}
    }

  st_file.close();
}


// function to read a binary tree state from a file
btree_state btree::read_state(string st_filename)
{
  fstream st_file;

  btree_state st;
  fork_rho f_r;
  
  string line;

  int delim;
  bool size_found = false;

  st_file.open(st_filename, ios::in);

  // cout << st_filename << endl;

  if (!st_file)
    {
      st.size = -1;
      cout << "ERROR: file not opened in read_state" << endl;
    }
  else
    {
      while (1)
	{
	  st_file >> line;
	  if (st_file.eof()) break;
	  
	  // cout << line << endl;

	  if (size_found == false)
	    {
	      
	      delim = line.find("=");
	      // cout << line.substr(0,delim) << endl;
	      // cout << line.substr(delim+1,line.length()) << endl;
	      
	      if (line.substr(0,delim) == "size")
		{
		  st.size = stoi(line.substr(delim+1,line.length()));
		  size_found = true;
		}
	      
	    }
	  else
	    {
	      
	      delim = line.find("_");
	      if (delim != -1)
		{
		  f_r.fork = line.substr(0,delim);
		  // cout << f_r.fork << endl;
		  f_r.rho = stoi(line.substr(delim+1,line.length()));
		  // cout << f_r.rho << endl;
		  st.fork_rhos.push_back(f_r);
		}
	      
	    }
      
	}
  
    }
  
  st_file.close();

  return st;
}


// function to read a binary tree state from a file
btree_transforms btree::read_transforms(string tr_filename)
{
  fstream tr_file;

  btree_transforms tr;
  fork_rho f_r;
  
  string line;

  int delim;

  tr_file.open(tr_filename, ios::in);

  // cout << tr_filename << endl;

  if (!tr_file)
    {
      cout << "ERROR: file not opened in read_transforms" << endl;
    }
  else
    {
      while (1)
	{
	  tr_file >> line;
	  if (tr_file.eof()) break;
	  
	  // cout << line << endl;
	      
	  delim = line.find("_");
	  if (delim != -1)
	    {
	      f_r.fork = line.substr(0,delim);
	      // cout << f_r.fork << endl;
	      f_r.rho = stoi(line.substr(delim+1,line.length()));
	      // cout << f_r.rho << endl;
	      tr.fork_rhos.push_back(f_r);
	    }

     	}
  
    }
  
  tr_file.close();

  return tr;
}


// function to initialize a branch
void btree::initialize_branch(node *branch, int g, int s)
{
  branch->gen = g;
  branch->size = s;
  branch->rho = 0;
  branch->leaf = true;
  branch->parent = nullptr;
  branch->left = nullptr;
  branch->right = nullptr;
  branch->topo.start = -1;
  branch->topo.end = -1;
  branch->topo.mid = -1;
  branch->topo.start_link = -1;
  branch->topo.end_link = -1;
}


// function to begin a branching procedure if it is valid
int btree::branch(string loc)
{
  node *branch;

  branch = get_branch(loc);

  if (branch == nullptr)
    {
      cout << "invalid branch ("
	   << loc
	   << ")"
	   << endl;
      
      return -1;
    }
  else
    {
      
      // cout << "branching at ("
      // 	   << loc
      // 	   << ")"
      // 	   << endl;
      
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

      if (g_branch->parent == nullptr)
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
  while ((i < n) && (tar != nullptr))
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
  return nullptr;
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


// calculate the total size of the system
int btree::total_size()
{
  return branch_size(get_branch("m"));
}

// calculate size of branch
int btree::branch_size(node *branch)
{
  int s = 0;

  if (branch->parent == nullptr)
    {
      s += branch->size;
    }
  
  if (branch->leaf == false)
    {
      s += branch->rho;
      s += branch_size(branch->left);
      s += branch_size(branch->right);
    }
  
  return s;
}


// function to partition units among branches
void btree::solve_topology()
{
  int temp_start, target_size;
  string free_leaf;
  string query_leaf;
  node *free_branch;
  node *query_branch;

  int N_leaves = count_total_leaves();
  vector<string> leaves = get_leaves();

  bool match;
  int j_q, j_f;

  int mid, d_start, d_end;

  temp_start = 0;

  for (int i=0; i<N_leaves; i++)
    {

      free_leaf = leaves[i];
      free_branch = get_branch(free_leaf);

      // first free branch is always circular
      if (i == 0)
	{
	  
	  target_size = free_branch->size;
	  
	}
      // find query branch supporting free branch
      else
	{

	  match = false;
	  query_leaf = leaves[i-1];

	  j_f = free_leaf.length() - 1;
	  while ((match == false) && (j_f >= 0))
	    {

	      j_q = query_leaf.length() - 1;
	      while ((match == false) && (j_q >= 0))
		{

		  if (free_leaf.substr(0,j_f) == query_leaf.substr(0,j_q))
		    {
		      query_branch = get_branch(query_leaf.substr(0,j_q));
		      match = true;
		    }

		  j_q -= 1;
	      
		}

	      j_f -= 1;
	      
	    }

	  // set the target size
	  target_size = query_branch->rho;

	  // descend to leftmost leaf of query branch
	  while (query_branch->leaf == false)
	    {
	      query_branch = query_branch->left;
	    }
	  
	}

      // set the start, end, and mid
      free_branch->topo.start = temp_start;
      free_branch->topo.end = temp_start + target_size - 1;
      free_branch->topo.mid = (free_branch->topo.start + free_branch->topo.end)/2;

      // increase the starting location by the added size
      temp_start += target_size;

      // create circular configuration if replication is complete
      if (target_size == free_branch->size)
	{
	  free_branch->topo.start_link = free_branch->topo.end;
	  free_branch->topo.end_link = free_branch->topo.start;
	}
      // create a theta configuration
      else
	{

	  d_start = target_size/2;
	  d_end = target_size - d_start;

	  mid = query_branch->topo.mid;

	  free_branch->topo.start_link = mid - d_start;
	  free_branch->topo.end_link = mid + d_end + 1;
	  
	}
      
    }
  
}


// function used to dump the topology to a file
void btree::dump_topology(string topo_filename, int idx)
{
  fstream topo_file;
  
  node *topo_branch;

  topo_file.open(topo_filename, ios::out);

  if (!topo_file)
    {
      cout << "ERROR: file not opened in dump_topology" << endl;
    }
  else
    {

      topo_file << "size=" << total_size() <<  endl;      
  
      for (string leaf: get_leaves())
	{
	  
	  topo_branch = get_branch(leaf);

	  topo_file << leaf
		    << "("
		    << (topo_branch->topo.end-topo_branch->topo.start+1)
		    << ")"
		    << "," << (topo_branch->topo.start_link + idx)
		    << "," << (topo_branch->topo.start + idx)
		    << "," << (topo_branch->topo.mid + idx)
		    << "," << (topo_branch->topo.end + idx)
		    << "," << (topo_branch->topo.end_link + idx)
		    << endl;
      
	}

    }

  topo_file.close();
}


// function to read chromosome regions from file
vector<chromo_region> btree::read_regions(string rg_filename, int idx)
{
  fstream rg_file;

  chromo_region c_r;
  vector<chromo_region> c_rs;

  string line;
  
  string rg_delim = ",";
  int delim;
  
  rg_file.open(rg_filename, ios::in);

  if (!rg_file)
    {
      cout << "ERROR: file not opened in read_regions" << endl;
    }
  else
    {
      while (1)
	{
	  rg_file >> line;
	  if (rg_file.eof()) break;

	  // region name
	  delim = line.find(rg_delim);
	  c_r.name = line.substr(0,delim);
	  line.erase(0,delim+1);

	  // start of region (inclusive)
	  delim = line.find(rg_delim);
	  c_r.start = stoi(line.substr(0,delim)) - idx;
	  line.erase(0,delim+1);

	  // end of region (inclusive)
	  c_r.end = stoi(line) - idx;

	  // initialize count to zero
	  c_r.count = 0;

	  c_rs.push_back(c_r);
	  
	}
  
    }
  
  rg_file.close();

  return c_rs;
  
}


// function to update region counts
void btree::update_region_counts(vector<chromo_region> &c_rs)
{

  int offset;
  int circ_size, exists_size;
  int start_circ, end_circ;
  node *count_branch;

  // zero counts before summation over leaves
  for (chromo_region &c_r: c_rs)
    {
      c_r.count = 0;
    }

  // iterate over leaves
  for (string leaf: get_leaves())
    {
	  
      count_branch = get_branch(leaf);

      // calculate appropriate offset from topology
      circ_size = count_branch->size;
      // exists_size = count_branch->topo.end - count_branch->topo.start;
      // offset = count_branch->topo.mid - count_branch->topo.start;
      exists_size = count_branch->parent->rho;
      offset = exists_size/2;	

      // cout << exists_size << endl;
      // cout << offset << endl;

      // iterate over regions
      for (chromo_region &c_r: c_rs)
	{

	  start_circ = (c_r.start + offset)%circ_size;
	  end_circ = (c_r.end + offset)%circ_size;

	  // cout << c_r.name << " " << start_circ << " " << end_circ << endl;

	  if ((start_circ >= 0) &&
	      (end_circ >= 0) &&
	      (start_circ < exists_size) &&
	      (end_circ < exists_size))
	    {

	      c_r.count += 1;

	    }

	}

    }
  
}


// function used to dump the regions and counts to a file
void btree::dump_regions(string rg_filename, vector<chromo_region> c_rs)
{
  fstream rg_file;

  rg_file.open(rg_filename, ios::out);

  if (!rg_file)
    {
      cout << "ERROR: file not opened in dump_regions" << endl;
    }
  else
    {
  
      for (chromo_region c_r: c_rs)
	{

	  rg_file << c_r.name
		  << "," << c_r.start
		  << "," << c_r.end
		  << ":" << c_r.count
		  << endl;
      
	}

    }

  rg_file.close();
}

void btree::dump_CG_map(string CG_filename, int idx, CG_map &m)
{
  fstream CG_file;

  CG_file.open(CG_filename, ios::out);

  if (!CG_file)
    {
      cout << "ERROR: file not opened in dump_CG_map" << endl;
    }
  else
    {

      CG_file << "N_base = "
	      << m.N_base
	      << "\nN = "
	      << m.N
	      << "\nf_CG = "
	      << m.f_CG
	      << "\nN_base_CG = "
	      << m.N_base_CG
	      << "\nN_CG = "
	      << m.N_CG
	      << endl;

      CG_file << "\nN_leaves = "
	      << m.CG_leaves.size()
	      << endl;

      for (CG_leaf leaf: m.CG_leaves)
	{

	  CG_file << leaf.leaf
		  << ","
		  << leaf.start + idx
		  << ","
		  << leaf.end + idx
		  << endl;
      
	}

      CG_file << "\nID, base-ID, min, max" << endl;
      
  
      for (CG_locus l: m.loci)
	{

	  CG_file << l.CG + idx
		  << ","
		  << l.bCG + idx
		  << ","
		  << l.start + idx
		  << ","
		  << l.end + idx
		  << endl;
      
	}

    }

  CG_file.close();
}


// function to create a coarse-graining for the entire system
CG_map btree::update_CG_map(int f_CG)
{

  CG_map m;
  CG_leaf temp_CG_leaf;

  int N_leaves = count_total_leaves();
  int init_size;//, N_CG_leaf;
  vector<string> leaves = get_leaves();

  m.f_CG = f_CG;
  m.N_base = root->size;
  m.N = total_size();

  for (int i=0; i<N_leaves; i++)
    {

      init_size = m.loci.size();

      centered_CG_map(m.loci,get_branch(leaves[i]),f_CG);

      temp_CG_leaf.leaf = leaves[i];
      temp_CG_leaf.start = init_size;
      temp_CG_leaf.end = m.loci.size() - 1;
      // N_CG_leaf = temp_CG_leaf.end - temp_CG_leaf.start + 1;

      if (i == 0)
	{
	  m.N_base_CG = m.loci.size();
	}
      else
	{

	  for (int j=temp_CG_leaf.start; j<temp_CG_leaf.end+1; j++)
	    {
	      m.loci[j].CG += (m.loci[temp_CG_leaf.start-1].CG + 1);
	      m.loci[j].start += (m.loci[temp_CG_leaf.start-1].end + 1);
	      m.loci[j].end += (m.loci[temp_CG_leaf.start-1].end + 1);
	    }
	  
	}
      
      m.CG_leaves.push_back(temp_CG_leaf);
      
    }

  m.N_CG = m.loci.size();

  return m;
}


// function to create a centered coarse-graining for a single branch
void btree::centered_CG_map(vector<CG_locus> &loci, node *branch, int f_CG)
{

  CG_locus l;
  int N = branch->topo.end - branch->topo.start + 1;
  int ori_idx_m = branch->topo.mid - branch->topo.start;
  int ori_idx_p = ori_idx_m + 1;

  int N_cum, N_CG_cum;

  int CG_rem_m, CG_rem_p, N_CG_m, N_CG_p, N_CG;

  CG_rem_m = (ori_idx_m+1)%f_CG;
  CG_rem_p = (N - ori_idx_p)%f_CG;

  N_CG_m = (ori_idx_m+1)/f_CG;
  if (CG_rem_m > 0) N_CG_m++;
  N_CG_p = (N - ori_idx_p)/f_CG;
  if (CG_rem_p > 0) N_CG_p++;

  N_CG = N_CG_m + N_CG_p;

  N_cum = 0;
  N_CG_cum = 0;

  // cout << "N=" << N << endl;
  // cout << "ori_idx_m=" <<  ori_idx_m << endl;
  // cout << "CG_rem_m=" <<  CG_rem_m << endl;
  // cout << "ori_idx_p=" <<  ori_idx_p << endl;
  // cout << "CG_rem_p=" <<  CG_rem_p << endl;

  while (N_CG_cum < N_CG)
    {

      // conditional for remainder beads at start
      if ((N_CG_cum == 0) &&
	  (CG_rem_m > 0))
	{

	  l.CG = N_CG_cum;
	  l.bCG = N_CG_cum;
	  l.start = N_cum;
	  l.end = N_cum + CG_rem_m - 1;

	  N_CG_cum += 1;
	  N_cum += CG_rem_m;

	  loci.push_back(l);

	  continue;
	  
	}

      // conditional for remainder beads at end
      if ((N_CG_cum == N_CG - 1) &&
	  (CG_rem_p > 0))
	{

	  l.CG = N_CG_cum;
	  l.bCG = N_CG_cum;
	  l.start = N_cum;
	  l.end = N_cum + CG_rem_p - 1;

	  N_CG_cum += 1;
	  N_cum += CG_rem_p;

	  loci.push_back(l);

	  continue;
	  
	}

      l.CG = N_CG_cum;
      l.bCG = N_CG_cum;
      l.start = N_cum;
      l.end = N_cum + f_CG - 1;

      N_CG_cum += 1;
      N_cum += f_CG;

      loci.push_back(l);
      

    }
  
}


// function used by destructor
void btree::destroy_tree(node *branch)
{
  if (branch->leaf == false)
    {
      destroy_tree(branch->left);
      destroy_tree(branch->right);
    }
  delete branch;
}

void btree::destroy_tree()
{
  if (root != nullptr)
    {
      destroy_tree(root);
      reset_root();
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
      if (branch->parent != nullptr)
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


  cout << gen_offset
       << "start = "
       << branch->topo.start
       << ", mid = "
       << branch->topo.mid
       << ", end = "
       << branch->topo.end
       << endl;
  cout << gen_offset
       << "start_link = "
       << branch->topo.start_link
       << ", end_link = "
       << branch->topo.end_link
       << endl;

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
  
  if (root != nullptr)
    {

      int N_leaves = count_total_leaves();
      int N_forks = count_total_forks();
      int N_completed_forks = count_completed_forks();
      int N_active_forks = count_active_forks();
      
      cout << "\nprinting tree with "
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

      cout << "total_size = "
	   << total_size()
	   << endl;
      
      print_branch(root);
    }
  else
    {
      cout << "\n\ntree does not exist" << endl;
    }
  cout << "\n\n" << endl;
}

void btree::foo()
{
  cout << "btree\n";
}
