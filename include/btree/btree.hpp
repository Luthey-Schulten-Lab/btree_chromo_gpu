#ifndef INCLUDE_BTREE_HPP
#define INCLUDE_BTREE_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>
#include <vector>
#include <memory>
#include <random>

using namespace std;

struct theta_topo
{
  int start, end, start_link, end_link, mid;
};

struct node
{
  int gen;
  int size;
  int rho_t, rho_cw, rho_ccw;
  bool leaf, complete;
  theta_topo topo;
  node *parent;
  node *left;
  node *right;
};

struct fork_rho
{
  string fork;
  int rho_cw, rho_ccw;
};

typedef vector<fork_rho> btree_transforms;

struct btree_state
{
  int size;
  //vector<fork_rho> fork_rhos;
  btree_transforms transforms;
};

struct chromo_region
{
  string name;
  int start, end, count;
};

struct CG_locus
{
  int CG, bCG, start, end;
};

struct CG_leaf
{
  string leaf;
  int start, end, mid;
};

struct CG_map
{
  int N_base, N;
  int N_base_CG, N_CG;
  int f_CG;

  vector<CG_locus> loci;
  vector<CG_leaf> CG_leaves;
};

struct mono_range
{
  bool wrapped;
  int N, ll, ul, mid_ll, mid_ul;
};

struct fork_partition
{
  string fork;
  vector<mono_range> left_monos, right_monos;
};

class btree
{
public:

  // constructor and destructor
  btree();
  ~btree();

  // reset root
  void reset_root();

  // recursively destroys tree
  void destroy_tree();

  // prng seeding
  void prng_seed(int s);

  // prints entire tree and information
  void print_tree();

  // initialize and branch
  void initialize_tree(int s);
  int branch(string loc);

  // prepare and dump state
  void prepare_state(btree_state st);
  btree_state dump_state();

  // read and write the state
  void write_state(string st_filename, btree_state st);
  btree_state read_state(string st_filename);
  btree_transforms read_transforms(string tr_filename);

  // parse transform of form "(branch)_cw(rho_cw)_ccw(rho_ccw)"
  fork_rho parse_transform(string s);

  // apply transformations to btree
  void apply_transforms(btree_transforms tr);
  void single_transform(fork_rho f_r);
  void random_transforms(int r);

  // queries about tree state
  int count_total_leaves();
  int count_total_forks();
  int count_completed_forks();
  int count_active_forks();
  int total_size();
  int max_size();

  // solve the theta structure topology
  void solve_topology();
  void dump_topology(string topo_filename, int idx);
  theta_topo get_leaf_topo(string loc);
  void dump_fork_partitions(string fork_partitions_filename, int idx);

  // update the coarse-graining map based on the current state
  CG_map update_CG_map(int f_CG);
  void dump_CG_map(string CG_filename, int idx, CG_map &m);

  // read, update, and dump chromo_regions
  vector<chromo_region> read_regions(string rg_filename, int idx);
  void update_region_counts(vector<chromo_region> &c_rs);
  void dump_regions(string rg_filename, vector<chromo_region> c_rs);

  // get details of tree
  vector<string> get_completed_forks();
  vector<string> get_active_forks();
  vector<string> get_leaves();

  // create the topologies
  void prepare_bonds(int **&c, int *&t, int &N, int idx);
  void prepare_angles(int **&c, int *&t, int &N, int idx);
  void prepare_types(int *&t, int &N, int base_type);
  
  void foo();

private:

  // recursively destroy tree
  void destroy_tree(node *branch);

  // recursively print branch topology
  void print_branch(node *branch);

  // used for calculating tree state
  int count_leaves(string loc);
  int count_total_forks(string loc);
  int count_completed_forks(string loc);
  int count_active_forks(string loc);
  int leaf_counter(node *branch);
  int total_fork_counter(node *branch);
  int completed_fork_counter(node *branch);
  int active_fork_counter(node *branch);
  int branch_size(node *branch);

  // used for calculating growth
  int grow_at_branch_sym(string loc, int r);
  int grow_at_branch_asym(string loc, int r_cw, int r_ccw);
  int get_max_growth_cw(node *branch);
  int get_max_growth_ccw(node *branch);
  array<int,2> partition_growths_sym(node *branch, int proposed_r_cw, int proposed_r_ccw);

  // get the fork partitions about a fork
  vector<fork_partition> get_all_fork_partitions();
  fork_partition get_fork_partition(string loc);

  // create centered CG maps per branch
  void centered_CG_map(int &mid, vector<CG_locus> &loci, node *branch, int f_CG);

  // traverse leaves and forks to determine identities
  vector<string> traverse_leaves(vector<string> leaves, node *branch);
  vector<string> traverse_completed_forks(vector<string> forks, int i, node *branch);
  vector<string> traverse_active_forks(vector<string> forks, string f, node *branch);

  // branch manipulation routines
  node *get_branch(string loc);
  void initialize_branch(node *branch, int g, int s);
  void split_branch(node *branch);
  node *parse_dir(node *branch, char d);


  // uniform_real_distribution<> u_dist;
  mt19937 rand_eng;
  node *root;
  

};

#endif
