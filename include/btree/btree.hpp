#ifndef INCLUDE_BTREE_HPP
#define INCLUDE_BTREE_HPP

#include <iostream>
#include <string>
#include <cstring>
#include <algorithm>
#include <vector>

using namespace std;

struct node
{
  int gen;
  int size;
  int rho;
  bool leaf;
  node *parent;
  node *left;
  node *right;
};

struct fork_rho
{
  string fork;
  int rho;
};

struct btree_state
{
  int size;
  vector<fork_rho> fork_rhos;
};

struct btree_transforms
{
  vector<fork_rho> fork_rhos;
};

class btree
{
public:

  // constructor and destructor
  btree();
  ~btree();

  // recursively destroys tree
  void destroy_tree();

  // prints entire tree and information
  void print_tree();

  // initialize, branch, and grow at branches
  void initialize_tree(int s);
  int branch(string loc);
  int grow_at_branch(string loc, int r);

  // prepare and dump state
  void prepare_state(btree_state st);
  btree_state dump_state();

  // apply transformations to btree
  void apply_transforms(btree_transforms tr);

  // queries about tree state
  int count_total_leaves();
  int count_total_forks();
  int count_completed_forks();
  int count_active_forks();
  int count_leaves(string loc);
  int count_total_forks(string loc);
  int count_completed_forks(string loc);
  int count_active_forks(string loc);

  // get details of tree
  vector<string> get_completed_forks();
  vector<string> get_active_forks();
  vector<string> get_leaves();
  
  void foo();

private:

  // recursively destroy tree
  void destroy_tree(node *branch);

  // recursively print branch topology
  void print_branch(node *branch);

  // used for calculating tree state
  int leaf_counter(node *branch);
  int total_fork_counter(node *branch);
  int completed_fork_counter(node *branch);
  int active_fork_counter(node *branch);

  // traverse leaves and forks to determine identities
  vector<string> traverse_leaves(vector<string> leaves, node *branch);
  vector<string> traverse_completed_forks(vector<string> forks, int i, node *branch);
  vector<string> traverse_active_forks(vector<string> forks, string f, node *branch);

  // branch manipulation routines
  node *get_branch(string loc);
  void initialize_branch(node *branch, int g, int s);
  void split_branch(node *branch);
  node *parse_dir(node *branch, char d);

  node *root;
  

};

#endif
