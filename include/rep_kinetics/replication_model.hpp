#ifndef INCLUDE_REPLICATION_MODEL_HPP
#define INCLUDE_REPLICATION_MODEL_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>
#include <memory>
#include <cmath>
#include <iterator>

#include <rep_kinetics/rxn_manipulator.hpp>

using namespace std;

struct rep_model_params
{
  double V;
  double k_rep;
  double k_c, k_d;
  double k_hi, k_lo, k_on, k_off;
  int N_hi, N_lo, N_fil;
  int max_replisomes, N_init_DnaA;
};

struct init_loc
{
  string loc;
  int N;
};

class replication_model
{
public:

  // constructor and destructor
  replication_model();
  ~replication_model();

  // read the replication model
  void read_rep_model(string rep_model_filename);

  // update the propensities
  void propensities(int *x, double *W);

  // getters to access state of replication model
  int get_N_init_DnaA();
  int get_max_replisomes();
  double get_k_rep();
  int get_N_species();
  int get_M_rxns();

  // setter for number of leaves
  void set_N_leaves(int n);

  // convert initiator distribution to species counts and vice-versa
  vector<species_count> id_to_sc(vector<init_loc> &init_dist);
  void update_id_from_sc(vector<init_loc> &init_dist, vector<species_count> s_cs);

  // create the FPT species counts
  vector<species_count> create_xFPT();

  // get the reactions
  vector<reaction> get_reactions();

  // functions for DnaA reaction model
  void number_rep_species();
  void number_rep_rxns();
  
private:

  rep_model_params r_m_p;
  int N_species, M_rxns;
  int N_leaves;
  int N_per_leaf;

  rxn_manipulator rxn_manip;
  
};

#endif
