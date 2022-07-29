#ifndef INCLUDE_REPLICATION_MODEL_HPP
#define INCLUDE_REPLICATION_MODEL_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>
#include <vector>
#include <memory>
#include <cmath>
#include <iterator>

#include <rep_kinetics/gillespie_solver.hpp>

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
  rep_model_params read_rep_model(string rep_model_filename);

  // update the propensities
  static void propensities(int *x, double *W);

  // simulate the system until t_max
  void run_replicate_FPT(rep_model_params &rep_model, vector<init_loc> &init_dist, double &t, double &t_max);
  
private:

  // convert initiator distribution to species counts and vice-versa
  vector<species_count> id_to_sc(vector<init_loc> &init_dist);
  void update_id_from_sc(vector<init_loc> &init_dist, vector<species_count> s_cs);

  // create the FPT species counts
  vector<species_count> create_xFPT();

  // prepare the system
  void prepare_system(rep_model_params rep_model, int n);

  // get the reactions
  void get_reactions(vector<reaction> &rxns);
  void reset_reaction(reaction *r);

  // functions for DnaA reaction model
  void number_rep_species();
  void number_rep_rxns();

  gillespie_solver solver;
  rep_model_params r_m_p;
  int N_species, M_rxns, N_leaves;
  int N_per_leaf;
  
};

#endif
