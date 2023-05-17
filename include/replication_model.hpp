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

class replication_model;

// typedef void (replication_model::* ptr)(int *, double*);
// typedef ptr (*pm)();

struct rep_model_params
{
  double V;
  double k_rep;
  double k_c, k_d;
  double k_hi, k_lo, k_on, k_off;
  double tau_SA, k_SA;
  int N_hi, N_lo, N_fil;
  int max_replisomes;
  int N_init_DnaA, N_init_SA;
};

struct init_loc
{
  string loc;
  int N;
};

class replication_model
{
  
public:

  // typedef void (replication_model::* ptr)(int *, double*);
  // typedef ptr (*pm)();

  // constructor and destructor
  replication_model();
  ~replication_model();

  // read the replication model
  void read_rep_model(string rep_model_filename);

  // update the propensities
  void propensities(int *x, double *W);
  // ptr get_propensity_fxn();

  // getters to access state of replication model
  int get_N_init_DnaA();
  int get_N_leaves();
  int get_max_replisomes();
  double get_k_rep();
  int get_N_species();
  int get_N_per_leaf();
  int get_N_non_leaf();
  int get_M_rxns();

  // setter for number of leaves
  void set_N_leaves(int n);

  // return reset species counts of noninitiator species
  void reset_noninit_s(int * &noninit_s);
  void update_noninit_s(int *noninit_s);
  void update_noninit_s_from_solver_s_cs(int *noninit_s,
					 vector<species_count> &solver_s_cs);

  // convert initiator distribution to species counts and vice-versa
  vector<species_count> id_to_sc(vector<init_loc> &init_dist);
  void update_init_from_solver_s_cs(vector<init_loc> &init_dist,
				    vector<species_count> &solver_s_cs);

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
  int N_non_leaf, N_per_leaf;

  rxn_manipulator rxn_manip;
  
};

#endif
