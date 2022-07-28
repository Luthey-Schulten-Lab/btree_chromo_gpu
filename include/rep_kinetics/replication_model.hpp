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

  // simulate the system until t_max
  void run_replicate_FPT(rep_model_params r_m_p, vector<init_loc> &init_dist, double &t, double &t_max);
  
private:

  // prepare the system
  void prepare_system(rep_model_params r_m_p, int N_leaves);

  // get the reactions
  void get_reactions(rep_model_params r_m_p, int N_leaves, vector<reaction> &rxns);
  void reset_reaction(reaction *r);

  // functions for DnaA reaction model
  int number_rep_species(rep_model_params r_m_p, int N_leaves);
  int number_rep_rxns(rep_model_params r_m_p, int N_leaves);

  gillespie_solver solver;
  int N_species, M_rxns, N_leaves;
  
};

#endif
