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

#include <rxn_manipulator.hpp>

class replication_model;

struct rep_model_params
{
  double k_rep;
  double k_c, k_d;
  double k_hi, k_lo, k_on, k_off, k_bubble;
  int N_hi, N_lo, N_fil;
};

struct init_loc
{
  std::string loc;
  int N;
};

class replication_model
{
  
public:

  // constructor and destructor
  replication_model();
  ~replication_model();

  // read the replication model
  void load_model(std::string rep_model_filename);

  // setter for number of leaves
  void set_N_leaves(int n);

  // getters to access replication model details
  double get_k_rep();
  int get_init_requirement();
  int get_N_species();
  int get_N_per_leaf();
  int get_N_non_leaf();
  int get_M_rxns();

  // create the FPT species counts
  std::vector<species_count> create_xFPT();

  // get the reactions
  std::vector<reaction> get_reactions();
  
private:

  // functions for creation of reaction model
  void number_rep_species();
  void number_rep_rxns();

  rep_model_params r_m_p;
  int N_species, M_rxns;
  int N_leaves;
  int N_non_leaf, N_per_leaf;
  int M_non_leaf, M_per_leaf;

  rxn_manipulator rxn_manip;
  
};

#endif
