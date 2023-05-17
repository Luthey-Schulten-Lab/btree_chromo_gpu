#ifndef INCLUDE_GILLESPIE_SOLVER_HPP
#define INCLUDE_GILLESPIE_SOLVER_HPP

#include <algorithm>
#include <memory>
#include <random>
#include <cmath>
#include <iostream>

#include <rxn_manipulator.hpp>
#include <replication_model.hpp>

using namespace std;

class gillespie_solver
{
public:

  // constructor and destructor
  gillespie_solver();
  ~gillespie_solver();

  // set the prng seed
  void prng_seed(int s);

  // run the system until max time or first-passage occurs
  void run_FPT(double &t, double &t_max);

  // initialize and destroy the reaction system
  void initialize_reaction_system(int N, int M);
  void destroy_reaction_system();

  // print the reaction system
  void print_reaction_system();
  
  // set arrays
  void set_x(vector<species_count> &s_cs);
  void set_xFPT(vector<species_count> s_cs);
  void set_S(vector<reaction> &rxns);
  void set_replication_model(replication_model &r_m);

  // convert the state to species counts
  vector<species_count> state_to_sc();

  void update_propensities();
  
private:

  // set N and M
  void set_N(int N);
  void set_M(int M);

  // initialize and destroy state vector
  void initialize_x();
  void destroy_x();

  // initialize and destroy FPT state vector
  void initialize_xFPT();
  void destroy_xFPT();

  // initialize and destroy the rate vector
  void initialize_W();
  void destroy_W();
  
  // initialize and destroy stoichiometry matrix
  void initialize_S();
  void destroy_S();

  // functions for Gillespie loop
  int select_rxn(double r_rxn, double total_propensity);
  double calc_total_propensity();
  void update_state(int j_rxn);
  int test_FPT();

  // private variables
  int N, M; // N species, M reactions
  int *x, *xFPT, **S; // state vector x, FPT state vector xFPT, and stoichiometry matrix S
  double *W; // rate vector
  mt19937 rand_eng;
  uniform_real_distribution<double> u_rand;
  replication_model rep_model;

};

#endif
