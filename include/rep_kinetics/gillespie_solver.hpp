#ifndef INCLUDE_GILLESPIE_SOLVER_HPP
#define INCLUDE_GILLESPIE_SOLVER_HPP

#include <algorithm>
#include <vector>
#include <memory>
#include <random>
#include <cmath>
#include <iostream>

using namespace std;

struct species_count
{
  int id;
  int n;
};

struct reaction
{
  vector<species_count> inputs, outputs;
};

class gillespie_solver
{
public:

  // constructor and destructor
  gillespie_solver();
  ~gillespie_solver();

  // set the prng seed
  void prng_seed(int s);

  // run the system until max time or first-passage occurs
  void run_FPT(double t, double t_max);

  // prepare the reaction system
  void prepare_reaction_system(int N_species, int N_rxns);

  // print the reaction system
  void print_reaction_system();

  void reset_reaction(reaction &r);
  void add_reaction_input(reaction &r, int id, int n);
  void add_reaction_output(reaction &r, int id, int n);
  
  // set arrays
  void set_x(vector<species_count> &s_cs);
  void set_xFPT(vector<species_count> s_cs);
  void set_S(vector<reaction> &rxns);
  void set_propensity_fxn(void (*func)(int *, double *));

  // conver the state to species counts
  vector<species_count> state_to_sc();
  
private:

  // set N and M
  void set_N(int N_species);
  void set_M(int M_rxns);

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
  void update_propensities();
  void update_state(int j_rxn);
  int test_FPT();

  // private variables
  int N, M; // N species, M reactions
  int *x, *xFPT, **S; // state vector x, FPT state vector xFPT, and stoichiometry matrix S
  double *W; // rate vector
  mt19937 rand_eng;
  uniform_real_distribution<double> u_rand;
  void (*rate_func)(int *, double *);

};

#endif
