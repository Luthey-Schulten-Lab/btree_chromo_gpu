#ifndef INCLUDE_GILLESPIE_SOLVER_HPP
#define INCLUDE_GILLESPIE_SOLVER_HPP

#include <algorithm>
#include <memory>
#include <random>
#include <cmath>
#include <iostream>

//#include <boost/multi_array.hpp>

#include <rxn_manipulator.hpp>
#include <replication_model.hpp>

class gillespie_solver
{
public:

  // constructor and destructor
  gillespie_solver();
  ~gillespie_solver();

  // set the prng seed
  void prng_seed(int s);

  // run the system until max time or first-passage occurs
  void run_FPT(double t, double &dt, double &dt_max);

  // initialize and destroy the reaction system
  void initialize_reaction_system(int N, int M);

  // print the reaction system
  void print_reaction_system();

  // set the volume
  void set_V(double V);
  
  // set arrays
  void set_x(std::vector<species_count> &s_cs);
  void set_xFPT(std::vector<species_count> s_cs);
  void set_reaction_system(std::vector<reaction> &rxns);

  // convert the state to species counts
  std::vector<species_count> state_to_sc();
  
private:

  void destroy_reaction_system();

  // functions for Gillespie loop
  int select_rxn(double r_rxn, double total_propensity);
  double calc_total_propensity();
  void update_propensities();
  void update_state(int j_rxn);
  int test_FPT();

  // private variables

  const int max_nproducts = 4;
  double V;
  int N, M; // N species, M reactions
  int *x; // particle counts
  int *nreactants, *nproducts; // number reactants and products per reaction
  int **reactants, **products; // reactants and products per reaction
  double *k, *W; // rates and propensities

  std::vector<species_count> xFPT;
  
  std::mt19937 rand_eng;
  std::uniform_real_distribution<double> u_rand;

};

#endif
