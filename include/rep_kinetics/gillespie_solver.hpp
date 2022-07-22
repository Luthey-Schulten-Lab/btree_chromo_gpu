#ifndef INCLUDE_GILLESPIE_SOLVER_HPP
#define INCLUDE_GILLESPIE_SOLVER_HPP

#include <algorithm>
#include <vector>
#include <memory>
#include <random>
#include <cmath>

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
  void run_FPT(double t, double t_max);

  
private:

  // functions for Gillespie loop
  int select_rxn(double r_rxn, double total_propensity);
  void update_propensities();
  void update_state(int j_rxn);

  // private variables
  mt19937 rand_eng;
  uniform_real_distribution<double> u_rand;

};

#endif
