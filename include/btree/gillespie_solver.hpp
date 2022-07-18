#ifndef INCLUDE_GILLESPIE_SOLVER_HPP
#define INCLUDE_GILLESPIE_SOLVER_HPP

#include <btree/btree.hpp>

using namespace std;


class gillespie_solver
{
public:

  // constructor and destructor
  gillespie_solver();
  ~gillespie_solver();

  // set the prng seed
  void prng_seed(int s);
  

private:

  void update_propensities();
  void update_state();

  mt19937 rand_eng;

};

#endif
