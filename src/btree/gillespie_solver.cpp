#include <btree/gillespie_solver.hpp>

// constructor
gillespie_solver::gillespie_solver()
{
  rand_eng.seed(0);
}

// destructor
gillespie_solver::~gillespie_solver()
{

}

// prng seed
void gillespie_solver::prng_seed(int s)
{
  rand_eng.seed(s);
}


