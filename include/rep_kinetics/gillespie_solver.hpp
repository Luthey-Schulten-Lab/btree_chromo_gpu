#ifndef INCLUDE_GILLESPIE_SOLVER_HPP
#define INCLUDE_GILLESPIE_SOLVER_HPP

#include <btree/btree.hpp>

#include <math>

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

class gillespie_solver
{
public:

  // constructor and destructor
  gillespie_solver();
  ~gillespie_solver();

  // set the prng seed
  void prng_seed(int s);

  // read the replication model
  rep_model_params read_rep_model(string rep_model_filename);

  // run the system until max time or first-passage occurs
  void run_rep_FPT(rep_model_params r_m_p, int N_forks, double t_f);

  
private:

  // prepare the system
  void prepare_system(rep_model_params r_m_p, int N_forks);

  // functions for Gillespie loop
  int select_rxn(double r_rxn, double total_propensity, );
  void update_propensities(rep_model_params r_m_p);
  void update_state(int j_rxn);

  // functions for DnaA reaction model
  int number_rep_species(rep_model_params r_m_p, int N_forks);
  int number_rep_rxns(rep_model_params r_m_p, int N_forks);

  // private variables
  mt19937 rand_eng;
  uniform_real_distribution<double> u_rand;

};

#endif
