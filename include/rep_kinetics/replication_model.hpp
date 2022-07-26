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
  void run_replicate_FPT(rep_model_params r_m_p, int N_forks, double t_max);
  
private:

  // prepare the system
  void prepare_system(rep_model_params r_m_p, int N_forks);

  // functions for DnaA reaction model
  int number_rep_species(rep_model_params r_m_p, int N_forks);
  int number_rep_rxns(rep_model_params r_m_p, int N_forks);

  gillespie_solver solver;
  int N_species, M_rxns, N_forks;
  
};

#endif
