#ifndef INCLUDE_REPLICATOR_HPP
#define INCLUDE_REPLICATOR_HPP

#include <rep_kinetics/replication_model.hpp>
#include <rep_kinetics/gillespie_solver.hpp>

using namespace std;


class replicator
{
public:

  // constructor and destructor
  replicator();
  ~replicator();

  // set the prng seed
  void prng_seed(int s);

  // read the replication model
  void read_rep_model(string rep_model_filename);

  // reset the noninitiator species
  vector<species_count> get_reset_noninit();
  
  // prepare the reaction system
  void prepare_system(vector<species_count> &noninit_s_cs, vector<init_loc> &init_dist);

  // simulate the system until t_max
  void run_replicate_FPT(vector<species_count> &noninit_s_cs,
			 vector<init_loc> &init_dist,
			 double &t, double &t_max);

  // getters to access replication model
  int get_N_init_DnaA();
  int get_max_replisomes();
  double get_k_rep();
  
private:

  gillespie_solver solver;
  replication_model rep_model;

};

#endif
