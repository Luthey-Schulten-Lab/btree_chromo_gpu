#ifndef INCLUDE_REPLICATOR_HPP
#define INCLUDE_REPLICATOR_HPP

#include <replication_model.hpp>
#include <gillespie_solver.hpp>

class replicator
{
public:

  // constructor and destructor
  replicator();
  ~replicator();

  // set the prng seed
  void prng_seed(int s);

  // read the replication model
  void read_rep_model(std::string rep_model_filename);

  // reset the noninitiator species
  void reset_noninit_s(int * &noninit_s);
  
  // prepare the reaction system
  void prepare_system(int *noninit_s, std::vector<init_loc> &init_dist);

  // simulate the system until t_max
  void run_replicate_FPT(int *noninit_s,
			 std::vector<init_loc> &init_dist,
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
