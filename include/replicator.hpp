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

  // set the time
  void set_time(double t);
  // set the free DnaA number
  void set_DnaA(int n_DnaA);
  // set the DnaA gene copy number
  void set_DnaA_genes(int n_DnaA_genes);
  // set the maximum number of replisomes
  void set_max_replisomes(int max_rep);
  // set the volume or load a volume protocol
  void set_volume(double V);
  void load_volume_protocol(std::string volume_protocol_filename);

  // read the replication model
  void load_model(std::string rep_model_filename);

  // run the replicator
  void run(std::vector<init_loc> &init_dist,
	   double &dt, double &dt_target);
  
private:

  bool V_protocol;

  int n_DnaA, max_rep, n_DnaA_genes;
  
  double t;
  double V;

  gillespie_solver solver;
  replication_model rep_model;

};

#endif
