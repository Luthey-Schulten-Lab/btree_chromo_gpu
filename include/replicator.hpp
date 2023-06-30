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
  // set the DnaA gene copy number
  void set_max_DnaA_genes(int max_DnaA_genes);
  // set the maximum number of replisomes
  void set_replisomes(int n_replisomes);
  // set the volume or load a volume protocol
  void set_volume(double V);
  void load_volume_protocol(std::string volume_protocol_filename);

  // get the time
  double get_time();
  // get the free DnaA
  int get_DnaA();
  int partition_DnaA();
  // get the replication rate
  double get_k_rep();
  // get the replisomes
  int get_replisomes();
  int get_free_replisomes();
  int partition_free_replisomes();

  // get the initiator distribution
  std::vector<init_loc> get_init_dist();
  // reset the initiator distribution
  void reset_init_dist();
  // update the initiator distribution
  void update_init_dist(std::vector<std::string> leaves);
  // add initiators on leaf
  void bind_init(std::string leaf, int n);
  // test for initiation event and update
  std::string initiation_test();
  
  // get the set of replicating forks
  std::vector<std::string> get_replicating_forks();
  // unbind the replisomes from completed forks
  void unbind_replisomes(std::vector<std::string> completed_forks);
  // reset the replisomes
  void reset_replisomes();
  // bind replisome on fork
  void bind_replisome(std::string fork);

  // read the replication model
  void load_model(std::string rep_model_filename);

  // run the replicator
  void run(double &dt, double dt_target);

  // dump the relicator state
  void dump_state(std::string state_filename);
  // load the relicator state
  void load_state(std::string state_filename);
  
private:

  // print the initiator distribution
  void print_state();

  // create the vector of species counts
  std::vector<species_count> prepare_species_counts();
  // create the FPT vector
  std::vector<species_count> prepare_xFPT();

  // convert species counts to the replicator state
  void species_counts_to_replicator_state(std::vector<species_count> s_cs);

  bool V_protocol;

  int N_leaves;
  int n_DnaA, n_replisomes, n_free_replisomes;
  int n_DnaA_genes, max_DnaA_genes;
  
  double t;
  double V;
  std::vector<volume_point> V_points;

  std::vector<init_loc> init_dist;
  std::vector<std::string> replicating_forks;

  gillespie_solver solver;
  replication_model rep_model;

  std::mt19937 rand_eng;

};

#endif
