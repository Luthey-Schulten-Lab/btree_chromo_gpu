#include <rep_kinetics/replicator.hpp>

// constructor
replicator::replicator()
{
}


// destructor
replicator::~replicator()
{
}


// prng seed
void replicator::prng_seed(int s)
{
  solver.prng_seed(s);
}


// load replication model parameters
void replicator::read_rep_model(string rep_model_filename)
{
  rep_model.read_rep_model(rep_model_filename);
}


// get the initial number of DnaA
int replicator::get_N_init_DnaA()
{
  return rep_model.get_N_init_DnaA();
}


// get the maximum number of replisomes
int replicator::get_max_replisomes()
{
  return rep_model.get_max_replisomes();
}


// get the replication rate
double replicator::get_k_rep()
{
  return rep_model.get_k_rep();
}


// prepare the system state using the replication model
void replicator::prepare_system(vector<init_loc> &init_dist)
{

  // set tje number of leaves in the replication model
  rep_model.set_N_leaves(static_cast<int>(init_dist.size())-1);

  // set the number of replication species and reactions in the replication model
  rep_model.number_rep_species();
  rep_model.number_rep_rxns();

  // initialize the solver with N and M from rep_model
  solver.initialize_reaction_system(rep_model.get_N_species(),
				    rep_model.get_M_rxns());

  // get reaction stoichiometries from the replication model
  vector<reaction> rxns = rep_model.get_reactions();

  // set the reaction stoichiometries in the solver
  solver.set_S(rxns);

  // convert the initiator distribution to species counts
  vector<species_count> s_cs = rep_model.id_to_sc(init_dist);
  
  // set the state vector
  solver.set_x(s_cs);

  // set the FPT state vector
  solver.set_xFPT(rep_model.create_xFPT());

  // set the propensity function
  solver.set_replication_model(rep_model);

  // print the system as a sanity check
  // solver.print_reaction_system();
  
}



void replicator::run_replicate_FPT(vector<init_loc> &init_dist, double &t, double &t_max)
{

  // prepare the reaction system based on the reaction model
  prepare_system(init_dist);

  // run the system
  solver.run_FPT(t,t_max);

  // redistribute initiators based on FPT result
  rep_model.update_id_from_sc(init_dist,solver.state_to_sc());

  int fil_trigger = rep_model.get_N_per_leaf() - 1;

  for (size_t i=1; i<init_dist.size(); i++)
    {
      if (init_dist[i].N == fil_trigger)
	{
	  init_dist[i].N = -1;
	  init_dist[0].N += fil_trigger;
	  break;
	}
    }

  solver.destroy_reaction_system();
}
