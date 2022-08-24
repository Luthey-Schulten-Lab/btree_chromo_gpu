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
  this->solver.prng_seed(s);
}


// load replication model parameters
void replicator::read_rep_model(string rep_model_filename)
{
  this->rep_model.read_rep_model(rep_model_filename);
}


// get the initial number of DnaA
int replicator::get_N_init_DnaA()
{
  return this->rep_model.get_N_init_DnaA();
}


// get the maximum number of replisomes
int replicator::get_max_replisomes()
{
  return this->rep_model.get_max_replisomes();
}


// get the replication rate
double replicator::get_k_rep()
{
  return this->rep_model.get_k_rep();
}


void replicator::reset_noninit_s(int * &noninit_s)
{
  this->rep_model.reset_noninit_s(noninit_s);
}


// prepare the system state using the replication model
void replicator::prepare_system(int *noninit_s,
				vector<init_loc> &init_dist)
{

  // set the number of leaves in the replication model
  this->rep_model.set_N_leaves(static_cast<int>(init_dist.size())-1);

  // set the number of replication species and reactions in the replication model
  this->rep_model.number_rep_species();
  this->rep_model.number_rep_rxns();

  // initialize the solver with N and M from rep_model
  this->solver.initialize_reaction_system(rep_model.get_N_species(),
					  rep_model.get_M_rxns());

  // get reaction stoichiometries from the replication model
  vector<reaction> rxns = this->rep_model.get_reactions();

  // set the reaction stoichiometries in the solver
  this->solver.set_S(rxns);

  // convert the initiator distribution to species counts
  vector<species_count> init_s_cs = this->rep_model.id_to_sc(init_dist);
  // set the state vector values from the initiator distribution
  this->solver.set_x(init_s_cs);

  // update gene counts based on the number of leaves
  this->rep_model.update_noninit_s(noninit_s);
  // set the state vector values from the noninitiator species
  vector<species_count> noninit_s_cs;
  species_count s_c;
  for (int i=0; i<(this->rep_model.get_N_non_leaf()-1); i++)
    {
      s_c.id = i;
      s_c.N = noninit_s[i];
      noninit_s_cs.push_back(s_c);
    }
  this->solver.set_x(noninit_s_cs);

  // set the FPT state vector
  this->solver.set_xFPT(rep_model.create_xFPT());

  // set the propensity function
  this->solver.set_replication_model(rep_model);

  // print the system as a sanity check
  // this->solver.print_reaction_system();
  
}



void replicator::run_replicate_FPT(int *noninit_s,
				   vector<init_loc> &init_dist,
				   double &t, double &t_max)
{

  // prepare the reaction system based on the reaction model
  prepare_system(noninit_s,
		 init_dist);

  // run the system
  this->solver.run_FPT(t,t_max);

  // redistribute initiators based on FPT result
  vector<species_count> solver_s_cs = this->solver.state_to_sc();
  this->rep_model.update_init_from_solver_s_cs(init_dist,solver_s_cs);
  this->rep_model.update_noninit_s_from_solver_s_cs(noninit_s,solver_s_cs);

  int fil_trigger = this->rep_model.get_N_per_leaf() - 1;
  for (size_t i=1; i<init_dist.size(); i++)
    {
      if (init_dist[i].N == fil_trigger)
	{
	  init_dist[i].N = -1;
	  init_dist[0].N += fil_trigger;
	  break;
	}
    }

  this->solver.destroy_reaction_system();
}
