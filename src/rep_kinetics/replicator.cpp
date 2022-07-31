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
  // solver.set_propensity_fxn(propensities);

  // print the system as a sanity check
  solver.print_reaction_system();
  
}



void replicator::run_replicate_FPT(vector<init_loc> &init_dist, double &t, double &t_max)
{

  // prepare the reaction system based on the reaction model
  prepare_system(init_dist);

  // run the system

  // redistribute initiators based on FPT result
  rep_model.update_id_from_sc(init_dist,solver.state_to_sc());

  for (init_loc temp_i_l : init_dist)
    {
      cout << temp_i_l.loc << " = " << temp_i_l.N << endl;
    }
  
  t += t_max/2.0;
  init_dist.back().N = -1;

  solver.destroy_reaction_system();
}
