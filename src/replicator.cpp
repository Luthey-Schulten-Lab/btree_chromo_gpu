#include <replicator.hpp>

// constructor
replicator::replicator()
{
  V_protocol = false;
  n_DnaA = 0;
  t = 0.0;
  V = 1.0; 
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


// set the time
void replicator::set_time(double t)
{
  this->t = t;
}


// set the DnaA count
void replicator::set_DnaA(int n_DnaA)
{
  this->n_DnaA = n_DnaA;
}


// set the DnaA gene copy number
void replicator::set_DnaA_genes(int n_DnaA_genes)
{
  this->n_DnaA_genes = n_DnaA_genes;
}


// set the maximum number of replisomes
void replicator::set_max_replisomes(int max_rep)
{
  this->max_rep = max_rep;
}


// set the volume
void replicator::set_volume(double V)
{
  this->V = V;
  V_protocol = false;
}


//load the replication model
void replicator::load_volume_protocol(std::string volume_protocol_filename)
{

  std::cout << volume_protocol_filename << std::endl;
  // rep_model.load_model(rep_model_filename);
}


//load the replication model
void replicator::load_model(std::string rep_model_filename)
{

  std::cout << rep_model_filename << std::endl;
  // rep_model.load_model(rep_model_filename);
}


//run the replicator for a time increment or until a first-passage
void replicator::run(std::vector<init_loc> &init_dist,
		     double &dt, double &dt_target)
{
  
  std::cout << "dt_target = " << dt_target << std::endl;

  // determine the number of Oris
  int N_leaves = static_cast<int>(init_dist.size()-1);
  // set the number of DnaA genes equal to the number of Oris
  set_DnaA_genes(N_leaves);
  // set the number of Oris for the reaction model
  rep_model.set_N_leaves(N_leaves);

  
  // prepare a vector of the initial species counts with the rep model

  // prepare a vector of FPT events

  // prepare the reaction system

  
  // initialize the solver

  // set the volume for the solver or assign a volume protocol
  if (V_protocol == false)
    {
      solver.set_V(V);
    }

  // resize the solver arrays given N_species and M_rxns

  // set the solver's state vector to the species counts

  // set the solver's FPT conditions

  // set the solver's reaction events

  
  // run the solver
  solver.run_FPT(t,dt,dt_target);


  // get the species counts from the solver
  
  // update the initiator distribution from the species counts
  
}
