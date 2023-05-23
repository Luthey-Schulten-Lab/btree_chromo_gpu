#include <replicator.hpp>

// constructor
replicator::replicator()
{
  V_protocol = false;
  n_DnaA = 0;
  N_leaves = 1;
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
  solver.prng_seed(s);
}


// reset the distribution of initiators
void replicator::reset_init_dist()
{
  init_dist.clear();
}


// update the distribution of initiators across leaves
void replicator::update_init_dist(std::vector<std::string> leaves)
{
  // update the initiator distribution if there are new leaves
  if (init_dist.size() != leaves.size())
    {
      
      init_loc temp_i_l;
      std::vector<init_loc> old_init_dist = init_dist;

      init_dist.clear();
      for (size_t i=0; i<leaves.size(); i++)
	{
	  temp_i_l.loc = leaves[i];
	  temp_i_l.N = 0;
	  init_dist.push_back(temp_i_l);

	  for (size_t j=0; j<old_init_dist.size(); j++)
	    {
	      if (init_dist[i].loc == old_init_dist[j].loc)
		{
		  init_dist[i].N = old_init_dist[j].N;
		}
	    }
	}

      // determine the number of Oris
      N_leaves = static_cast<int>(init_dist.size());
      // set the number of Oris for the reaction model
      rep_model.set_N_leaves(N_leaves);
      
    }
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

  V_points.clear();
  V_protocol = true;

  std::fstream volume_protocol_file;

  std::string param_delim;
  int delim;
  
  std::string line;

  volume_point v_p;

  param_delim = ",";

  volume_protocol_file.open(volume_protocol_filename, std::ios::in);

  if (!volume_protocol_file.is_open())
    {
      std::cout << "ERROR: file not opened in load_model" << std::endl;
    }
  else
    {
      while (1)
	{
	  volume_protocol_file >> line;
	  if (volume_protocol_file.eof()) break;
	  

	  if ((line.length() > 0) &&
	      (line.find("#") != 0))
	    {

	      delim = line.find(param_delim);

	      if (delim != -1)
		{

		  v_p.t = stod(line.substr(0,delim));
		  v_p.V = stod(line.substr(delim+1,line.length()));

		  V_points.push_back(v_p);

		}
	      	      
	    }
	      
     	} // end while loop

      volume_protocol_file.close();
  
    }

}


//load the replication model
void replicator::load_model(std::string rep_model_filename)
{

  std::cout << rep_model_filename << std::endl;
  rep_model.load_model(rep_model_filename);
}


// get the time
double replicator::get_time()
{
  return t;
}


// get the replication rate
double replicator::get_k_rep()
{
  return rep_model.get_k_rep();
}


// get the maximum number of replisomes
int replicator::get_max_replisomes()
{
  return max_rep;
}


// get the number of free DnaA
int replicator::get_DnaA()
{
  return n_DnaA;
}


// test for an initiation event
std::string replicator::initiation_test()
{
  // assume by default that there is no initiation
  std::string leaf_init = "no";

  // get the initiation requirements from the model
  std::vector<std::array<int,2>> init_reqs = rep_model.get_init_requirements();

  // loop over the locations
  for (size_t i=0; i<init_dist.size(); i++)
    {
      // loop over the possible bubble formation outcomes
      for (size_t j=0; j<init_reqs.size(); j++)
	{
	  if (init_dist[i].N == init_reqs[j][0])
	    {
	      // set the intiated leaf to the location in the distribution
	      leaf_init = init_dist[i].loc;
	      // return bound DnaA to pool of free DnaA
	      n_DnaA += init_reqs[j][1];
	      return leaf_init;
	    }
	}
    }

  return leaf_init;
}


// prepare the species counts
std::vector<species_count> replicator::prepare_species_counts()
{

  std::vector<species_count> s_cs;
  species_count s_c;
  int N_non_leaf = rep_model.get_N_non_leaf();
  int N_per_leaf = rep_model.get_N_per_leaf();

  // DnaA genes
  s_c.id = 0;
  s_c.N = n_DnaA_genes;
  s_cs.push_back(s_c);

  // free DnaA
  s_c.id = 1;
  s_c.N = n_DnaA;
  s_cs.push_back(s_c);

  // convert the initiator distributions
  for (size_t i=0; i<init_dist.size(); i++)
    {
      s_c.id = N_non_leaf + i*N_per_leaf + init_dist[i].N;
      s_c.N = 1;
      s_cs.push_back(s_c);
    }

  return s_cs;
}


// create the FPT vector
std::vector<species_count> replicator::prepare_xFPT()
{
  std::vector<species_count> s_cs;
  species_count s_c;
  int N_non_leaf = rep_model.get_N_non_leaf();
  int N_per_leaf = rep_model.get_N_per_leaf();
  std::vector<std::array<int,2>> init_reqs = rep_model.get_init_requirements();

  for (int i=0; i<N_leaves; i++)
    {
      for (size_t j=0; j<init_reqs.size(); j++)
	{
	  s_c.id = N_non_leaf + i*N_per_leaf;
	  s_c.id += init_reqs[j][0];
	  s_c.N = 1;
	  s_cs.push_back(s_c);
	  std::cout << "FPT[" << i << "],"
		    << "id = " << s_c.id << " at "
		    << "n = " << s_c.N << std::endl;
	}
    }

  return s_cs;
}


// convert the species counts to the replicator state
void replicator::species_counts_to_replicator_state(std::vector<species_count> s_cs)
{
  int N_non_leaf = rep_model.get_N_non_leaf();
  int N_per_leaf = rep_model.get_N_per_leaf();

  int id;
  int j_loc, j_count;

  for (size_t i=0; i<s_cs.size(); i++)
    {

      id = s_cs[i].id;
      
      if (id < N_non_leaf)
	{
	  if (id == 0)
	    {
	      n_DnaA_genes = s_cs[i].N;
	    }
	  else if (id == 1)
	    {
	      n_DnaA = s_cs[i].N;
	    }
	}
      else
	{
	  id -= N_non_leaf;

	  j_loc = id/N_per_leaf;
	  j_count = id%N_per_leaf;

	  init_dist[j_loc].N = j_count;

	}

    }
}


// print the replicator's state
void replicator::print_state()
{
  std::cout << "n_DnaA_genes = " << n_DnaA_genes << std::endl;
  std::cout << "n_DnaA = " << n_DnaA << std::endl;
  for (size_t i=0; i<init_dist.size(); i++)
    {
      std::cout << " " << init_dist[i].loc
		<< "," << init_dist[i].N << std::endl;
    }
}


// dump the replicator's state
void replicator::dump_state(std::string state_filename)
{
  std::fstream state_file;

  state_file.open(state_filename, std::ios::out);

  if (!state_file.is_open())
    {
      std::cout << "ERROR: file not opened in dump_state" << std::endl;
    }
  else
    {
      state_file << "n_DnaA=" << n_DnaA << std::endl;
      state_file << "N_leaves=" << N_leaves << std::endl;

      for (size_t i=0; i<init_dist.size(); i++)
	{
	  state_file << init_dist[i].loc << ","
		     << init_dist[i].N << std::endl;
	}
      state_file.close();
    }
  
}


// load the replicator's state
void replicator::load_state(std::string state_filename)
{
  std::cout << state_filename << std::endl;
}


//run the replicator for a time increment or until a first-passage
void replicator::run(double &dt, double dt_target)
{
  
  std::cout << "dt_target = " << dt_target << std::endl;

  // set the number of DnaA genes equal to the number of Oris
  set_DnaA_genes(N_leaves);

  std::cout << "preparing reaction system" << std::endl;
  
  // prepare a vector of the initial species counts with the rep model
  std::vector<species_count> x0 = prepare_species_counts();
  
  // prepare a vector of FPT events
  std::vector<species_count> xFPT = prepare_xFPT();

  // prepare the reaction system
  std::vector<reaction> rxns = rep_model.get_reactions();
  
  // initialize the solver

  // set the volume for the solver or assign a volume protocol
  if (V_protocol == false)
    {
      solver.set_V_constant(V);
    }
  else
    {
      solver.set_V_points(V_points);
    }

  std::cout << "preparing solver" << std::endl;
  
  // resize the solver arrays given N_species and M_rxns
  solver.initialize_reaction_system(rep_model.get_N_species(),
				    rep_model.get_M_rxns());
  
  // set the solver's state vector to the species counts
  solver.set_x(x0);
  
  // set the solver's FPT conditions
  solver.set_xFPT(xFPT);
  
  // set the solver's reaction events
  solver.set_reaction_system(rxns);

  std::cout << "running solver" << std::endl;
  
  // run the solver
  solver.run_FPT(t,dt,dt_target);

  std::cout << "analyzing solution" << std::endl;

  // get the species counts from the solver
  std::vector<species_count> x = solver.state_to_sc();

  // update the initiator distribution from the species counts
  species_counts_to_replicator_state(x);

  print_state();

  // increment the time
  t += dt;
  
}
