#include <rep_kinetics/replication_model.hpp>

// constructor
replication_model::replication_model()
{
}

// destructor
replication_model::~replication_model()
{
}

// read the replication model
void replication_model::read_rep_model(string rep_model_filename)
{

  fstream rep_model_file;

  string param_delim, param, val;
  int delim;
  
  string line;

  param_delim = "=";

  rep_model_file.open(rep_model_filename, ios::in);

  if (!rep_model_file)
    {
      cout << "ERROR: file not opened in read_rep_model" << endl;
    }
  else
    {
      while (1)
	{
	  rep_model_file >> line;
	  if (rep_model_file.eof()) break;
	  

	  if ((line.length() > 0) &&
	      (line.find("#") != 0))
	    {

	      delim = line.find(param_delim);

	      if (delim != -1)
		{

		  param = line.substr(0,delim);
		  val = line.substr(delim+1,line.length());

		  // cout << param << "=" << val << endl;

		  if (param == "V")
		    {
		      r_m_p.V = stod(val);
		    }

		  else if (param == "max_replisomes")
		    {
		      r_m_p.max_replisomes = stoi(val);
		    }

		  else if (param == "N_init_DnaA")
		    {
		      r_m_p.N_init_DnaA = stoi(val);
		    }

		  else if (param == "k_rep")
		    {
		      r_m_p.k_rep = stod(val);
		    }
		  
		  else if (param == "k_c")
		    {
		      r_m_p.k_c = stod(val);
		    }

		  else if (param == "k_d")
		    {
		      r_m_p.k_d = stod(val);
		    }

		  else if (param == "N_hi")
		    {
		      r_m_p.N_hi = stoi(val);
		    }

		  else if (param == "k_hi")
		    {
		      r_m_p.k_hi = stod(val);
		    }

		  else if (param == "N_lo")
		    {
		      r_m_p.N_lo = stoi(val);
		    }

		  else if (param == "k_lo")
		    {
		      r_m_p.k_lo = stod(val);
		    }

		  else if (param == "N_fil")
		    {
		      r_m_p.N_fil = stoi(val);
		    }

		  else if (param == "k_on")
		    {
		      r_m_p.k_on = stod(val);
		    }

		  else if (param == "k_off")
		    {
		      r_m_p.k_off = stod(val);
		    }

		  else if (param == "tau_SA")
		    {
		      r_m_p.tau_SA = stod(val);
		    }

		  else if (param == "N_init_SA")
		    {
		      r_m_p.N_init_SA = stoi(val);
		    }

		}
	      	      
	    }
	      
     	} // end while loop
  
    }

  r_m_p.k_SA = r_m_p.N_init_SA/r_m_p.tau_SA;
  
  rep_model_file.close();

}

// get N_init_DnaA
int replication_model::get_N_init_DnaA()
{
  return r_m_p.N_init_DnaA;
}

// get N_leaves
int replication_model::get_N_leaves()
{
  return N_leaves;
}

// get max_replisomes
int replication_model::get_max_replisomes()
{
  return r_m_p.max_replisomes;
}

// get k_rep
double replication_model::get_k_rep()
{
  return r_m_p.k_rep;
}

// get N_species
int replication_model::get_N_species()
{
  return N_species;
}

// get M_rxns
int replication_model::get_M_rxns()
{
  return M_rxns;
}


// get N_per_leaf
int replication_model::get_N_per_leaf()
{
  return N_per_leaf;
}

// get N_per_leaf
int replication_model::get_N_non_leaf()
{
  return N_non_leaf;
}


// set the number of leaves
void replication_model::set_N_leaves(int N_leaves)
{
  this->N_leaves = N_leaves;
}


// set the number of species based on the replication model
void replication_model::number_rep_species()
{
  N_species = 0;
  
  N_per_leaf = 1; // empty origin

  N_per_leaf += r_m_p.N_hi; // high affinity sites
  N_per_leaf += r_m_p.N_lo; // low affinity sites
  N_per_leaf += r_m_p.N_fil; // filament sites

  N_species = N_leaves*N_per_leaf;

  N_species += N_non_leaf;
  
}


// set the number of reactions based on the replication model
void replication_model::number_rep_rxns()
{
  M_rxns = 0;

  M_rxns += r_m_p.N_hi; // high affinity site binding beginning with empty origin
  M_rxns += r_m_p.N_lo; // low affinity site binding beginning with last high affinity site
  M_rxns += 2*(r_m_p.N_fil - 1); // filament addition beginning with last low affinity site
  M_rxns += 1; // filament to bubble

  M_rxns = N_leaves*M_rxns;

  M_rxns += 2; // free DnaA creation and destruction
  M_rxns += 1; // SA particle creation
}

// prepare a vector of the reactions
vector<reaction> replication_model::get_reactions()
{

  vector<reaction> rxns;
  reaction r;
  int df = N_non_leaf;
  int free_idx = N_non_leaf - 1;
  int c, db;

  // SA particle creation
  rxn_manip.add_reaction_output(r,0,1); // 1 free SA particle out
  rxns.push_back(r);
  rxn_manip.reset_reaction(r);

  // DnaA creation
  rxn_manip.add_reaction_input(r,1,1); // 1 gene in
  rxn_manip.add_reaction_output(r,1,1); // 1 gene out
  rxn_manip.add_reaction_output(r,free_idx,1); // 1 free DnaA out
  rxns.push_back(r);
  rxn_manip.reset_reaction(r);

  // DnaA destruction
  rxn_manip.add_reaction_input(r,free_idx,1); // 1 free DnaA in
  rxns.push_back(r);
  rxn_manip.reset_reaction(r);

  for (int i=0; i<N_leaves; i++)
    {
      c = 0;
            
      // add reactions for high-affinity sites
      db = 0;
      for (int j=0; j<r_m_p.N_hi; j++)
	{
	  // input is empty hi site and free DnaA
	  rxn_manip.add_reaction_input(r,free_idx,1);
	  rxn_manip.add_reaction_input(r,df+i*N_per_leaf+c+j,1);
	  rxn_manip.add_reaction_output(r,df+i*N_per_leaf+c+j+1,1);
	  rxns.push_back(r);
	  rxn_manip.reset_reaction(r);
	  db += 1;
	}

      c += db;

      // add reactions for low-affinity sites
      db = 0;
      for (int j=0; j<r_m_p.N_lo; j++)
	{
	  // input is empty lo site and free DnaA
	  rxn_manip.add_reaction_input(r,free_idx,1);
	  rxn_manip.add_reaction_input(r,df+i*N_per_leaf+c+j,1);
	  rxn_manip.add_reaction_output(r,df+i*N_per_leaf+c+j+1,1);
	  rxns.push_back(r);
	  rxn_manip.reset_reaction(r);
	  db += 1;
	}

      c += db;

      // add reactions for filamentation
      for (int j=0; j<r_m_p.N_fil; j++)
	{
	  // input is empty hi site and free DnaA
	  rxn_manip.add_reaction_input(r,free_idx,1);
	  rxn_manip.add_reaction_input(r,df+i*N_per_leaf+c+j,1);
	  rxn_manip.add_reaction_output(r,df+i*N_per_leaf+c+j+1,1);
	  rxns.push_back(r);
	  rxn_manip.reset_reaction(r);
	  // defilamentation reaction
	  if (j < r_m_p.N_fil-1)
	    {
	      rxn_manip.add_reaction_output(r,free_idx,1);
	      rxn_manip.add_reaction_output(r,df+i*N_per_leaf+c+j,1);
	      rxn_manip.add_reaction_input(r,df+i*N_per_leaf+c+j+1,1);
	      rxns.push_back(r);
	      rxn_manip.reset_reaction(r);
	    }
	}
      
    }
  return rxns;
}


// initialize array of noninitiator species
void replication_model::reset_noninit_s(int * &noninit_s)
{

  if (noninit_s != nullptr)
    {
      delete[] noninit_s;
      noninit_s = nullptr;
    }

  N_non_leaf = 1; // free DnaA
  N_non_leaf += 1; // SA particles
  N_non_leaf += 1; // DnaA genes

  noninit_s = new int[N_non_leaf-1];

  for (int i=0; i<N_non_leaf-1; i++)
    {
      noninit_s[i] = 0;
    }

  // SA particles
  noninit_s[0] = r_m_p.N_init_SA;

  // DnaA genes, there has to be at least 1 gene copy by default
  noninit_s[1] = 1;
  
}


// update the noninitiator species
void replication_model::update_noninit_s(int *noninit_s)
{
  noninit_s[1] = N_leaves;
}


vector<species_count> replication_model::id_to_sc(vector<init_loc> &init_dist)
{
  vector<species_count> s_cs;
  species_count s_c;

  s_c.id = N_non_leaf - 1;
  s_c.N = init_dist[0].N;

  s_cs.push_back(s_c);

  for (int i=0; i<N_leaves; i++)
    {
      s_c.id = N_non_leaf + i*N_per_leaf + init_dist[i+1].N;
      s_c.N = 1;
      s_cs.push_back(s_c);
    }

  return s_cs;
}

void replication_model::update_init_from_solver_s_cs(vector<init_loc> &init_dist,
							vector<species_count> &solver_s_cs)
{
  for (species_count s_c : solver_s_cs)
    {
      if (s_c.id > N_non_leaf - 1)
	{
	  // cout << (s_c.id-N_non_leaf)/N_per_leaf+1 << " " << (s_c.id-N_non_leaf)%N_per_leaf << endl;
	  init_dist[(s_c.id-N_non_leaf)/N_per_leaf+1].N = (s_c.id-N_non_leaf)%N_per_leaf;
	}
      else if(s_c.id == N_non_leaf - 1)
	{
	  init_dist[0].N = s_c.N;
	}
    }
}

void replication_model::update_noninit_s_from_solver_s_cs(int *noninit_s,
							  vector<species_count> &solver_s_cs)
{
  for (int i=0; i<N_non_leaf-1; i++)
    {
      noninit_s[i] = 0;
    }
  for (species_count &s_c : solver_s_cs)
    {
      if (s_c.id < N_non_leaf - 1)
	{
	  noninit_s[s_c.id] = s_c.N;
	}
    }
}

vector<species_count> replication_model::create_xFPT()
{
  vector<species_count> s_cs;
  species_count s_c;

  for (int i=0; i<N_leaves; i++)
    {
      s_c.id = (N_non_leaf - 1) + (i+1)*N_per_leaf;
      s_c.N = 1;
      s_cs.push_back(s_c);
    }
  return s_cs;
}

void replication_model::propensities(int *xf, double *Wf)
{
  int df = N_non_leaf;
  int free_idx = N_non_leaf - 1;
  int c, db;
  int k;
  double inv_V = 1.0/r_m_p.N_init_SA;

  inv_V = xf[0]*inv_V;
  inv_V = pow(inv_V,3.0/2.0);
  inv_V = 1.0/(min(inv_V,2.0)*r_m_p.V); // calculate the inverse volume based on the SA particle change

  // create SA particles at rate given by SA doubling time (tau_SA)
  if (xf[0] < 2*r_m_p.N_init_SA)
    {
      Wf[0] = r_m_p.k_SA;
    }
  else
    {
      Wf[0] = 0.0;
    }
    
  Wf[1] = r_m_p.k_c*xf[1]; // create DnaA at rate proportional to number DnaA genes, i.e. leaves
  Wf[2] = r_m_p.k_d*xf[2]; // degrade DnaA at rate based on presumed doubling time

  k = 3;

  for (int i=0; i<N_leaves; i++)
    {
      c = 0;
            
      // add reactions for high-affinity sites
      db = 0;
      for (int j=0; j<r_m_p.N_hi; j++)
	{
	  // input is empty hi site and free DnaA
	  Wf[k] = r_m_p.k_hi*xf[free_idx]*xf[df+i*N_per_leaf+c+j]*inv_V;
	  k += 1;
	  db += 1;
	}

      c += db;

      // add reactions for low-affinity sites
      db = 0;
      for (int j=0; j<r_m_p.N_lo; j++)
	{
	  // input is empty lo site and free DnaA
	  Wf[k] = r_m_p.k_lo*xf[free_idx]*xf[df+i*N_per_leaf+c+j]*inv_V;
	  k += 1;
	  db += 1;
	}

      c += db;

      // add reactions for filamentation
      for (int j=0; j<r_m_p.N_fil; j++)
	{
	  // input is empty hi site and free DnaA
	  Wf[k] = r_m_p.k_on*xf[free_idx]*xf[df+i*N_per_leaf+c+j]*inv_V;
	  k += 1;
	  if (j < r_m_p.N_fil-1)
	    {
	      Wf[k] = r_m_p.k_off*xf[df+i*N_per_leaf+c+j+1];
	      k += 1;
	    }
	}
      
    }

  
}

// ptr replication_model::get_propensity_fxn()
// {
//   return &replication_model::propensities;
// }

