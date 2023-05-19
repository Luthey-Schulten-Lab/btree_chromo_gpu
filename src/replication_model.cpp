#include <replication_model.hpp>

// constructor
replication_model::replication_model()
{
  r_m_p.k_rep = 0.0;
  r_m_p.k_c = 0.0;
  r_m_p.k_d = 0.0;
  r_m_p.k_hi = 0.0;
  r_m_p.k_lo = 0.0;
  r_m_p.k_on = 0.0;
  r_m_p.k_off = 0.0;
  r_m_p.k_bubble = 0.0;

  r_m_p.N_hi = 1;
  r_m_p.N_lo = 1;
  r_m_p.N_fil = 1;
}

// destructor
replication_model::~replication_model()
{
}

// read the replication model
void replication_model::load_model(std::string rep_model_filename)
{

  std::fstream rep_model_file;

  std::string param_delim, param, val;
  int delim;
  
  std::string line;

  param_delim = "=";

  rep_model_file.open(rep_model_filename, std::ios::in);

  if (!rep_model_file.is_open())
    {
      std::cout << "ERROR: file not opened in read_rep_model" << std::endl;
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

		  // std::cout << param << "=" << val << std::endl;

		  if (param == "k_rep")
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

		  else if (param == "k_bubble")
		    {
		      r_m_p.k_bubble = stod(val);
		    }

		}
	      	      
	    }
	      
     	} // end while loop

      rep_model_file.close();
  
    }

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
  number_rep_species();
  number_rep_rxns();
}


// set the number of species based on the replication model
void replication_model::number_rep_species()
{
  N_species = 0;

  N_non_leaf = 2; // free DnaA and DnaA genes
  
  N_per_leaf = 1; // empty origin
  N_per_leaf += r_m_p.N_hi; // high affinity sites
  N_per_leaf += r_m_p.N_lo; // low affinity sites
  N_per_leaf += r_m_p.N_fil; // filament sites

  N_species = N_leaves*N_per_leaf + N_non_leaf;
}


// set the number of reactions based on the replication model
void replication_model::number_rep_rxns()
{
  M_rxns = 0;

  M_non_leaf = 2; // free DnaA creation and destruction
  
  M_per_leaf = 0;
  M_per_leaf += r_m_p.N_hi; // high affinity site binding beginning with empty origin
  M_per_leaf += r_m_p.N_lo; // low affinity site binding beginning with last high affinity site
  M_per_leaf += 2*(r_m_p.N_fil - 1); // filament addition beginning with last low affinity site
  M_per_leaf += 1; // filament to bubble

  M_rxns = N_leaves*M_per_leaf + M_non_leaf;
}

// prepare a vector of the reactions
std::vector<reaction> replication_model::get_reactions()
{

  std::vector<reaction> rxns;
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



std::vector<species_count> replication_model::create_xFPT()
{
  std::vector<species_count> s_cs;
  species_count s_c;

  for (int i=0; i<N_leaves; i++)
    {
      s_c.id = (N_non_leaf - 1) + (i+1)*N_per_leaf;
      s_c.N = 1;
      s_cs.push_back(s_c);
    }
  return s_cs;
}
