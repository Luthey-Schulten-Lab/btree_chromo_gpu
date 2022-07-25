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
rep_model_params replication_model::read_rep_model(string rep_model_filename)
{
  rep_model_params r_i_p;

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
		      r_i_p.V = stod(val);
		    }

		  else if (param == "max_replisomes")
		    {
		      r_i_p.max_replisomes = stoi(val);
		    }

		  else if (param == "N_init_DnaA")
		    {
		      r_i_p.N_init_DnaA = stoi(val);
		    }

		  else if (param == "k_rep")
		    {
		      r_i_p.k_rep = stod(val);
		    }
		  
		  else if (param == "k_c")
		    {
		      r_i_p.k_c = stod(val);
		    }

		  else if (param == "k_d")
		    {
		      r_i_p.k_d = stod(val);
		    }

		  else if (param == "N_hi")
		    {
		      r_i_p.N_hi = stoi(val);
		    }

		  else if (param == "k_hi")
		    {
		      r_i_p.k_hi = stod(val);
		    }

		  else if (param == "N_lo")
		    {
		      r_i_p.N_lo = stoi(val);
		    }

		  else if (param == "k_lo")
		    {
		      r_i_p.k_lo = stod(val);
		    }

		  else if (param == "N_fil")
		    {
		      r_i_p.N_fil = stoi(val);
		    }

		  else if (param == "k_on")
		    {
		      r_i_p.k_on = stod(val);
		    }

		  else if (param == "k_off")
		    {
		      r_i_p.k_off = stod(val);
		    }

		}
	      	      
	    }
	      
     	} // end while loop
  
    }
  
  rep_model_file.close();


  return r_i_p;
}

// prepare the system state
void replication_model::prepare_system(rep_model_params r_m_p, int N_forks)
{
  int N_species = number_rep_species(r_m_p,N_forks);
  int N_rxns = number_rep_rxns(r_m_p,N_forks);

  cout << "N_species = " << N_species << endl;
  cout << "N_rxns = " << N_rxns << endl;

  solver.prepare_reaction_system(N_species,N_rxns);
  solver.print_reaction_system();
  
}

// get the number of species based on the replication model
int replication_model::number_rep_species(rep_model_params r_m_p, int N_forks)
{
  int N_species = 0;

  N_species += 2; // empty origin and origin with bubble

  N_species += r_m_p.N_hi; // high affinity sites
  N_species += r_m_p.N_lo; // low affinity sites
  N_species += r_m_p.N_fil; // filament sites

  N_species = N_forks*N_species;

  N_species += 1; // free DnaA

  return N_species;
  
}


// get the number of reactions based on the replication model
int replication_model::number_rep_rxns(rep_model_params r_m_p, int N_forks)
{
  int N_rxns = 0;

  N_rxns += r_m_p.N_hi; // high affinity site binding beginning with empty origin
  N_rxns += r_m_p.N_lo; // low affinity site binding beginning with last high affinity site
  N_rxns += r_m_p.N_fil; // filament addition beginning with last low affinity site
  N_rxns += 1; // filament to bubble

  N_rxns = N_forks*N_rxns;

  N_rxns += 2; // free DnaA creation and destruction
  
  return N_rxns;
}

void replication_model::run_replicate_FPT(rep_model_params r_m_p, int N_forks, double t_max)
{
  prepare_system(r_m_p,N_forks);
  cout << t_max << endl;
}
