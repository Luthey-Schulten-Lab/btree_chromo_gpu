#include <btree/gillespie_solver.hpp>

// constructor
gillespie_solver::gillespie_solver()
{
  rand_eng.seed(0);
  u_rand = uniform_real_distribution<double>(0.0,1.0);
}

// destructor
gillespie_solver::~gillespie_solver()
{
}

// prng seed
void gillespie_solver::prng_seed(int s)
{
  rand_eng.seed(s);
}


// read the replication model
rep_model_params gillespie_solver::read_rep_model(string rep_model_filename)
{
  rep_model_params r_i_p;

  fstream rep_model_file;

  string param_delim, param, val;
  int delim;
  
  string line;

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
void gillespie_solver::prepare_system(rep_model_params r_m_p, int N_forks)
{
  int N_species = number_species(r_m_p,N_forks);
  int N_rxns = number_rxns(r_m_p,N_forks);

  cout << "N_species = " << N_species << endl;
  cout << "N_rxns = " << N_rxns << endl;
  
}

// get the number of species based on the replication model
int gillespie_solver::number_species(rep_model_params r_m_p, int N_forks)
{
  int N_species = 0;

  N_species += 1; // free DnaA

  N_species += N_forks*(2+r_m_p.N_hi+r_m_p.N_lo+r_m_p.N_fil);

  return N_species;
}

// get the number of reactions based on the replication model
int gillespie_solver::number_rxns(rep_model_params r_m_p, int N_forks)
{
  int N_rxns = 0;

  N_rxns += 1; // free DnaA

  N_rxns += N_forks*(2+r_m_p.N_hi+r_m_p.N_lo+r_m_p.N_fil);
  
  return N_rxns;
}
