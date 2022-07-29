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
void replication_model::prepare_system(rep_model_params rep_model, int n)
{

  r_m_p = rep_model;
  N_leaves = n;
  N_species = number_rep_species();
  M_rxns = number_rep_rxns();
  vector<reaction> rxns;

  cout << "N_species = " << N_species << endl;
  cout << "N_rxns = " << M_rxns << endl;

  solver.prepare_reaction_system(N_species,M_rxns);

  get_reactions(rxns);
  
  solver.set_S(rxns);
  solver.print_reaction_system();
  
}

// get the number of species based on the replication model
int replication_model::number_rep_species()
{
  int N_species = 0;

  N_species += 1; // empty origin and origin with bubble

  N_species += r_m_p.N_hi; // high affinity sites
  N_species += r_m_p.N_lo; // low affinity sites
  N_species += r_m_p.N_fil; // filament sites

  N_species = N_leaves*N_species;

  N_species += 1; // free DnaA

  return N_species;
  
}


// get the number of reactions based on the replication model
int replication_model::number_rep_rxns()
{
  int N_rxns = 0;

  N_rxns += r_m_p.N_hi; // high affinity site binding beginning with empty origin
  N_rxns += r_m_p.N_lo; // low affinity site binding beginning with last high affinity site
  N_rxns += 2*(r_m_p.N_fil - 1); // filament addition beginning with last low affinity site
  N_rxns += 1; // filament to bubble

  N_rxns = N_leaves*N_rxns;

  N_rxns += 2; // free DnaA creation and destruction
  
  return N_rxns;
}

void replication_model::get_reactions(vector<reaction> &rxns)
{

  reaction r;
  int df = 1;
  int dl = r_m_p.N_hi + r_m_p.N_lo + r_m_p.N_fil + 1;
  int c, db;

  rxns.clear();
  
  solver.reset_reaction(r);

  // DnaA creation
  solver.add_reaction_output(r,0,1);
  rxns.push_back(r);
  
  solver.reset_reaction(r);

  // DnaA destruction
  solver.add_reaction_input(r,0,1);
  rxns.push_back(r);
  solver.reset_reaction(r);

  for (int i=0; i<N_leaves; i++)
    {
      c = 0;
            
      // add reactions for high-affinity sites
      db = 0;
      for (int j=0; j<r_m_p.N_hi; j++)
	{
	  // input is empty hi site and free DnaA
	  solver.add_reaction_input(r,0,1);
	  solver.add_reaction_input(r,df+i*dl+c+j,1);
	  solver.add_reaction_output(r,df+i*dl+c+j+1,1);
	  rxns.push_back(r);
	  solver.reset_reaction(r);
	  db += 1;
	}

      c += db;

      // add reactions for low-affinity sites
      db = 0;
      for (int j=0; j<r_m_p.N_lo; j++)
	{
	  // input is empty hi site and free DnaA
	  solver.add_reaction_input(r,0,1);
	  solver.add_reaction_input(r,df+i*dl+c+j,1);
	  solver.add_reaction_output(r,df+i*dl+c+j+1,1);
	  rxns.push_back(r);
	  solver.reset_reaction(r);
	  db += 1;
	}

      c += db;

      // add reactions for filamentation
      db = 0;
      for (int j=0; j<r_m_p.N_fil; j++)
	{
	  // input is empty hi site and free DnaA
	  solver.add_reaction_input(r,0,1);
	  solver.add_reaction_input(r,df+i*dl+c+j,1);
	  solver.add_reaction_output(r,df+i*dl+c+j+1,1);
	  rxns.push_back(r);
	  solver.reset_reaction(r);
	  if (j < r_m_p.N_fil-1)
	    {
	      solver.add_reaction_output(r,0,1);
	      solver.add_reaction_output(r,df+i*dl+c+j,1);
	      solver.add_reaction_input(r,df+i*dl+c+j+1,1);
	      rxns.push_back(r);
	      solver.reset_reaction(r);
	      db += 2;
	    }
	  else
	    {
	      db += 1;
	    }
	}
      
    }

  cout << rxns.size() << endl;
}

void replication_model::run_replicate_FPT(rep_model_params &rep_model, vector<init_loc> &init_dist, double &t, double &t_max)
{
  prepare_system(rep_model,static_cast<int>(init_dist.size())-1);
  t += t_max/2.0;
  init_dist.back().N = -1;
}
