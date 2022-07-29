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
  number_rep_species();
  number_rep_rxns();
  vector<reaction> rxns;

  cout << "N_species = " << N_species << endl;
  cout << "N_rxns = " << M_rxns << endl;

  solver.prepare_reaction_system(N_species,M_rxns);

  get_reactions(rxns);
  
  solver.set_S(rxns);
  
}

// get the number of species based on the replication model
void replication_model::number_rep_species()
{
  N_species = 0;
  
  N_per_leaf = 1; // empty origin

  N_per_leaf += r_m_p.N_hi; // high affinity sites
  N_per_leaf += r_m_p.N_lo; // low affinity sites
  N_per_leaf += r_m_p.N_fil; // filament sites

  N_species = N_leaves*N_per_leaf;

  N_species += 1; // free DnaA
}


// get the number of reactions based on the replication model
void replication_model::number_rep_rxns()
{
  M_rxns = 0;

  M_rxns += r_m_p.N_hi; // high affinity site binding beginning with empty origin
  M_rxns += r_m_p.N_lo; // low affinity site binding beginning with last high affinity site
  M_rxns += 2*(r_m_p.N_fil - 1); // filament addition beginning with last low affinity site
  M_rxns += 1; // filament to bubble

  M_rxns = N_leaves*M_rxns;

  M_rxns += 2; // free DnaA creation and destruction
}

// prepare a vector of the reactions
void replication_model::get_reactions(vector<reaction> &rxns)
{

  reaction r;
  int df = 1;
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
	  solver.add_reaction_input(r,df+i*N_per_leaf+c+j,1);
	  solver.add_reaction_output(r,df+i*N_per_leaf+c+j+1,1);
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
	  solver.add_reaction_input(r,df+i*N_per_leaf+c+j,1);
	  solver.add_reaction_output(r,df+i*N_per_leaf+c+j+1,1);
	  rxns.push_back(r);
	  solver.reset_reaction(r);
	  db += 1;
	}

      c += db;

      // add reactions for filamentation
      for (int j=0; j<r_m_p.N_fil; j++)
	{
	  // input is empty hi site and free DnaA
	  solver.add_reaction_input(r,0,1);
	  solver.add_reaction_input(r,df+i*N_per_leaf+c+j,1);
	  solver.add_reaction_output(r,df+i*N_per_leaf+c+j+1,1);
	  rxns.push_back(r);
	  solver.reset_reaction(r);
	  if (j < r_m_p.N_fil-1)
	    {
	      solver.add_reaction_output(r,0,1);
	      solver.add_reaction_output(r,df+i*N_per_leaf+c+j,1);
	      solver.add_reaction_input(r,df+i*N_per_leaf+c+j+1,1);
	      rxns.push_back(r);
	      solver.reset_reaction(r);
	    }
	}
      
    }

}

vector<species_count> replication_model::id_to_sc(vector<init_loc> &init_dist)
{
  vector<species_count> s_cs;
  species_count s_c;

  s_c.id = 0;
  s_c.n = init_dist[0].N;

  s_cs.push_back(s_c);

  for (int i=0; i<N_leaves; i++)
    {
      s_c.id = 1 + i*N_per_leaf + init_dist[i+1].N;
      s_c.n = 1;
      s_cs.push_back(s_c);
    }

  return s_cs;
}

void replication_model::update_id_from_sc(vector<init_loc> &init_dist, vector<species_count> s_cs)
{
  for (species_count s_c : s_cs)
    {
      if (s_c.id > 0)
	{
	  init_dist[(s_c.id-1)/N_per_leaf+1].N = (s_c.id-1)%N_per_leaf;
	}
      else
	{
	  init_dist[0].N = s_c.n;
	}
    }
}

vector<species_count> replication_model::create_xFPT()
{
  vector<species_count> s_cs;
  species_count s_c;

  for (int i=0; i<N_leaves; i++)
    {
      s_c.id = (i+1)*N_per_leaf;
      s_c.n = 1;
      s_cs.push_back(s_c);
    }
  return s_cs;
}

void replication_model::propensities(int *xf, double *Wf)
{
  int df = 1;
  int c, db;
  int k;
  
  Wf[0] = r_m_p.k_c;
  Wf[1] = r_m_p.k_d*xf[0];

  k = 2;

  for (int i=0; i<N_leaves; i++)
    {
      c = 0;
            
      // add reactions for high-affinity sites
      db = 0;
      for (int j=0; j<r_m_p.N_hi; j++)
	{
	  // input is empty hi site and free DnaA
	  Wf[k] = r_m_p.k_hi*xf[0]*xf[df+i*N_per_leaf+c+j];
	  k += 1;
	  db += 1;
	}

      c += db;

      // add reactions for low-affinity sites
      db = 0;
      for (int j=0; j<r_m_p.N_lo; j++)
	{
	  // input is empty lo site and free DnaA
	  Wf[k] = r_m_p.k_lo*xf[0]*xf[df+i*N_per_leaf+c+j];
	  k += 1;
	  db += 1;
	}

      c += db;

      // add reactions for filamentation
      for (int j=0; j<r_m_p.N_fil; j++)
	{
	  // input is empty hi site and free DnaA
	  Wf[k] = r_m_p.k_on*xf[0]*xf[df+i*N_per_leaf+c+j];
	  k += 1;
	  if (j < r_m_p.N_fil-1)
	    {
	      Wf[k] = r_m_p.k_off*xf[df+i*N_per_leaf+c+j+1];
	      k += 1;
	    }
	}
      
    }

  
}

void replication_model::run_replicate_FPT(rep_model_params &rep_model, vector<init_loc> &init_dist, double &t, double &t_max)
{

  // prepare the reaction system based on the reaction model
  prepare_system(rep_model,static_cast<int>(init_dist.size())-1);

  // convert the initiator distribution to species counts
  vector<species_count> s_cs = id_to_sc(init_dist);
  
  // set the state vector
  solver.set_x(s_cs);

  // set the FPT state vector
  solver.set_xFPT(create_xFPT());

  // set the propensity function
  solver.set_propensity_fxn(propensities);

  // print the system as a sanity check
  solver.print_reaction_system();

  // run the system

  // redistribute initiators based on FPT result
  update_id_from_sc(init_dist,solver.state_to_sc());

  // for (init_loc temp_i_l : init_dist)
  //   {
  //     cout << temp_i_l.loc << " = " << temp_i_l.N << endl;
  //   }
  
  t += t_max/2.0;
  init_dist.back().N = -1;
}
