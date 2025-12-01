#include <replication_model.hpp>

// constructor
replication_model::replication_model()
{
  r_m_p.k_rep = 0.0;
  r_m_p.k_c = 0.0;
  r_m_p.k_d = 0.0;
  r_m_p.k_hi = 0.0;
  r_m_p.k_lo = 0.0;
  r_m_p.k_melt_bind = 0.0;
  r_m_p.k_unbind_unmelt = 0.0;
  r_m_p.k_on = 0.0;
  r_m_p.k_off = 0.0;
  r_m_p.k_bubble = 0.0;

  r_m_p.N_hi = 1;
  r_m_p.N_lo = 1;
  r_m_p.N_fil = 1;
  r_m_p.bubble_min_fil = r_m_p.N_hi + r_m_p.N_lo + r_m_p.N_fil;

  G = 1.0;
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
      std::cout << "ERROR: file not opened in load_model" << std::endl;
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

		  else if (param == "k_melt_bind")
		    {
		      r_m_p.k_melt_bind = stod(val);
		    }

		  else if (param == "k_unbind_unmelt")
		    {
		      r_m_p.k_unbind_unmelt = stod(val);
		    }

		  else if (param == "k_on")
		    {
		      r_m_p.k_on = stod(val);
		    }

		  else if (param == "k_off")
		    {
		      r_m_p.k_off = stod(val);
		    }

		  else if (param == "bubble_min_fil")
		    {
		      r_m_p.bubble_min_fil = stoi(val);
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


// get N_bubble
int replication_model::get_N_bubble()
{
  return N_bubble;
}


// get N_binding
int replication_model::get_N_binding()
{
  return N_binding;
}


// get N_non_leaf
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
  prepare_reactions();
  prepare_init_requirements();
}


// set the scaled size
void replication_model::set_G(double G)
{
  this->G = G;
}


// set the number of species based on the replication model
void replication_model::number_rep_species()
{
  N_species = 0;

  N_non_leaf = 3; // free DnaA, DnaA genes, and free replisomes
  
  N_per_leaf = 1; // empty origin

  N_binding = 0; // total number of binding sites
  N_binding += r_m_p.N_hi; // high affinity sites
  N_binding += r_m_p.N_lo; // low affinity sites
  N_binding += r_m_p.N_fil; // filament sites
  
  N_per_leaf += N_binding; // binding sites
  N_bubble = (N_binding - r_m_p.bubble_min_fil + 1); // possible bubbles
 
  N_per_leaf += N_bubble; // bubble species

  N_species = N_leaves*N_per_leaf + N_non_leaf; // total species
}


// set the number of reactions based on the replication model
void replication_model::number_rep_rxns()
{
  M_rxns = 0;

  M_non_leaf = 2; // free DnaA creation and destruction
  
  M_per_leaf = 0;
  M_per_leaf += r_m_p.N_hi; // high affinity site binding beginning with empty origin
  M_per_leaf += r_m_p.N_lo; // low affinity site binding beginning with last high affinity site
  M_per_leaf += 2*r_m_p.N_fil; // filament addition beginning with last low affinity site
  M_per_leaf += N_bubble; // filament to bubble

  M_rxns = N_leaves*M_per_leaf + M_non_leaf;
}

// prepare a vector of the reactions
void replication_model::prepare_reactions()
{
  rxns.clear();
  reaction r;
  int df = N_non_leaf;
  int free_idx = N_non_leaf - 1;
  int c, db;

  // DnaA creation
  rxn_manip.add_reaction_input(r,0); // 1 gene in
  rxn_manip.add_reaction_output(r,0); // 1 gene out
  rxn_manip.add_reaction_output(r,free_idx); // 1 free DnaA out
  rxn_manip.add_reaction_rate(r,r_m_p.k_c);
  // rxn_manip.add_reaction_rate(r,r_m_p.k_c/G);
  rxns.push_back(r);
  rxn_manip.reset_reaction(r);

  // DnaA destruction
  rxn_manip.add_reaction_input(r,free_idx); // 1 free DnaA in
  rxn_manip.add_reaction_rate(r,r_m_p.k_d);
  // rxn_manip.add_reaction_rate(r,r_m_p.k_d*G);
  rxns.push_back(r);
  rxn_manip.reset_reaction(r);

  // create the set of initiation reactions for each Ori
  for (int i=0; i<N_leaves; i++)
    {
      c = 0;
            
      // add reactions for high-affinity sites
      db = 0;
      for (int j=0; j<r_m_p.N_hi; j++)
	{
	  // input is empty hi site and free DnaA
	  rxn_manip.add_reaction_input(r,free_idx);
	  rxn_manip.add_reaction_input(r,df+i*N_per_leaf+c+j);
	  rxn_manip.add_reaction_output(r,df+i*N_per_leaf+c+j+1);
	  rxn_manip.add_reaction_rate(r,r_m_p.k_hi);
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
	  rxn_manip.add_reaction_input(r,free_idx);
	  rxn_manip.add_reaction_input(r,df+i*N_per_leaf+c+j);
	  rxn_manip.add_reaction_output(r,df+i*N_per_leaf+c+j+1);
	  rxn_manip.add_reaction_rate(r,r_m_p.k_lo);
	  rxns.push_back(r);
	  rxn_manip.reset_reaction(r);
	  db += 1;
	}

      c += db;

      // add reactions for filamentation
      db = 0;
      
      // first filament is formed from melting and binding at low site
      // input is L_(N_lo) site and free DnaA
      rxn_manip.add_reaction_input(r,free_idx);
      rxn_manip.add_reaction_input(r,df+i*N_per_leaf+c);
      // output is F_(1)
      rxn_manip.add_reaction_output(r,df+i*N_per_leaf+c+1);
      rxn_manip.add_reaction_rate(r,r_m_p.k_melt_bind);
      rxns.push_back(r);
      rxn_manip.reset_reaction(r);
      db += 1;

      // filament is lost from unbinding and unmelting at low site
      // input is F_(1)
      rxn_manip.add_reaction_input(r,df+i*N_per_leaf+c+1);
      // output is L_(N_lo) and free DnaA
      rxn_manip.add_reaction_output(r,free_idx);
      rxn_manip.add_reaction_output(r,df+i*N_per_leaf+c);
      rxn_manip.add_reaction_rate(r,r_m_p.k_unbind_unmelt);
      rxns.push_back(r);
      rxn_manip.reset_reaction(r);

      
      // filaments grow and degrade
      for (int j=1; j<r_m_p.N_fil; j++)
	{

	  // filamentation reaction
	  // input is F_(i-1) site and free DnaA
	  rxn_manip.add_reaction_input(r,free_idx);
	  rxn_manip.add_reaction_input(r,df+i*N_per_leaf+c+j);
	  // output is F(i)
	  rxn_manip.add_reaction_output(r,df+i*N_per_leaf+c+j+1);
	  rxn_manip.add_reaction_rate(r,r_m_p.k_on);
	  rxns.push_back(r);
	  rxn_manip.reset_reaction(r);
	  db += 1;
	  
	  // defilamentation reaction
	  // input is F_(i+1)
	  rxn_manip.add_reaction_input(r,df+i*N_per_leaf+c+j+1);
	  // output is F_(i) and free DnaA
	  rxn_manip.add_reaction_output(r,free_idx);
	  rxn_manip.add_reaction_output(r,df+i*N_per_leaf+c+j);
	  rxn_manip.add_reaction_rate(r,r_m_p.k_off);
	  rxns.push_back(r);
	  rxn_manip.reset_reaction(r);

	}

      c += db;

      // add reactions for conversion of filaments to bubbles
      db = 0;
      for (int j=0; j<N_bubble; j++)
	{
	  // filament in
	  rxn_manip.add_reaction_input(r,df+i*N_per_leaf+c-N_bubble+j+1);
	  // replisome in
	  rxn_manip.add_reaction_input(r,1);
	  // bubble out
	  rxn_manip.add_reaction_output(r,df+i*N_per_leaf+c+j+1);
	  rxn_manip.add_reaction_rate(r,r_m_p.k_bubble);
	  rxns.push_back(r);
	  rxn_manip.reset_reaction(r);
	}
      c += db;
      
    }
}


// get the reaction model
std::vector<reaction> replication_model::get_reactions()
{
  return rxns;
}


// prepare the initiator requirements
void replication_model::prepare_init_requirements()
{
  init_requirements.clear();

  std::array<int,2> init_req;

  std::cout << "init_requirements" << std::endl;
  for (int i=r_m_p.bubble_min_fil; i<(N_binding+1); i++)
    {
      init_req[0] = i + N_bubble;
      init_req[1] = i;
      init_requirements.push_back(init_req);
      std::cout << init_req[0] << ","
		<< init_req[1] << std::endl;
    }
  
}


// get the initiation requirements
std::vector<std::array<int,2>> replication_model::get_init_requirements()
{
  return init_requirements;
}
