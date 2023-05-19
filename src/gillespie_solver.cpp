#include <gillespie_solver.hpp>

// constructor
gillespie_solver::gillespie_solver()
{
  this->N = -1;
  this->M = -1;
  this->x = nullptr;
  this->nreactants = nullptr;
  this->nproducts = nullptr;
  this->reactants = nullptr;
  this->products = nullptr;
  this->k = nullptr;
  this->W = nullptr;
  this->rand_eng.seed(0);
  this->u_rand = std::uniform_real_distribution<double>(0.0,1.0);
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


// reaction volume
void gillespie_solver::set_V(double V)
{
  this->V = V;
}


// set the state vector
void gillespie_solver::set_x(std::vector<species_count> &s_cs)
{
  for (species_count s_c : s_cs)
    {
      x[s_c.id] = s_c.N;
    }
}


// set the FPT vector
void gillespie_solver::set_xFPT(std::vector<species_count> s_cs)
{
  xFPT = s_cs;
}


// conver the system state to species counts
std::vector<species_count> gillespie_solver::state_to_sc()
{
  std::vector<species_count> s_cs;
  species_count s_c;

  for (int i=0; i<N; i++)
    {
      if (x[i] > 0)
	{
	  s_c.id = i;
	  s_c.N = x[i];
	  s_cs.push_back(s_c);
	}
    }
  return s_cs;
}


// initialize the reaction system
void gillespie_solver::initialize_reaction_system(int N, int M)
{
  destroy_reaction_system();
  
  this->N = N;
  this->M = M;

  if (N > 0)
    {
      x = new int[N];

      for (int i=0; i<N; i++)
	{
	  x[i] = 0;
	}
    }

  if (M > 0)
    {
      nreactants = new int[M];
      nproducts = new int[M];

      reactants = new int*[M];
      products = new int*[M];

      k = new double[M];
      W = new double[M];

      for (int j=0; j<M; j++)
	{
	  nreactants[j] = 0;
	  nproducts[j] = 0;
	  k[j] = 0.0;
	  W[j] = 0.0;
	  
	  reactants[j] = new int[2];
	  products[j] = new int[max_nproducts];

	  for (int i=0; i<2; i++)
	    {
	      reactants[j][i] = 0;
	    }
	  for (int i=0; i<max_nproducts; i++)
	    {
	      products[j][i] = 0;
	    }
	  
	}
      
    }

  
  
}


// destroy the reaction system
void gillespie_solver::destroy_reaction_system()
{
  if (x != nullptr)
    {
      delete[] x;
      x = nullptr;
    }
  
  if (nreactants != nullptr)
    {
      delete[] nreactants;
      nreactants = nullptr;
    }
  if (nproducts != nullptr)
    {
      delete[] nproducts;
      nproducts = nullptr;
    }

  if (reactants != nullptr)
    {
      for (int j=0; j<M; j++)
	{
	  delete[] reactants[j];
	}
      delete[] reactants;
      reactants = nullptr;
    }
  if (products != nullptr)
    {
      for (int j=0; j<M; j++)
	{
	  delete[] products[j];
	}
      delete[] products;
      products = nullptr;
    }
  
  if (k != nullptr)
    {
      delete[] k;
      k = nullptr;
    }
  if (W != nullptr)
    {
      delete[] W;
      W = nullptr;
    }
}


// select a reaction based on the propensities
int gillespie_solver::select_rxn(double r_rxn, double total_propensity)
{
  int j;
  double utp = r_rxn*total_propensity;
  double p;

  j = 0;
  p = W[j];
  while (p <= utp)
    {
      j += 1;
      p += W[j];
    }
  return j;
}

// test the xFPT
int gillespie_solver::test_FPT()
{
  for (size_t i=0; i<xFPT.size(); i++)
    {
      if (x[xFPT[i].id] == xFPT[i].N) return xFPT[i].id;
    }
  return -1;
}

// update the system state based on the chosen reaction
void gillespie_solver::update_state(int j_rxn)
{
  // subtract the reactants
  for (int i=0; i<nreactants[j_rxn]; i++)
    {
      x[reactants[j_rxn][i]] -= 1;
    }
  // add the products
  for (int i=0; i<nproducts[j_rxn]; i++)
    {
      x[products[j_rxn][i]] += 1;
    }
}


// update the propensities
void gillespie_solver::update_propensities()
{
  
}



// calculate the total propensity
double gillespie_solver::calc_total_propensity()
{
  double t_p = 0.0;
  for (int j=0; j<M; j++)
    {
      t_p += W[j];
    }
  return t_p;
}


// run the system until max time or first-passage occurs
void gillespie_solver::run_FPT(double t,
			       double &dt,
			       double &dt_max)
{

  // double t_V;
  double r_t, r_rxn;
  double ds, total_propensity;
  int FPT_index = -1;
  int j_rxn;

  dt = 0.0;
  
  r_rxn = 0.0;
  total_propensity = 0.0;
  j_rxn = 0;
  
  while (1)
    {

      // std::cout << "t=" << t << std::endl;
      // sample random numbers for the reaction and time
      r_t = u_rand(rand_eng);
      r_rxn = u_rand(rand_eng);

      // calculate time-dependent volume
      // t_V = t + dt;
      
      // update the propensities
      update_propensities();

      // calculate the total propensity
      total_propensity = calc_total_propensity();

      // sample time based on total propensity
      ds = -log(r_t)/total_propensity;
      // std::cout << "dt=" << dt << std::endl;
      
      // test of proposed time exceeds maximum
      if (dt + ds > dt_max)
	{
	  dt = dt_max;
	  break;
	}
      else
	{
	  dt += ds;
	}

      // sample the reactions
      j_rxn = select_rxn(r_rxn,total_propensity);
      
      // update the state
      update_state(j_rxn);

      // test for FPT species and break the loop if they are found
      FPT_index = test_FPT();

      // std::cout << "FPT_index = " << FPT_index << std::endl;

      if (FPT_index != -1) break;
      
    } // end while loop

  std::cout << "\nxf, final state vector at t = " << t << " (s)" << std::endl;
  for (int i=0; i<N-1; i++)
    {
      std::cout << x[i] << ",";
    }
  std::cout << x[N-1] << std::endl;
  std::cout << "\n" << std::endl;
}



