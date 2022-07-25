#include <rep_kinetics/gillespie_solver.hpp>

// constructor
gillespie_solver::gillespie_solver()
{
  N = -1;
  M = -1;
  S = nullptr;
  x = nullptr;
  xFPT = nullptr;
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

void gillespie_solver::initialize_x()
{
  destroy_x();
  
  if (N > 0)
    {
      
      x = new int[N];

      for (int i=0; i<N; i++)
	{
	  x[i] = 0;
	}
      
    }
  
}

void gillespie_solver::destroy_x()
{
  if (x != nullptr)
    {
      delete[] x;
      xFPT = nullptr;
    }
}

void gillespie_solver::initialize_xFPT()
{
  destroy_xFPT();
  
  if (N > 0)
    {
      
      xFPT = new int[N];

      for (int i=0; i<N; i++)
	{
	  x[i] = 0;
	}
      
    }
  
}

void gillespie_solver::destroy_xFPT()
{
  if (xFPT != nullptr)
    {
      delete[] xFPT;
      xFPT = nullptr;
    }
}

void gillespie_solver::initialize_S()
{
  destroy_S();
  if ((N > 0) &&
      (M > 0))
    {
      
      S = new int*[M];
      for (int j=0; j<M; j++)
	{
	  S[j] = new int[N];
	}

      for (int j=0; j<M; j++)
	{
	  for (int i=0; i<N; i++)
	    {
	      S[j][i] = 0;
	    }
	}
      
    }
}

void gillespie_solver::destroy_S()
{
  if (S != nullptr)
    {
      for (int j=0; j<M; j++)
	{
	  delete[] S[j];
	}
      delete[] S;
      S = nullptr;
    }
}

void gillespie_solver::set_N(int N_species)
{
  if (N_species > 0)
    {
      N = N_species;
    }
}

void gillespie_solver::set_M(int N_rxns)
{
  if (N_rxns > 0)
    {
      M = N_rxns;
    }
}

void gillespie_solver::prepare_reaction_system(int N_species, int N_rxns)
{
  cout << "setting M and N" << endl;
  set_N(N_species);
  set_M(N_rxns);
  cout << "initializing x" << endl;
  initialize_x();
  cout << "initializing xFPT" << endl;
  initialize_xFPT();
  cout << "initializing S" << endl;
  initialize_S();
}

void gillespie_solver::print_reaction_system()
{
  
  cout << "x, state vector" << endl;
  for (int i=0; i<N-1; i++)
    {
      cout << x[i] << ",";
    }
  cout << x[N-1] << endl;
  
  cout << "S, stoichiometric matrix" << endl;
  for (int j=0; j<M-1; j++)
    {
      for (int i=0; i<N-1; i++)
	{
	  cout << S[j][i] << ",";
	}
      cout << S[j][N-1] << endl;
    }
  for (int i=0; i<N-1; i++)
    {
      cout << S[M-1][i] << ",";
    }
  cout << S[M-1][N-1] << endl;
}


// select a reaction based on the propensities
// int gillespie_solver::select_rxn(double r_rxn, double total_propensity)
// {

// }

// update the propensities based on the current system state
// void gillespie_solver::update_propensities()
// {

// }

// update the system state based on the chosen reaction
void gillespie_solver::update_state(int j_rxn)
{
  for (int i=0; i<N; i++)
    {
      x[i] += S[j_rxn][i];
    }
}

// run the system until max time or first-passage occurs
// void gillespie_solver::run_FPT(double t, double t_max)
// {
//   double r_t, r_rxn;
//   double dt, total_propensity;
//   int FPT_index = -1;
//   int j_rxn;

//   r_rxn = 0.0;
//   total_propensity = 0.0;
//   j_rxn = 0;
  
//   while (1)
//     {
//       // sample random numbers for the reaction and time
//       r_t = u_rand(rand_eng);
//       r_rxn = u_rand(rand_eng);

//       // update the propensities
//       update_propensities();

//       // calculate the total propensity
      

//       // sample time based on total propensity
//       dt = -log(r_t)/total_propensity;
      
//       // test of proposed time exceeds maximum
//       if (t+dt > t_max)
// 	{
// 	  t = t_max;
// 	  break;
// 	}

//       // sample the reactions
//       j_rxn = select_rxn(r_rxn,total_propensity);
      
//       // update the state
//       update_state(j_rxn);

//       // test for FPT species and break the loop if they are found

//       if (FPT_index != -1)
// 	{
// 	  t += dt;
// 	  break;
// 	}
      
//     } // end while loop
// }



