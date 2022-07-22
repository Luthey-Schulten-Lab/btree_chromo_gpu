#include <rep_kinetics/gillespie_solver.hpp>

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

// int gillespie_solver::select_rxn(double r_rxn, double total_propensity)
// {

// }

// void gillespie_solver::update_propensities()
// {

// }

// void gillespie_solver::update_state(int j_rxn)
// {

// }

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



