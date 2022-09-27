#include <LAMMPS_sys/LAMMPS_simulator.hpp>

// constructor
LAMMPS_simulator::LAMMPS_simulator()
{
  sim_MPI_initialized = 0;
  sim_MPI_finalized = 0;
  sim_MPI_size = 0;
  sim_MPI_rank = 0;
  lmp = nullptr;

  nProc = 1;
  prng_seed = 0;
}


// destructor
LAMMPS_simulator::~LAMMPS_simulator()
{
  // destroy the LAMMPS object
  LAMMPS_destroy();
  
  // finalize MPI instance
  MPI_Finalized(&sim_MPI_finalized);
  if (!sim_MPI_finalized) MPI_Finalize();
}


// initialize LAMMPS object
void LAMMPS_simulator::LAMMPS_initialize(string logfile)
{
  
  // int argc;
  // char **argv;
  // argc = 0;
  // argv = nullptr;
  
  // set up MPI instance
  // MPI_Init(&argc,&argv);
  MPI_Initialized(&sim_MPI_initialized);
  if (!sim_MPI_initialized) MPI_Init(nullptr,nullptr);
  MPI_Comm_size(MPI_COMM_WORLD,&sim_MPI_size); // MPI size
  MPI_Comm_rank(MPI_COMM_WORLD,&sim_MPI_rank); // current MPI rank

  // custom argument vector for LAMMPS library
  const char *lmpargv[] {"liblammps", "-log", logfile.c_str()};
  int lmpargc = sizeof(lmpargv)/sizeof(const char *);

  lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv, MPI_COMM_WORLD);
  
}


// destroy LAMMPS object
void LAMMPS_simulator::LAMMPS_destroy()
{
  if (lmp != nullptr)
    {
      delete lmp;
      lmp = nullptr;
    }
}


// set the lmp_sys object
void LAMMPS_simulator::set_lmp_sys(LAMMPS_sys *lmp_sys)
{
  this->lmp_sys = lmp_sys;
}


// feed an include file to LAMMPS simulation object
void LAMMPS_simulator::include_file(string filename)
{
  lmp->input->one(("include " + filename).c_str());
}


// read data into LAMMPS simulation object
void LAMMPS_simulator::read_data(string data_file)
{
  lmp->input->one(("read_data " + data_file).c_str());
}


// run the global setup with packages, units, atom_style, boundary, and atom_modify
void LAMMPS_simulator::global_setup()
{
  lmp->input->one("include ${DNA_model_dir}/protocol_subroutines/subroutine.global_setup");
}


// clear the LAMMPS system state
void LAMMPS_simulator::clear()
{
  lmp->input->one("clear");
}


// set number of processors
void LAMMPS_simulator::set_nProc(int nProc)
{
  // variable with number of processors for OpenMP
  this->nProc = nProc;
  lmp->input->one(("variable nProc internal " + to_string(this->nProc)).c_str());
}


// set DNA model
void LAMMPS_simulator::set_DNA_model_dir(string DNA_model_dir)
{
  // directory with DNA model files (properties, parameters, and basic routines)
  this->DNA_model_dir = DNA_model_dir;
  lmp->input->one(("variable DNA_model_dir string " + this->DNA_model_dir).c_str());
}


// set output details
void LAMMPS_simulator::set_output_details(string output_dir, string output_file_label)
{
  // output directory and file label
  this->output_dir = output_dir;
  this->output_file_label = output_file_label;
  lmp->input->one(("variable output_dir string " + this->output_dir).c_str());
  lmp->input->one(("variable output_file_label string " + this->output_file_label).c_str());
}


// set the PRNG seed
void LAMMPS_simulator::set_prng_seed(int s)
{
  // PRNG seed for LAMMPS object
  this->prng_seed = s;
  string old_str = to_string(this->prng_seed);
  auto seed_str = string(6 - min(6,static_cast<int>(old_str.length())),'0') + old_str;
  lmp->input->one(("variable rng_seed internal " + seed_str).c_str());
}


// set the timestep
void LAMMPS_simulator::set_delta_t(double delta_t)
{
  // timestep
  this->delta_t = delta_t;
  lmp->input->one(("variable delta_t internal " + to_string(delta_t)).c_str());
}


// reset variables specifying the simulation protocol
void LAMMPS_simulator::reset_protocol_variables()
{
  // variable with number of processors for OpenMP
  lmp->input->one(("variable nProc internal " + to_string(nProc)).c_str());

  // directory with DNA model files (properties, parameters, and basic routines)
  lmp->input->one(("variable DNA_model_dir string " + DNA_model_dir).c_str());

  // output files
  lmp->input->one(("variable output_dir string " + output_dir).c_str());
  lmp->input->one(("variable output_file_label string " + output_file_label).c_str());

  // PRNG seed for LAMMPS object
  string old_str = to_string(prng_seed);
  auto seed_str = string(6 - min(6,static_cast<int>(old_str.length())),'0') + old_str;
  lmp->input->one(("variable rng_seed internal " + seed_str).c_str());

  // timestep
  lmp->input->one(("variable delta_t internal " + to_string(delta_t)).c_str());
}


// minimize with soft potentials and harmonic bonds
void LAMMPS_simulator::minimize_soft_harmonic(thermo_dump_parameters t_d_p)
{
  lmp->input->one("include ${DNA_model_dir}/minimize_subroutines/subroutine.min_soft_harmonic");
}


// dump the simulation state to the system state
void LAMMPS_simulator::sim_to_sys()
{
  int N = lmp_sys->get_N_total();
  int N_mono = lmp_sys->get_N_mono();
  int N_ribo = lmp_sys->get_N_ribo();
  // int N_bdry = lmp_sys->get_N_bdry();

  int N_mono_ribo = N_mono + N_ribo;

  if (N > 0)
    {

      // copy the coordinates to the system state
      
      double *coords = new double[3*N];

      // 1 for per-atom type, 3 for size of data (x_i,x_j,x_k)
      lammps_gather_atoms(lmp, const_cast<char*>("x"), 1, 3, coords);

      lmp_sys->set_coords_arr_total(coords,"row");

      delete[] coords;

      // copy the quaternions to the system state

      // get the quats from a compute
      void *quats_p;
      // 1 for LMP_STYLE_ATOM, 2 for LMP_TYPE_ARRAY
      quats_p = lammps_extract_compute(lmp,const_cast<char*>("quat"),1,2);
      double **quats_2d{static_cast<double**>(quats_p)};

      // get the ids to match up the quats
      void *ids_p;
      // 1 for LMP_STYLE_ATOM, 1 for LMP_TYPE_VECTOR
      ids_p = lammps_extract_compute(lmp,const_cast<char*>("id_track"),1,1);
      double *ids{static_cast<double*>(ids_p)};

      double *quats = new double[4*N_mono_ribo];

      int id;
      for (int i=0; i<N; i++)
	{
	  id = int(ids[i]);
	  if (id <= N_mono_ribo)
	    {
	      for (int j=0; j<4; j++)
		{
		  quats[4*(id-1)+j] = quats_2d[i][j];
		}
	    }
	}

      lmp_sys->set_quats_arr_total(quats,"row");

      delete[] quats;

      // sync the subarrays with the total array
      
      lmp_sys->sync_subarrays();
      
    }
}
