#include <LAMMPS_sys/LAMMPS_simulator.hpp>

// constructor
LAMMPS_simulator::LAMMPS_simulator()
{
  sim_MPI_initialized = 0;
  sim_MPI_finalized = 0;
  sim_MPI_size = 0;
  sim_MPI_rank = 0;
  lmp = nullptr;
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

  cout << logfile << endl;
  
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
  lmp_sys->write_mono_xyz("/home/ben/Workspace/btree_chromo/test_case/test_mono_internal.xyz");
  lmp->input->one(("include " + filename).c_str());
}


// dump the simulation state to the system state
void LAMMPS_simulator::sim_to_sys()
{
  int N = lmp_sys->get_N_total();

  if (N > 0)
    {
      cout << "N_atoms = " << N << endl;
      double *coords = new double[3*N];
  
      lammps_gather_atoms(lmp, const_cast<char*>("x"), 1, 3, coords);

      cout << "finished gathering atoms" << endl;
      lmp_sys->set_coords_arr_total(coords,"row");

      // need to repeat procedure with quaternions for ellipsoids

      // need to separate atoms into mono_atoms, ribo_atoms, and bdry_atoms

      // need to separate ellipsoids into mono_ellipsoids and ribo_ellipsoids
      
      lmp_sys->write_mono_xyz("/home/ben/Workspace/btree_chromo/test_case/test_mono_internal2.xyz");

      delete[] coords;
      
    }
}
