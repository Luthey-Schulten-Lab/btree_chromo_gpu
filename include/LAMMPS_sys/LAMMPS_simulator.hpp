#ifndef INCLUDE_LAMMPS_SIMULATOR_HPP
#define INCLUDE_LAMMPS_SIMULATOR_HPP

// standard library include files
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <memory>

// OpenMPI include files
#include <mpi.h>

// LAMMPS include files
#include <lammps.h>
#include "domain.h"
#include "group.h"
#include "input.h"
#include "atom.h"
#include "library.h"

// btree_chromo include files
#include <LAMMPS_sys/LAMMPS_sys.hpp>

using namespace std;

class LAMMPS_simulator
{
public:

  // constructor and destructor
  LAMMPS_simulator();
  ~LAMMPS_simulator();

  void LAMMPS_initialize(string logfile);
  void LAMMPS_destroy();

  void include_file(string filename);
  
private:

  int sim_MPI_initialized, sim_MPI_finalized;
  int sim_MPI_size; // MPI size
  int sim_MPI_rank; // current MPI rank

  LAMMPS_NS::LAMMPS *lmp;
  

};

#endif
