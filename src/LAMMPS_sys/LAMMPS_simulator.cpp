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
  Nt = 0;
  stored_Nt = 0;
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


// feed a single command to the LAMMPS simulation object
void LAMMPS_simulator::command(string command)
{
  lmp->input->one(command.c_str());
}


// read data into LAMMPS simulation object
void LAMMPS_simulator::read_data(string data_file)
{
  // read the data
  lmp->input->one(("read_data " + data_file + " extra/bond/per/atom 1").c_str());

  // include the physical parameterization of the DNA polymer model
  lmp->input->one("include ${DNA_model_dir}/lmp.DNA_physical_params");

  // fix the boundary particles in a static position
  lmp->input->one("include ${DNA_model_dir}/potentials/lmp.bdry_static");
}


// run the global setup with packages, units, atom_style, boundary, and atom_modify
void LAMMPS_simulator::global_setup()
{

  computes_active.quats = false;
  computes_active.ids = false;
  computes_active.types = false;
  computes_active.MSD = false;
  dumps_active.lammpstrj = false;
  T_freq_specified = false;
  D_freq_specified = false;
  
  lmp->input->one("include ${DNA_model_dir}/protocol_subroutines/subroutine.global_setup");
}


// initialize the standard computes
void LAMMPS_simulator::standard_computes()
{
  // include compute for ids
  compute_trigger("ids");

  // include compute for types
  compute_trigger("types");
  
  // include compute for quats
  compute_trigger("quats");
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
  lmp->input->one("variable output_file string ${output_dir}${output_file_label}");
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
// void LAMMPS_simulator::set_delta_t(string delta_t)
void LAMMPS_simulator::set_delta_t(double delta_t)
{
  // timestep
  this->delta_t = delta_t;
  lmp->input->one(("variable delta_t internal " + to_string(delta_t)).c_str());
  // lmp->input->one(("variable delta_t internal " + delta_t).c_str());
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
  lmp->input->one("variable output_file string ${output_dir}${output_file_label}");

  // PRNG seed for LAMMPS object
  string old_str = to_string(prng_seed);
  auto seed_str = string(6 - min(6,static_cast<int>(old_str.length())),'0') + old_str;
  lmp->input->one(("variable rng_seed internal " + seed_str).c_str());

  // timestep
  lmp->input->one(("variable delta_t internal " + to_string(delta_t)).c_str());
  // lmp->input->one(("variable delta_t internal " + delta_t).c_str());
}


// setup for minimize routines
void LAMMPS_simulator::setup_minimize(thermo_dump_parameters &t_d_p)
{

  // set thermo frequency
  set_T_freq(t_d_p.thermo_freq);

  // set dump frequency
  set_D_freq(t_d_p.dump_freq);

  // include compute for ids
  compute_trigger("ids");

  // include compute for quats
  compute_trigger("quats");

  // include default thermo
  lmp->input->one("include ${DNA_model_dir}/dump_subroutines/subroutine.default_thermo");

  // include dump
  prepare_dump(t_d_p);
  
}


// minimize with soft potentials and harmonic bonds
void LAMMPS_simulator::minimize_soft_harmonic(thermo_dump_parameters t_d_p)
{

  // cout << "---[ minimizing SOFT_HARMONIC ]---" << endl;

  // setup the minimization
  setup_minimize(t_d_p);

  // include minimization subroutine
  lmp->input->one("include ${DNA_model_dir}/minimize_subroutines/subroutine.min_soft_harmonic");

  // reset the number of timesteps to Nt
  reset_timestep_to_Nt();
}



// minimize with hard potentials and harmonic bonds
void LAMMPS_simulator::minimize_hard_harmonic(thermo_dump_parameters t_d_p)
{

  // cout << "---[ minimizing HARD_HARMONIC ]---" << endl;
  
  // setup the minimization
  setup_minimize(t_d_p);

  // include minimization subroutine
  lmp->input->one("include ${DNA_model_dir}/minimize_subroutines/subroutine.min_hard_harmonic");

  // reset the number of timesteps to Nt
  reset_timestep_to_Nt();
}


// minimize with soft (topoisomerase) potentials and harmonic bonds
void LAMMPS_simulator::minimize_topoDNA_harmonic(thermo_dump_parameters t_d_p)
{

  // cout << "---[ minimizing topoDNA_HARMONIC ]---" << endl;
  
  // setup the minimization
  setup_minimize(t_d_p);

  // include minimization subroutine
  lmp->input->one("include ${DNA_model_dir}/minimize_subroutines/subroutine.min_topoDNA_harmonic");

  // reset the number of timesteps to Nt
  reset_timestep_to_Nt();
}

 
// minimize with soft potentials and FENE bonds
void LAMMPS_simulator::minimize_soft_FENE(thermo_dump_parameters t_d_p)
{

  // cout << "---[ minimizing SOFT_FENE ]---" << endl;
  
  // setup the minimization
  setup_minimize(t_d_p);

  // include minimization subroutine
  lmp->input->one("include ${DNA_model_dir}/minimize_subroutines/subroutine.min_soft_FENE");

  // reset the number of timesteps to Nt
  reset_timestep_to_Nt();
}


// minimize with hard potentials and FENE bonds
void LAMMPS_simulator::minimize_hard_FENE(thermo_dump_parameters t_d_p)
{

  // cout << "---[ minimizing HARD_FENE ]---" << endl;
  
  // setup the minimization
  setup_minimize(t_d_p);

  // include minimization subroutine
  lmp->input->one("include ${DNA_model_dir}/minimize_subroutines/subroutine.min_hard_FENE");

  // reset the number of timesteps to Nt
  reset_timestep_to_Nt();
}


// minimize with soft (topoisomerase) potentials and FENE bonds
void LAMMPS_simulator::minimize_topoDNA_FENE(thermo_dump_parameters t_d_p)
{

  // cout << "---[ minimizing topoDNA_FENE ]---" << endl;
  
  // setup the minimization
  setup_minimize(t_d_p);

  // include minimization subroutine
  lmp->input->one("include ${DNA_model_dir}/minimize_subroutines/subroutine.min_topoDNA_FENE");

  // reset the number of timesteps to Nt
  reset_timestep_to_Nt();
}


// setup for run routines
void LAMMPS_simulator::setup_run(thermo_dump_parameters &t_d_p)
{
  // set thermo frequency
  set_T_freq(t_d_p.thermo_freq);

  // set dump frequency
  set_D_freq(t_d_p.dump_freq);

  // include compute for ids
  compute_trigger("ids");

  // include compute for quats
  compute_trigger("quats");

  // include compute for MSD
  compute_trigger("MSD");

  // include default thermo
  lmp->input->one("include ${DNA_model_dir}/dump_subroutines/subroutine.MSD_thermo");

  // include dump
  prepare_dump(t_d_p);
}


// run with soft potentials and harmonic bonds
void LAMMPS_simulator::run_soft_harmonic(unsigned long N_steps, thermo_dump_parameters t_d_p)
{
  // setup for run
  setup_run(t_d_p);

  // include run subroutine
  lmp->input->one("include ${DNA_model_dir}/run_subroutines/subroutine.run_soft_harmonic");

  // set the timestep
  lmp->input->one("timestep ${delta_t}");
  
  // run for N_steps
  lmp->input->one(("run " + to_string(N_steps)).c_str());

  // increment Nt
  Nt += N_steps;
}


// run with hard potentials and harmonic bonds
void LAMMPS_simulator::run_hard_harmonic(unsigned long N_steps, thermo_dump_parameters t_d_p)
{
  // setup for run
  setup_run(t_d_p);

  // include run subroutine
  lmp->input->one("include ${DNA_model_dir}/run_subroutines/subroutine.run_hard_harmonic");

  // set the timestep
  lmp->input->one("timestep ${delta_t}");

  // run for N_steps
  lmp->input->one(("run " + to_string(N_steps)).c_str());

  // increment Nt
  Nt += N_steps;
}


// run with soft (topoisomerase) potentials and harmonic bonds
void LAMMPS_simulator::run_topoDNA_harmonic(unsigned long N_steps, thermo_dump_parameters t_d_p)
{
  // setup for run
  setup_run(t_d_p);

  // include run subroutine
  lmp->input->one("include ${DNA_model_dir}/run_subroutines/subroutine.run_topoDNA_harmonic");

  // set the timestep
  lmp->input->one("timestep ${delta_t}");

  // run for N_steps
  lmp->input->one(("run " + to_string(N_steps)).c_str());

  // increment Nt
  Nt += N_steps;
}


// run with soft potentials and FENE bonds
void LAMMPS_simulator::run_soft_FENE(unsigned long N_steps, thermo_dump_parameters t_d_p)
{
  // setup for run
  setup_run(t_d_p);

  // include run subroutine
  lmp->input->one("include ${DNA_model_dir}/run_subroutines/subroutine.run_soft_FENE");

  // set the timestep
  lmp->input->one("timestep ${delta_t}");

  // run for N_steps
  lmp->input->one(("run " + to_string(N_steps)).c_str());

  // increment Nt
  Nt += N_steps;
}


// run with hard potentials and FENE bonds
void LAMMPS_simulator::run_hard_FENE(unsigned long N_steps, thermo_dump_parameters t_d_p)
{  
  // setup for run
  setup_run(t_d_p);

  // set the timestep
  lmp->input->one("timestep ${delta_t}");

  // include run subroutine
  lmp->input->one("include ${DNA_model_dir}/run_subroutines/subroutine.run_hard_FENE");

  // run for N_steps
  lmp->input->one(("run " + to_string(N_steps)).c_str());

  // increment Nt
  Nt += N_steps;
}


// run with soft (topoisomerase) potentials and FENE bonds
void LAMMPS_simulator::run_topoDNA_FENE(unsigned long N_steps, thermo_dump_parameters t_d_p)
{  
  // setup for run
  setup_run(t_d_p);

  // set the timestep
  lmp->input->one("timestep ${delta_t}");

  // include run subroutine
  lmp->input->one("include ${DNA_model_dir}/run_subroutines/subroutine.run_topoDNA_FENE");

  // run for N_steps
  lmp->input->one(("run " + to_string(N_steps)).c_str());

  // increment Nt
  Nt += N_steps;
}


// prepare the dump
void LAMMPS_simulator::prepare_dump(thermo_dump_parameters &t_d_p)
{
  if (dumps_active.lammpstrj == true)
    {
      // undump the lammpstrj
      lmp->input->one("undump dumplammpstrj");
      dumps_active.lammpstrj = false;
    }


  cout << "dump_freq = " << t_d_p.dump_freq << endl;
  lmp->input->one("print ${D_freq}");
  if ((dumps_active.lammpstrj == false) && (t_d_p.dump_freq > 0))
    {
      // reset the number of timesteps to Nt
      lmp->input->one(("reset_timestep " + to_string(Nt)).c_str());

      // lmp->input->one(("variable skip_condition equal \"step > "+ to_string(Nt) + "\"").c_str());
      lmp->input->one(("variable D_delay equal "+ to_string(Nt + t_d_p.dump_freq)).c_str());
  
      // include dump
      if (t_d_p.append == true)
	{
	  lmp->input->one("include ${DNA_model_dir}/dump_subroutines/subroutine.dump_append_nofirst");
	}
      else
	{
	  if (t_d_p.write_first == true)
	    {
	      lmp->input->one("include ${DNA_model_dir}/dump_subroutines/subroutine.dump_noappend_first");
	    }
	  else
	    {
	      lmp->input->one("include ${DNA_model_dir}/dump_subroutines/subroutine.dump_noappend_nofirst");
	    }
	}

      dumps_active.lammpstrj = true;
    }
}


// set the variable for the thermo frequency
void LAMMPS_simulator::set_T_freq(int T_freq)
{
  if (T_freq_specified == true)
    {
      lmp->input->one("variable T_freq delete");
    }
  else
    {
      T_freq_specified = true;
    }
  lmp->input->one(("variable T_freq internal " + to_string(T_freq)).c_str());
}


// set the variable for the dump frequency
void LAMMPS_simulator::set_D_freq(int D_freq)
{
  if (D_freq_specified == true)
    {
      lmp->input->one("variable D_freq delete");
    }
  else
    {
      D_freq_specified = true;
    }
  lmp->input->one(("variable D_freq internal " + to_string(D_freq)).c_str());
}


// reset the simulation timestep to Nt
void LAMMPS_simulator::reset_timestep_to_Nt()
{
  if (dumps_active.lammpstrj == true)
    {
      // undump the lammpstrj
      lmp->input->one("undump dumplammpstrj");
      dumps_active.lammpstrj = false;
    }
  // reset the number of timesteps to Nt
  lmp->input->one(("reset_timestep " + to_string(Nt)).c_str());
}


// reset the timestep counter, Nt
void LAMMPS_simulator::reset_Nt(int Nt)
{
  this->Nt = Nt;
  reset_timestep_to_Nt();
}


// store the timestep counter, Nt
void LAMMPS_simulator::store_Nt()
{
  stored_Nt = Nt;
}


// restore the timestep counter to the stored value, stored_Nt
void LAMMPS_simulator::restore_Nt()
{
  Nt = stored_Nt;
  reset_timestep_to_Nt();
}


// trigger for preventing reuse of computes
void LAMMPS_simulator::compute_trigger(string compute_label)
{
  if (compute_label == "ids")
    {
      if (computes_active.ids == false)
	{
	  lmp->input->one("include ${DNA_model_dir}/compute_subroutines/subroutine.compute_ids");
	  computes_active.ids = true;
	}
    }
  else if (compute_label == "types")
    {
      if (computes_active.types == false)
	{
	  lmp->input->one("include ${DNA_model_dir}/compute_subroutines/subroutine.compute_types");
	  computes_active.types = true;
	}
    }
  else if (compute_label == "quats")
    {
      if (computes_active.quats == false)
	{
	  lmp->input->one("include ${DNA_model_dir}/compute_subroutines/subroutine.compute_quats");
	  computes_active.quats = true;
	}
    }
  else if (compute_label == "MSD")
    {
      if (computes_active.MSD == false)
	{
	  lmp->input->one("include ${DNA_model_dir}/compute_subroutines/subroutine.compute_MSD");
	  computes_active.MSD = true;
	}
    }
}


// get atom counts from the simulator and resize the system
void LAMMPS_simulator::sim_to_sys_atom_counts()
{
  
  double Nd = lammps_get_natoms(lmp);
  int N = int(Nd);
  cout << "Nd = " << Nd << endl;
  cout << "N = " << N << endl;

  // determine the sizes of the subarrays
  // unsigned long int *types = new unsigned long int[N];

  // 1 for per-atom type, 1 for size of data (t)
  // lammps_gather_atoms(lmp, const_cast<char*>("type"), 1, 1, types);

  // get the types to match up the quats
  void *types_p;
  // 1 for LMP_STYLE_ATOM, 1 for LMP_TYPE_VECTOR
  types_p = lammps_extract_compute(lmp,const_cast<char*>("type_track"),1,1);
  double *types{static_cast<double*>(types_p)};

  int N_mono = 0;
  int N_ribo = 0;
  int N_bdry = 0;

  int t;
  for (int i=0; i<N; i++)
    {
      // t = static_cast<int>(types[i]);
      t = int(types[i]);
      if (i%50 == 0) cout << i << "\t" << t << endl;
      if (t == 1)
	{
	  N_bdry += 1;
	}
      else if (t == 2)
	{
	  N_ribo += 1;
	}
      else
	{
	  N_mono += 1;
	}
    }

  int N_mono_ribo = N_mono + N_ribo;
  cout << "N_mono = " << N_mono << endl;
  cout << "N_ribo = " << N_ribo << endl;
  cout << "N_bdry = " << N_bdry << endl;
  cout << "N_mono_ribo = " << N_mono_ribo << endl;

  delete[] types;
  cout << "past delete" << endl;

  // resize the system state based on the simulator
  lmp_sys->set_N_total(N);
  lmp_sys->set_N_mono(N_mono);
  lmp_sys->set_N_ribo(N_ribo);
  lmp_sys->set_N_bdry(N_bdry);
  lmp_sys->set_N_total_ellipsoids(N_mono_ribo);
  lmp_sys->set_N_mono_ellipsoids(N_mono);
  lmp_sys->set_N_ribo_ellipsoids(N_ribo);
  
}


// dump the simulation state to the system state
void LAMMPS_simulator::sim_to_sys()
{
  
  // sim_to_sys_atom_counts();
  
  int N = lmp_sys->get_N_total();
  int N_mono = lmp_sys->get_N_mono();
  int N_ribo = lmp_sys->get_N_ribo();
  int N_bdry = lmp_sys->get_N_bdry();
  int N_mono_ribo = N_mono + N_ribo;

  cout << "N_mono = " << N_mono << endl;
  cout << "N_ribo = " << N_ribo << endl;
  cout << "N_bdry = " << N_bdry << endl;
  cout << "N_mono_ribo = " << N_mono_ribo << endl;

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


// read the replication model
int LAMMPS_simulator::read_loop_params(string loop_param_filename)
{

  fstream loop_param_file;

  string param_delim, param, val;
  int delim;
  
  string line;

  loop_sim_params l_sim_p;
  loop_sys_params l_sys_p;

  param_delim = "=";

  loop_param_file.open(loop_param_filename, ios::in);

  if (!loop_param_file.is_open())
    {
      cout << "ERROR: file not opened in loop_param_file" << endl;
      return 1;
    }
  else
    {
      while (1)
	{
	  loop_param_file >> line;
	  if (loop_param_file.eof()) break;
	  

	  if ((line.length() > 0) &&
	      (line.find("#") != 0))
	    {

	      delim = line.find(param_delim);

	      if (delim != -1)
		{

		  param = line.substr(0,delim);
		  val = line.substr(delim+1,line.length());

		  // cout << param << "=" << val << endl;

		  if (param == "min_dist")
		    {
		      l_sys_p.min_dist = stoi(val);
		    }

		  else if (param == "family")
		    {
		      l_sys_p.family = val;
		    }

		  else if (param == "ext_avg")
		    {
		      l_sys_p.ext_avg = stod(val);
		    }

		  else if (param == "ext_max")
		    {
		      l_sys_p.ext_max = stoi(val);
		    }

		  else if (param == "p_unbinding")
		    {
		      l_sys_p.p_unbinding = stod(val);
		    }
		  
		  else if (param == "r_g")
		    {
		      l_sys_p.r_g = stod(val);
		    }

		  else if (param == "r_0")
		    {
		      l_sim_p.r_0 = stod(val);
		    }

		  else if (param == "k")
		    {
		      l_sim_p.k = stoi(val);
		    }

		  else if (param == "freq_loop")
		    {
		      l_sim_p.freq_loop = stoul(val);
		    }

		  else if (param == "freq_topo")
		    {
		      l_sim_p.freq_topo = stoul(val);
		    }

		  else if (param == "dNt_topo")
		    {
		      l_sim_p.dNt_topo = stoul(val);
		    }

		}
	      	      
	    }
	      
     	} // end while loop

      loop_param_file.close();

      set_loop_sim_params(l_sim_p);
      lmp_sys->set_loop_sys_params(l_sys_p);

      return 0;
  
    }

}


// set the loop sim parameters
void LAMMPS_simulator::set_loop_sim_params(loop_sim_params &l_sim_p)
{
  this->l_sim_p = l_sim_p;
}


// update the loop bonds
void LAMMPS_simulator::update_loop_bonds(bool new_bonds)
{

  // delete the existing loop bonds
  if (new_bonds == false)
    {
      lmp->input->one("delete_bonds DNA bond 2 remove");
    }

  // type array for scatter
  int N = lmp_sys->get_N_total();
  int *types = new int[N];

  // 1 for per-atom type, 1 for size of data (t)
  // lammps_gather_atoms(lmp, const_cast<char*>("type"), 1, 1, types);

  lmp_sys->get_types(types);

  if (new_bonds == false)
    {
      
      lmp_sys->update_loop_topo();
      
    }


  // get the loop bonds
  vector<bond> loop_bonds = lmp_sys->get_loop_bonds();
  
  // add the updated loop bonds
  string temp_bond_command = "create_bonds single/bond";
  string bond_command;
  for (size_t i_loop=0; i_loop<loop_bonds.size(); i_loop++)
    {
      bond_command = temp_bond_command;
      bond_command += (" " + to_string(loop_bonds[i_loop].type));
      bond_command += (" " + to_string(loop_bonds[i_loop].i));
      bond_command += (" " + to_string(loop_bonds[i_loop].j));

      if (i_loop == loop_bonds.size() - 1)
	{
	  bond_command += " special yes";
	}
      else
	{
	  bond_command += " special no";
	}

      cout << i_loop << " " << bond_command << endl;
      
      lmp->input->one(bond_command);

      if ((types[loop_bonds[i_loop].i-1] != 4) && (types[loop_bonds[i_loop].i-1] != 5))
	{
	  types[loop_bonds[i_loop].i-1] = 7; // anchor atom
	}
      if ((types[loop_bonds[i_loop].j-1] != 4) && (types[loop_bonds[i_loop].j-1] != 5))
	{
	  types[loop_bonds[i_loop].j-1] = 8; // hinge atom
	}
    }

  // scatter the now modified atom types
  // 0 for integer type, 1 for per-atom count
  lammps_scatter_atoms(lmp, const_cast<char*>("type"), 0, 1, types);
  
  delete[] types;
}


// run a system with loops
void LAMMPS_simulator::run_loops(int N_loops, unsigned long N_steps, thermo_dump_parameters t_d_p)
{

  unsigned long step_counter = 0;
  unsigned long step_increment;
  thermo_dump_parameters t_d_p_iter = t_d_p;
  thermo_dump_parameters t_d_p_topo = t_d_p;
  int Nt_pre_topo, step_prev_topo;

  sim_to_sys();
  lmp_sys->initialize_loop_topo(N_loops);

  t_d_p_iter.write_first = false;
  t_d_p_iter.append = true;
  t_d_p_topo.dump_freq = 0;

  step_prev_topo = step_counter;
  
  bool new_bonds = true;
  update_loop_bonds(new_bonds);

  // minimize the system
  minimize_hard_harmonic(t_d_p_topo);
  minimize_hard_FENE(t_d_p_topo);

  // run the system with hard pairs and FENE bonds
  step_increment = min(l_sim_p.freq_loop,N_steps-step_counter);

  // Nt_pre_topo = Nt;
  // run_hard_harmonic(step_increment/2,t_d_p_topo);
  // reset_Nt(Nt_pre_topo);

  // cout << "\n\nstep_counter = " << step_counter << endl;
  // cout << "step_increment = " << step_increment << "\n\n" << endl;
  
  run_hard_FENE(step_increment,t_d_p);

  step_counter += step_increment;

  new_bonds = false;

  while (step_counter < N_steps)
    {

      // change the bonds
      sim_to_sys();
      
      // update the loop bonds
      update_loop_bonds(new_bonds);

      // simulate topoisomerase action
      if (step_counter >= (l_sim_p.freq_topo + step_prev_topo))
	{
	  // store the timestep
	  Nt_pre_topo = Nt;
	  
	  // minimize with topoisomerase pair potentials
	  minimize_topoDNA_harmonic(t_d_p_topo);
	  minimize_topoDNA_FENE(t_d_p_topo);
	  run_topoDNA_FENE(l_sim_p.dNt_topo,t_d_p_topo);

	  // step the pair potentials back to full strength of hard pairs
	  minimize_soft_harmonic(t_d_p_topo);
	  minimize_soft_FENE(t_d_p_topo);
	  minimize_hard_harmonic(t_d_p_topo);
	  minimize_hard_FENE(t_d_p_topo);

	  // reset the timestep
	  reset_Nt(Nt_pre_topo);

	  // reset the previous topo step
	  step_prev_topo = step_counter;
	}
      else
	{
	  // minimize with hard pair potentials
	  minimize_hard_harmonic(t_d_p_topo);
	  minimize_hard_FENE(t_d_p_topo);
	}
      
      // run for loop freq
      step_increment = min(l_sim_p.freq_loop,N_steps-step_counter);

      // cout << "\n\nstep_counter = " << step_counter << endl;
      // cout << "step_increment = " << step_increment << "\n\n" << endl;

      // run the system
      run_hard_FENE(step_increment,t_d_p_iter);

      // advance the step counter
      step_counter += step_increment;
      
    }

}
