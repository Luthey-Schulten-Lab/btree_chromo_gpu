#include <btree/btree_driver.hpp>

// constructor
btree_driver::btree_driver()
{
  driver_bt.reset_root();
}


// destructor
btree_driver::~btree_driver()
{
  // cout << "btree_driver destructor after destroy" << endl;
  driver_bt.destroy_tree();
  driver_rg.clear();
  drctvs.clear();
  // cout << "btree_driver destructor after destroy" << endl;
} 


// function to read directives from file
void btree_driver::read_directives(string drctvs_filename)
{
  fstream drctvs_file;
  
  string line;

  drctvs_file.open(drctvs_filename, ios::in);

  cout << "\nREADING DIRECTIVES:\n" << endl;
  cout << "\t" << drctvs_filename << endl;

  if (!drctvs_file.is_open())
    {
      cout << "ERROR: file not opened in read_directives" << endl;
    }
  else
    {
      while (getline(drctvs_file,line))
	{
	  
	  if ((line.length() > 0) &&
	      (line.at(0) != '#'))
	    {
	      drctvs.push_back(line);
	    }
	      
     	}

      drctvs_file.close();
    }
  
}


// function print all directives
void btree_driver::print_directives()
{
  cout << "\n---PROGRAM DIRECTIVES---\n" << endl;
  for (string drctv: drctvs)
    {
      cout << drctv << endl;
    }
  cout << "\n------------------------\n" << endl;
}


// reset the command requirements
void btree_driver::reset_command_locks_and_updates()
{

  // locks
  lock_state["btree_initialized"] = false;
  lock_state["regions_present"] = false;
  lock_state["rep_model_present"] = false;
  lock_state["BD_lengths_present"] = false;
  lock_state["simulator_prepared"] = false;
  lock_state["DNA_model"] = false;
  lock_state["output_details"] = false;
  lock_state["delta_t"] = false;
  lock_state["lmp_data_present"] = false;
  lock_state["loop_params_present"] = false;

  // updates
  lock_state["topo_update"] = true;
  lock_state["CG_update"] = true;
  
}


// test the validity of the command given the lock state
int btree_driver::test_command_parameter_validity(string &command, vector<string> &params)
{
  
  // test the number of parameters
  if (N_param_reqs[command] != params.size())
    {
      cout << "ERROR: wrong number of parameters for (" << command << ")" << endl;
      cout << "\t" << params.size() << " were given but " << N_param_reqs[command] << " are required" << endl;
      return 1;
    }

  return 0;
}


// test the validity of the command given the lock state
int btree_driver::test_command_lock_validity(string &command)
{

  // test the lock state for compatability
  string key;
  bool s;
  for (size_t i_lock=0; i_lock<lock_tests[command].size(); i_lock++)
    {
      key = lock_tests[command][i_lock].key;
      s = lock_tests[command][i_lock].s;

      if (lock_state[key] != s)
	{
	  cout << "ERROR: incompatible lock state for (" << command << ")" << endl;
	  cout << "\t" << key << " = " << lock_state[key] << endl;
	  e += 1;
	}
    }

  // return the number of errors
  return e;
  
}


// update the locks after the completion of a command
void btree_driver::update_lock_state_post_command(string &command)
{
  // update the lock state given the command
  string key;
  bool s;
  for (size_t i_lock=0; i_lock<lock_updates[command].size(); i_lock++)
    {
      key = lock_updates[command][i_lock].key;
      s = lock_updates[command][i_lock].s;
      lock_state[key] = s;
    }
}


// create a new lock
lock btree_driver::new_lock(string key, bool s)
{
  lock l;
  l.key = key;
  l.s = s;
  return l;
}


// validate the command sequence lock states
int btree_driver::validate_command_sequence()
{
  int e = 0;

  size_t i = 0;

  while ((i < commands.size()) && (e == 0))
    {
      e += test_command_lock_validity(commands[i]);
      update_lock_state_post_command(commands[i]);
    }

  return e;
}


// validate the parameters in the command sequence
int btree_driver::validate_command_sequence_parameters()
{
  int e = 0;

  size_t i = 0;

  while ((i < commands.size()) && (e == 0))
    {
      e += test_command_lock_validity(commands[i],command_params[i]);
    }

  return e;
}


// prepare the set of lock tests
void btree_driver::prepare_command_requirements()
{
  lock t_l;
  vector<lock> t_ls;

  ////////////////////////////////
  // btree command requirements //
  ////////////////////////////////

  // new_chromo
  // number of required parameters
  N_param_reqs["new_chromo"] = 1;
  // lock tests
  t_ls.clear();
  lock_tests["new_chromo"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("topo_update",true));
  t_ls.push_back(new_lock("CG_update",true));
  lock_updates["new_chromo"] = t_ls;

  // input_state
  // number of required parameters
  N_param_reqs["input_state"] = 1;
  // lock tests
  t_ls.clear();
  lock_tests["input_state"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("topo_update",true));
  t_ls.push_back(new_lock("CG_update",true));
  lock_updates["input_state"] = t_ls;

  // output_state
  // number of required parameters
  N_param_reqs["output_state"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["output_state"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["output_state"] = t_ls;

  // transforms_file
  // number of required parameters
  N_param_reqs["transforms_file"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["transforms_file"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("topo_update",true));
  t_ls.push_back(new_lock("CG_update",true));
  lock_updates["transforms_file"] = t_ls;

  // transform
  // number of required parameters
  N_param_reqs["transform"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["transform"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("topo_update",true));
  t_ls.push_back(new_lock("CG_update",true));
  lock_updates["transform"] = t_ls;

  // random_transforms
  // number of required parameters
  N_param_reqs["random_transforms"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["random_transforms"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("topo_update",true));
  t_ls.push_back(new_lock("CG_update",true));
  lock_updates["random_transforms"] = t_ls;

  // regions_file
  // number of required parameters
  N_param_reqs["regions_file"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["regions_file"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("regions_present",true));
  lock_updates["regions_file"] = t_ls;

  // dump_regions
  // number of required parameters
  N_param_reqs["dump_regions"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  t_ls.push_back(new_lock("regions_present",true));
  lock_tests["dump_regions"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["regions_file"] = t_ls;

  // dump_topology
  // number of required parameters
  N_param_reqs["dump_topology"] = 2;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["dump_topology"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["dump_topology"] = t_ls;

  // update_topology
  // number of required parameters
  N_param_reqs["update_topology"] = 0;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["update_topology"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["update_topology"] = t_ls;

  // update_CG_map
  // number of required parameters
  N_param_reqs["update_CG_map"] = 1;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["update_CG_map"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["update_CG_map"] = t_ls;

  // dump_CG_map
  // number of required parameters
  N_param_reqs["dump_CG_map"] = 3;
  // lock tests
  t_ls.clear();
  t_ls.push_back(new_lock("btree_initialized",true));
  lock_tests["dump_CG_map"] = t_ls;
  // lock updates
  t_ls.clear();
  lock_updates["dump_CG_map"] = t_ls;

  
  /////////////////////////////////////
  // replicator command requirements //
  /////////////////////////////////////

  // load_rep_model
  // number of required parameters
  N_param_reqs["load_rep_model"] = 1;
  // lock tests
  t_ls.clear();
  lock_tests["load_rep_model"] = t_ls;
  // lock updates
  t_ls.clear();
  t_ls.push_back(new_lock("rep_model_present",true));
  lock_updates["load_rep_model"] = t_ls;
  
}


// function to execute directives
int btree_driver::execute_directives()
{

  string command, param, temp_params;
  vector<string> params;
  string cmd_delim = ":";
  string param_delim = ",";
  int delim;
  int error_code = 0;

  drctv_reqs reqs;



  cout << "\n---BEGIN EXECUTING DIRECTIVES---\n" << endl;

  // loop over directives
  for (string drctv: drctvs)
    {

      delim = drctv.find(cmd_delim);

      // get command and parameters
      if (delim != -1)
	{
	  
	  command = drctv.substr(0,delim);
	  temp_params = drctv.substr(delim+1,drctv.length());

	  params.clear();

	  delim = temp_params.find(param_delim);

	  while (delim != -1)
	    {
	      param = temp_params.substr(0,delim);
	      params.push_back(param);
	      temp_params.erase(0,delim+1);
	      delim = temp_params.find(param_delim);
	    }
	  params.push_back(temp_params);

	  cout << "\nCOMMAND: " << command << endl;
	  for (size_t i=0; i<params.size(); i++)
	    {
	      cout << "\tparam_" << i
		   << ": " << params[i]
		   << endl;
	    }

	  cout << "\n" << endl;
	  
	}
      // get a command without parameters
      else
	{
	  command = drctv;

	  cout << "COMMAND: " << command << endl;
	  cout << "\n" << endl;
	}

      

      //////////////////
      // COMMAND LIST //
      //////////////////


      ////////////////////////////
      // Binary Tree Directives //
      ////////////////////////////


      // seed the PRNG for the btree
      if (command == "btree_prng_seed")
	{
	  error_code = btree_prng_seed(params);
	}      

      
      // create a new chromosome
      else if (command == "new_chromo")
	{
	  error_code = new_chromo(params,reqs);
	}
      
      
      // read input state from a file
      else if (command == "input_state")
	{
	  error_code = input_state(params,reqs);
	}
      

      // write the state to an output file
      else if (command == "output_state")
	{
	  error_code = output_state(params);
	}

      
      // apply state transformations from a file
      else if (command == "transforms_file")
	{
	  error_code = transforms_file(params,reqs);
	}

      
      // apply state transformations from a file
      else if (command == "transform")
	{
	  error_code = transform(params,reqs);
	}

      
      // apply random transformations
      else if (command == "random_transforms")
	{
	  error_code = random_transforms(params,reqs);
	}

      
      // apply state transformations from a file
      else if (command == "regions_file")
	{
	  error_code = regions_file(params,reqs);
	}

      
      // apply state transformations from a file
      else if (command == "dump_regions")
	{
	  error_code = dump_regions(params,reqs);
	}

      
      // write the topology to an output file
      else if (command == "dump_topology")
	{
	  error_code = dump_topology(params,reqs);
	}

      
      // solve the topology of the current state
      else if (command == "update_topology")
	{
	  error_code = update_topology(reqs);
	}

      
      // update the CG map
      else if (command == "update_CG_map")
	{
	  error_code = update_CG_map(params,reqs);
	}

      
      // dump the CG map
      else if (command == "dump_CG_map")
	{
	  error_code = dump_CG_map(params,reqs);
	}


      // print the current state
      else if (command == "print")
	{
	  error_code = print_state(reqs);
	}

      
      ///////////////////////////
      // Replicator Directives //
      ///////////////////////////
 
      
      // read the replication model
      else if (command == "load_replication_model")
	{
	  error_code = load_rep_model(params,reqs);
	}


      // replicate until t_max
      else if (command == "replicate")
	{
	  error_code = replicate(params,reqs);
	}


      // seed the PRNG for the replicator
      else if (command == "replicator_prng_seed")
	{
	  error_code = replicator_prng_seed(params);
	}


      ///////////////////////////
      // LAMMPS_sys Directives //
      ///////////////////////////
      

      // load file containing the monomer coordinates
      else if (command == "load_mono_coords")
	{
	  error_code = load_mono_coords(params,reqs);
	}


      // load file containing the monomer quaternions
      else if (command == "load_mono_quats")
	{
	  error_code = load_mono_quats(params,reqs);
	}


      // load file containing the ribosome coordinates
      else if (command == "load_ribo_coords")
	{
	  error_code = load_ribo_coords(params,reqs);
	}


      // load file containing the ribosome quaternions
      else if (command == "load_ribo_quats")
	{
	  error_code = load_ribo_quats(params,reqs);
	}


      // load file containing the boundary coordinates
      else if (command == "load_bdry_coords")
	{
	  error_code = load_bdry_coords(params,reqs);
	}


      // write file containing the monomer coordinates
      else if (command == "write_mono_coords")
	{
	  error_code = write_mono_coords(params);
	}


      // write file containing the monomer quaternions
      else if (command == "write_mono_quats")
	{
	  error_code = write_mono_quats(params);
	}


      // write file containing the ribosome coordinates
      else if (command == "write_ribo_coords")
	{
	  error_code = write_ribo_coords(params);
	}


      // write file containing the ribosome quaternions
      else if (command == "write_ribo_quats")
	{
	  error_code = write_ribo_quats(params);
	}


      // write file containing the boundary coordinates
      else if (command == "write_bdry_coords")
	{
	  error_code = write_bdry_coords(params);
	}


      // write file containing the boundary coordinates
      else if (command == "spherical_bdry")
	{
	  error_code = spherical_bdry(params,reqs);
	}


      // load file containing length-scales for BD simulations
      else if (command == "load_BD_lengths")
	{
	  error_code = load_BD_lengths(params,reqs);
	}


      // switch bonds on/off
      else if (command == "switch_bonds")
	{
	  error_code = switch_bonds(params);
	}


      // switch bending angles on/off
      else if (command == "switch_bending_angles")
	{
	  error_code = switch_bending_angles(params);
	}


      // switch twisting angles on/off
      else if (command == "switch_twisting_angles")
	{
	  error_code = switch_twisting_angles(params);
	}


      // set initial state for mapper
      else if (command == "set_initial_state")
	{
	  error_code = set_initial_state(reqs);
	}


      // set final state for mapper
      else if (command == "set_final_state")
	{
	  error_code = set_final_state(reqs);
	}


      // map the replication based on the initial and final states
      else if (command == "map_replication")
	{
	  error_code = map_replication(reqs);
	}
      

      // write the LAMMPS system data
      else if (command == "write_LAMMPS_data")
	{
	  error_code = write_LAMMPS_data(params,reqs);
	}

      
      // write the LAMMPS system data
      else if (command == "write_mono_xyz")
	{
	  error_code = write_mono_xyz(params,reqs);
	}

      
      //////////////////////////
      // Simulator Directives //
      //////////////////////////
      

      // prepare the simulator
      else if (command == "prepare_simulator")
	{
	  error_code = prepare_simulator(params,reqs);
	}


      // sync the simulator and the system
      else if (command == "sync_simulator_and_system")
	{
	  error_code = sync_simulator_and_system(reqs);
	}


      // sync the simulator and the system
      else if (command == "clear_simulator")
	{
	  error_code = clear_simulator(reqs);
	}


      // run a file using the simulator
      else if (command == "simulator_include_file")
	{
	  error_code = simulator_include_file(params,reqs);
	}


      // read data (LAMMPS data.* file) into the simulator
      else if (command == "simulator_read_data")
	{
	  error_code = simulator_read_data(params,reqs);
	}


      // set the prng_seed for the simulator
      else if (command == "simulator_set_prng_seed")
	{
	  error_code = simulator_set_prng_seed(params,reqs);
	}


      // set the number of processors for the simulator
      else if (command == "simulator_set_nProc")
	{
	  error_code = simulator_set_nProc(params,reqs);
	}


      // set the DNA model for the simulator
      else if (command == "simulator_set_DNA_model")
	{
	  error_code = simulator_set_DNA_model(params,reqs);
	}


      // set the output details for the simulator
      else if (command == "simulator_set_output_details")
	{
	  error_code = simulator_set_output_details(params,reqs);
	}


      // set the timestep for the simulator
      else if (command == "simulator_set_delta_t")
	{
	  error_code = simulator_set_delta_t(params,reqs);
	}


      // store the simulator's timestep
      else if (command == "simulator_store_timestep")
	{
	  error_code = simulator_store_timestep(reqs);
	}


      // restore the simulator's timestep from the stored value
      else if (command == "simulator_restore_timestep")
	{
	  error_code = simulator_restore_timestep(reqs);
	}


      // minimize with soft potentials and harmonic bonds
      else if (command == "simulator_minimize_soft_harmonic")
	{
	  error_code = simulator_minimize<0,0>(params,reqs);
	}


      // minimize with hard potentials and harmonic bonds
      else if (command == "simulator_minimize_hard_harmonic")
	{
	  error_code = simulator_minimize<1,0>(params,reqs);
	}


      // minimize with soft potentials and FENE bonds
      else if (command == "simulator_minimize_soft_FENE")
	{
	  error_code = simulator_minimize<0,1>(params,reqs);
	}


      // minimize with hard potentials and FENE bonds
      else if (command == "simulator_minimize_hard_FENE")
	{
	  error_code = simulator_minimize<1,1>(params,reqs);
	}


      // run with soft potentials and harmonic bonds
      else if (command == "simulator_run_soft_harmonic")
	{
	  error_code = simulator_run<0,0>(params,reqs);
	}


      // run with hard potentials and harmonic bonds
      else if (command == "simulator_run_hard_harmonic")
	{
	  error_code = simulator_run<1,0>(params,reqs);
	}


      // run with soft potentials and FENE bonds
      else if (command == "simulator_run_soft_FENE")
	{
	  error_code = simulator_run<0,1>(params,reqs);
	}


      // run with hard potentials and FENE bonds
      else if (command == "simulator_run_hard_FENE")
	{
	  error_code = simulator_run<1,1>(params,reqs);
	}

      // load file containing loop parameters
      else if (command == "simulator_load_loop_params")
	{
	  error_code = simulator_load_loop_params(params,reqs);
	}


      // run a simulation with loops
      else if (command == "simulator_run_loops")
	{
	  error_code = simulator_run_loops(params,reqs);
	}


      if (error_code != 0) return 1;

    } // end loop over directives
  cout << "---END EXECUTING DIRECTIVES---\n" << endl;
  
  return 0;
}



////////////////////////////////////////////
// set of functions to perform directives //
////////////////////////////////////////////


int btree_driver::new_chromo(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  driver_st.size = stoi(params[0]);
  driver_st.transforms.clear();
  driver_bt.prepare_state(driver_st);
  // btree is now initialized
  reqs.btree_initialized = true;
  // require a topology update
  reqs.topo_update = true;
  // require a coarse-graining update
  reqs.CG_update = true;
  return 0;
}


int btree_driver::input_state(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  driver_st = driver_bt.read_state(params[0]);
  driver_bt.prepare_state(driver_st);
  // btree is now initialized
  reqs.btree_initialized = true;
  // require a topology update
  reqs.topo_update = true;
  // require a coarse-graining update
  reqs.CG_update = true;
  return 0;
}


int btree_driver::output_state(vector<string> &params)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  driver_bt.write_state(params[0],driver_bt.dump_state());
  return 0;
}


int btree_driver::transforms_file(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  driver_tr = driver_bt.read_transforms(params[0]);
  driver_bt.apply_transforms(driver_tr);
  // require a topology update
  reqs.topo_update = true;
  // require a coarse-graining update
  reqs.CG_update = true;
  return 0;
}


int btree_driver::transform(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  driver_bt.single_transform(driver_bt.parse_transform(params[0]));
  // require a topology update
  reqs.topo_update = true;
  // require a coarse-graining update
  reqs.CG_update = true;
  return 0;
}


int btree_driver::random_transforms(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  driver_bt.random_transforms(stoi(params[0]));
  // require a topology update
  reqs.topo_update = true;
  // require a coarse-graining update
  reqs.CG_update = true;
  return 0;
}


int btree_driver::regions_file(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 2)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  driver_rg.clear();
  driver_rg = driver_bt.read_regions(params[0],stoi(params[1]));
  // flag for regions being present
  reqs.regions_present = true;
  return 0;
}


int btree_driver::dump_regions(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.regions_present == true)
    {
	      
      // update topology before updating regions
      if (reqs.topo_update == true)
	{
	  driver_bt.solve_topology();
	  reqs.topo_update = false;
	}
 
      driver_bt.update_region_counts(driver_rg);
      driver_bt.dump_regions(params[0],driver_rg);
    }
  else
    {
      cout << "  missing (read_regions), no regions to dump" << endl;
    }
  return 0;
}


int btree_driver::dump_topology(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 2)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  // update topology before dumping
  if (reqs.topo_update == true)
    {
      driver_bt.solve_topology();
      reqs.topo_update = false;
    }
  driver_bt.dump_topology(params[0],stoi(params[1]));
  return 0;
}

int btree_driver::update_topology(drctv_reqs &reqs)
{
  driver_bt.solve_topology();
  // disable flag for topology updating after an update
  reqs.topo_update = false;
  return 0;
}


int btree_driver::update_CG_map(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  // update topology before dumping
  if (reqs.topo_update == true)
    {
      driver_bt.solve_topology();
      reqs.topo_update = false;
    }
  driver_CG = driver_bt.update_CG_map(stoi(params[0]));
  return 0;
}


int btree_driver::dump_CG_map(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 3)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  // update topology before updating CG_map
  if (reqs.topo_update == true)
    {
      driver_bt.solve_topology();
      reqs.topo_update = false;
    }
  // update CG_map before dumping
  if (reqs.CG_update == true)
    {
      driver_CG = driver_bt.update_CG_map(stoi(params[1]));
      reqs.CG_update = false;
    }
  driver_bt.dump_CG_map(params[0],stoi(params[2]),driver_CG);
  return 0;
}


int btree_driver::load_rep_model(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  driver_replicator.read_rep_model(params[0]);
  // replication model is now present
  reqs.rep_model_present = true;
  return 0;
}


int btree_driver::replicate(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  // test if replication model is present
  if (reqs.rep_model_present == false)
    {
      cout << "ERROR: missing replication model" << endl;
      return 1;
    }

  vector<init_loc> init_dist, new_leaves;
  string rep_leaf;
  init_loc i_l;
  double t, dt, t_max;
  int rep_amount, error_code, i_rep;
  int *noninit_s = nullptr;

  // create an initial distribution of initiator species with only free DnaA
  init_dist.clear();
  
  i_l.loc = "free";
  i_l.N = driver_replicator.get_N_init_DnaA();

  init_dist.push_back(i_l);
  
  for (string leaf : driver_bt.get_leaves())
    {
      i_l.loc = leaf;
      i_l.N = 0;

      init_dist.push_back(i_l);
    }

  // create the species counts of noninitiator species
  driver_replicator.reset_noninit_s(noninit_s);


  // set the current time and the maximum time
  t = 0.0;
  t_max = stod(params[0]);

  cout << "Performing random replications until t = " << t_max << "\n" << endl;

  while (t < t_max)
    {

      dt = t;

      // run Gillespie algorithm until an initiation event occurs
      driver_replicator.run_replicate_FPT(noninit_s,
					  init_dist,
					  t,
					  t_max);
      
      // calculate time difference
      dt = t - dt;

      // calculate amount of replicated DNA prior to new replication event
      // amount is proportional to time difference and number of active forks
      rep_amount = driver_replicator.get_k_rep()*min(2*driver_bt.count_active_forks(),
						     driver_replicator.get_max_replisomes())*dt;

      cout << "\n" << rep_amount << " units were replicated on active forks prior to event(/termination)" << endl;

      // perform random replications
      driver_bt.random_transforms(rep_amount);

      // clear the new leaves
      new_leaves.clear();
      
      // update the distribution to represent the new number of leaves and determine the branch leaf
      i_rep = 0;
      for (size_t i=0; i<init_dist.size(); i++)
	{
	  // cout << init_dist[i].loc << " = " << init_dist[i].N << endl;
	  if (init_dist[i].N == -1)
	    {
	      i_rep = i;
	      break;
	    }
	}

      if (i_rep > 0)
	{
	  cout << "\nreplication event at t = " << t << endl;
	  rep_leaf = init_dist[i_rep].loc;

	  // branch the btree at the replicating leaf
	  error_code = driver_bt.branch(rep_leaf);

	  if (error_code == 1)
	    {
	      return 1;
	    }
	  else
	    {
	      cout << "\tsplitting at initiated branch (" << rep_leaf << ") and updating initiator distribution\n" << endl;
	      init_dist.erase(init_dist.begin()+i_rep);

	      i_l.N = 0;
	      i_l.loc = rep_leaf + "l";
	      new_leaves.push_back(i_l);
	      i_l.loc = rep_leaf + "r";
	      new_leaves.push_back(i_l);

	      init_dist.insert(init_dist.begin()+i_rep,new_leaves.begin(),new_leaves.end());

	      for (init_loc temp_i_l : init_dist)
		{
		  cout << temp_i_l.loc << " = " << temp_i_l.N << endl;
		}

	    }
	  
	}
      
    }

  cout << "\nfinal initiator distribution" << endl;
  for (init_loc temp_i_l : init_dist)
    {
      cout << temp_i_l.loc << " = " << temp_i_l.N << endl;
    }
  cout << "\n" << endl;

  if (noninit_s != nullptr)
    {
      delete[] noninit_s;
      noninit_s = nullptr;
    }
  
  return 0;
}


int btree_driver::load_BD_lengths(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  driver_lmp_sys.read_BD_lengths(params[0]);
  // replication model is now present
  reqs.BD_lengths_present = true;
  return 0;
}


int btree_driver::load_mono_coords(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 2)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.btree_initialized == false)
    {
      cout << "ERROR: empty btree" << endl;
      return 1;
    }
  driver_lmp_sys.set_btree(driver_bt.dump_state());
  int e = driver_lmp_sys.read_mono_coords(params[0],params[1]);
  return e;
}


int btree_driver::load_mono_quats(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 2)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.btree_initialized == false)
    {
      cout << "ERROR: empty btree" << endl;
      return 1;
    }
  driver_lmp_sys.set_btree(driver_bt.dump_state());
  int e = driver_lmp_sys.read_mono_quats(params[0],params[1]);
  return e;
}


int btree_driver::load_ribo_coords(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 2)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.btree_initialized == false)
    {
      cout << "ERROR: empty btree" << endl;
      return 1;
    }
  int e = driver_lmp_sys.read_ribo_coords(params[0],params[1]);
  return e;
}


int btree_driver::load_ribo_quats(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 2)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.btree_initialized == false)
    {
      cout << "ERROR: empty btree" << endl;
      return 1;
    }
  int e = driver_lmp_sys.read_ribo_quats(params[0],params[1]);
  return e;
}


int btree_driver::load_bdry_coords(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 2)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.btree_initialized == false)
    {
      cout << "ERROR: empty btree" << endl;
      return 1;
    }
  int e = driver_lmp_sys.read_bdry_coords(params[0],params[1]);
  return e;
}


int btree_driver::write_mono_coords(vector<string> &params)
{
  if (params.size() != 2)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  int e = driver_lmp_sys.write_mono_coords(params[0],params[1]);
  return e;
}


int btree_driver::write_mono_quats(vector<string> &params)
{
  if (params.size() != 2)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  int e = driver_lmp_sys.write_mono_quats(params[0],params[1]);
  return e;
}


int btree_driver::write_ribo_coords(vector<string> &params)
{
  if (params.size() != 2)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  int e = driver_lmp_sys.write_ribo_coords(params[0],params[1]);
  return e;
}


int btree_driver::write_ribo_quats(vector<string> &params)
{
  if (params.size() != 2)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  int e = driver_lmp_sys.write_ribo_quats(params[0],params[1]);
  return e;
}


int btree_driver::write_bdry_coords(vector<string> &params)
{
  if (params.size() != 2)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  int e = driver_lmp_sys.write_bdry_coords(params[0],params[1]);
  return e;
}


int btree_driver::write_LAMMPS_data(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.btree_initialized == false)
    {
      cout << "ERROR: empty btree" << endl;
      return 1;
    }
  if (reqs.BD_lengths_present == false)
    {
      cout << "ERROR: missing BD lengths" << endl;
      return 1;
    }
  driver_lmp_sys.set_btree(driver_bt.dump_state());
  driver_lmp_sys.write_data(params[0]);
  return 0;
}


int btree_driver::spherical_bdry(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 4)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.BD_lengths_present == false)
    {
      cout << "ERROR: missing BD lengths" << endl;
      return 1;
    }
  driver_lmp_sys.generate_spherical_bdry(stod(params[0]),
					 stod(params[1]),
					 stod(params[2]),
					 stod(params[3]));
  return 0;
}


int btree_driver::switch_bonds(vector<string> &params)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }

  if (params[0] == "T")
    {
      driver_lmp_sys.switch_bonds(true);
    }
  else if (params[0] == "F")
    {
      driver_lmp_sys.switch_bonds(false);
    }
  else
    {
      cout << "ERROR: invalid switch" << endl;
      return 1;
    }
  return 0;
}


int btree_driver::switch_bending_angles(vector<string> &params)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }

  if (params[0] == "T")
    {
      driver_lmp_sys.switch_bending_angles(true);
    }
  else if (params[0] == "F")
    {
      driver_lmp_sys.switch_bending_angles(false);
    }
  else
    {
      cout << "ERROR: invalid switch" << endl;
      return 1;
    }
  return 0;
}


int btree_driver::switch_twisting_angles(vector<string> &params)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }

  if (params[0] == "T")
    {
      driver_lmp_sys.switch_twisting_angles(true);
    }
  else if (params[0] == "F")
    {
      driver_lmp_sys.switch_twisting_angles(false);
    }
  else
    {
      cout << "ERROR: invalid switch" << endl;
      return 1;
    }
  return 0;
}


int btree_driver::write_mono_xyz(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.btree_initialized == false)
    {
      cout << "ERROR: empty btree" << endl;
      return 1;
    }
  driver_lmp_sys.write_mono_xyz(params[0]);
  return 0;
}


int btree_driver::set_initial_state(drctv_reqs &reqs)
{
  if (reqs.btree_initialized == false)
    {
      cout << "ERROR: empty btree" << endl;
      return 1;
    }
  driver_mapper.set_initial_state(driver_bt.dump_state());
  reqs.map_initial_present = true;
  return 0;
}


int btree_driver::set_final_state(drctv_reqs &reqs)
{
  if (reqs.btree_initialized == false)
    {
      cout << "ERROR: empty btree" << endl;
      return 1;
    }
  driver_mapper.set_final_state(driver_bt.dump_state());
  reqs.map_final_present = true;
  return 0;
}


int btree_driver::map_replication(drctv_reqs &reqs)
{
  if (reqs.map_initial_present == false)
    {
      cout << "ERROR: missing initial state" << endl;
      return 1;
    }
  if (reqs.map_final_present == false)
    {
      cout << "ERROR: missing final state" << endl;
      return 1;
    }
  int e = driver_mapper.prepare_mapping();

  driver_lmp_sys.apply_mono_mapping(driver_mapper.get_map());
  return e;
}


int btree_driver::prepare_simulator(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  driver_lmp_simulator.LAMMPS_initialize(params[0]);
  driver_lmp_simulator.set_lmp_sys(&driver_lmp_sys);
  reqs.simulator_prepared = true;
  return 0;
}


int btree_driver::simulator_include_file(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.simulator_prepared == false)
    {
      cout << "ERROR: missing simulator" << endl;
      return 1;
    }
  driver_lmp_simulator.include_file(params[0]);
  return 0;
}


int btree_driver::sync_simulator_and_system(drctv_reqs &reqs)
{
  if (reqs.simulator_prepared == false)
    {
      cout << "ERROR: missing simulator" << endl;
      return 1;
    }
  driver_lmp_sys.set_btree(driver_bt.dump_state());
  driver_lmp_simulator.sim_to_sys();
  return 0;
}


int btree_driver::clear_simulator(drctv_reqs &reqs)
{
  if (reqs.simulator_prepared == false)
    {
      cout << "ERROR: missing simulator" << endl;
      return 1;
    }
  driver_lmp_simulator.clear();
  return 0;
}


int btree_driver::simulator_set_nProc(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.simulator_prepared == false)
    {
      cout << "ERROR: missing simulator" << endl;
      return 1;
    }
  driver_lmp_simulator.set_nProc(stoi(params[0]));  
  return 0;
}


int btree_driver::simulator_set_prng_seed(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.simulator_prepared == false)
    {
      cout << "ERROR: missing simulator" << endl;
      return 1;
    }
  driver_lmp_simulator.set_prng_seed(stoi(params[0]));
  driver_lmp_sys.prng_seed(stoi(params[0]));
  return 0;
}


int btree_driver::simulator_set_DNA_model(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.simulator_prepared == false)
    {
      cout << "ERROR: missing simulator" << endl;
      return 1;
    }
  driver_lmp_simulator.set_DNA_model_dir(params[0]);
  reqs.DNA_model = true;
  return 0;
}


int btree_driver::simulator_set_output_details(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 2)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.simulator_prepared == false)
    {
      cout << "ERROR: missing simulator" << endl;
      return 1;
    }
  driver_lmp_simulator.set_output_details(params[0],params[1]);
  reqs.output_details = true;
  return 0;
}


int btree_driver::simulator_set_delta_t(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  if (reqs.simulator_prepared == false)
    {
      cout << "ERROR: missing simulator" << endl;
      return 1;
    }
  driver_lmp_simulator.set_delta_t(stod(params[0]));
  // driver_lmp_simulator.set_delta_t(params[0]);
  reqs.delta_t = true;
  return 0;
}


int btree_driver::simulator_read_data(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  
  if (reqs.simulator_prepared == false)
    {
      cout << "ERROR: missing simulator" << endl;
      return 1;
    }
  else if (reqs.DNA_model == false)
    {
      cout << "ERROR: missing DNA model" << endl;
      return 1;
    }
  else if (reqs.output_details == false)
    {
      cout << "ERROR: missing output details" << endl;
      return 1;
    }
  else if (reqs.delta_t == false)
    {
      cout << "ERROR: missing timestep" << endl;
      return 1;
    }

  driver_lmp_simulator.clear();
  driver_lmp_simulator.reset_protocol_variables();
  driver_lmp_simulator.global_setup();
  driver_lmp_simulator.read_data(params[0]);
  driver_lmp_simulator.standard_computes();

  reqs.lmp_data_present = true;
  
  return 0;  
}


template<int SOFT_HARD, int HARMONIC_FENE>
int btree_driver::simulator_minimize(vector<string> &params, drctv_reqs &reqs)
{

  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }

  if (reqs.lmp_data_present == false)
    {
      cout << "ERROR: missing LAMMPS data for simulator" << endl;
      return 1;
    }

  thermo_dump_parameters t_d_p;

  t_d_p.append = false;
  t_d_p.write_first = true;
  t_d_p.dump_freq = 0;
  t_d_p.thermo_freq = stoi(params[0]);

  // run the minimization
  if (SOFT_HARD == 0)
    {
      if (HARMONIC_FENE == 0)
	{
	  driver_lmp_simulator.minimize_soft_harmonic(t_d_p);
	}
      else if (HARMONIC_FENE == 1)
	{
	  driver_lmp_simulator.minimize_soft_FENE(t_d_p);
	}
    }
  else if (SOFT_HARD == 1)
    {
      if (HARMONIC_FENE == 0)
	{
	  driver_lmp_simulator.minimize_hard_harmonic(t_d_p);
	}
      else if (HARMONIC_FENE == 1)
	{
	  driver_lmp_simulator.minimize_hard_FENE(t_d_p);
	}
    }

  return 0;
}


template<int SOFT_HARD, int HARMONIC_FENE>
int btree_driver::simulator_run(vector<string> &params, drctv_reqs &reqs)
{

  if (params.size() != 5)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }

  if (reqs.lmp_data_present == false)
    {
      cout << "ERROR: missing LAMMPS data for simulator" << endl;
      return 1;
    }

  thermo_dump_parameters t_d_p;

  t_d_p.append = false;
  t_d_p.write_first = true;
  if (params[3] == "append")
    {
      t_d_p.append = true;
      t_d_p.write_first = false;
    }
  if (params[4] == "skip_first") t_d_p.write_first = false;
  t_d_p.dump_freq = stoi(params[2]);
  t_d_p.thermo_freq = stoi(params[1]);

  // run the Brownian dynamics
  if (SOFT_HARD == 0)
    {
      if (HARMONIC_FENE == 0)
	{
	  driver_lmp_simulator.run_soft_harmonic(stoul(params[0]),t_d_p);
	}
      else if (HARMONIC_FENE == 1)
	{
	  driver_lmp_simulator.run_soft_FENE(stoul(params[0]),t_d_p);
	}
    }
  else if (SOFT_HARD == 1)
    {
      if (HARMONIC_FENE == 0)
	{
	  driver_lmp_simulator.run_hard_harmonic(stoul(params[0]),t_d_p);
	}
      else if (HARMONIC_FENE == 1)
	{
	  driver_lmp_simulator.run_hard_FENE(stoul(params[0]),t_d_p);
	}
    }

  return 0;
}


int btree_driver::simulator_load_loop_params(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }

  if (reqs.simulator_prepared == false)
    {
      cout << "ERROR: missing simulator" << endl;
      return 1;
    }

  int e = driver_lmp_simulator.read_loop_params(params[0]);
  reqs.loop_params_present = true;
  return e;
}


int btree_driver::simulator_run_loops(vector<string> &params, drctv_reqs &reqs)
{

  if (params.size() != 6)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }

  if (reqs.lmp_data_present == false)
    {
      cout << "ERROR: missing LAMMPS data for simulator" << endl;
      return 1;
    }

  if (reqs.loop_params_present == false)
    {
      cout << "ERROR: missing loop parameters for loop run" << endl;
      return 1;
    }

  thermo_dump_parameters t_d_p;

  t_d_p.append = false;
  t_d_p.write_first = true;
  if (params[4] == "append")
    {
      t_d_p.append = true;
      t_d_p.write_first = false;
    }
  if (params[5] == "skip_first") t_d_p.write_first = false;
  t_d_p.dump_freq = stoi(params[3]);
  t_d_p.thermo_freq = stoi(params[2]);

  driver_lmp_simulator.run_loops(stoi(params[0]),stoul(params[1]),t_d_p);

  return 0;
}


int btree_driver::simulator_store_timestep(drctv_reqs &reqs)
{

  if (reqs.simulator_prepared == false)
    {
      cout << "ERROR: missing simulator" << endl;
      return 1;
    }

  driver_lmp_simulator.store_Nt();
  return 0;
}


int btree_driver::simulator_restore_timestep(drctv_reqs &reqs)
{

  if (reqs.simulator_prepared == false)
    {
      cout << "ERROR: missing simulator" << endl;
      return 1;
    }

  driver_lmp_simulator.restore_Nt();
  return 0;
}


int btree_driver::btree_prng_seed(vector<string> &params)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  // seed the PRNG
  driver_bt.prng_seed(stoi(params[0]));
  return 0;
}


int btree_driver::replicator_prng_seed(vector<string> &params)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  // seed the PRNG
  driver_replicator.prng_seed(stoi(params[0]));
  return 0;
}


int btree_driver::print_state(drctv_reqs &reqs)
{
  // update topology before printing
  if (reqs.topo_update == true)
    {
      driver_bt.solve_topology();
      reqs.topo_update = false;
    }
  driver_bt.print_tree();
  return 0;
}
