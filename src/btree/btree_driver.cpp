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

  if (!drctvs_file)
    {
      cout << "ERROR: file not opened in read_directives" << endl;
    }
  else
    {
      while (1)
	{
	  drctvs_file >> line;
	  if (drctvs_file.eof()) break;
	  

	  if ((line.length() > 0) &&
	      (line.find("#") != 0))
	    {
	      drctvs.push_back(line);
	    }
	      
     	}
  
    }
  
  drctvs_file.close();
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

  reqs.topo_update = true;
  reqs.CG_update = true;
  reqs.regions_present = false;
  reqs.rep_model_present = false;

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

	  cout << "COMMAND: " << command << endl;
	  for (long unsigned int i=0; i<params.size(); i++)
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

      

      // read input state from a file
      if (command == "input_state")
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

      
      // seed the PRNG
      else if (command == "prng_seed")
	{
	  error_code = prng_seed(params);
	}

      
      // print the current state
      else if (command == "print")
	{
	  error_code = print_state(reqs);
	}

      if (error_code != 0) return 1;

    } // end loop over directives
  cout << "---END EXECUTING DIRECTIVES---\n" << endl;
  
  return 0;
}

////////////////////////////////////////////
// set of functions to perform directives //
////////////////////////////////////////////



int btree_driver::input_state(vector<string> &params, drctv_reqs &reqs)
{
  if (params.size() != 1)
    {
      cout << "ERROR: wrong number of parameters, correct input file" << endl;
      return 1;
    }
  driver_st = driver_bt.read_state(params[0]);
  driver_bt.prepare_state(driver_st);
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
  driver_rep_model = driver_replicator.read_rep_model(params[0]);
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
  driver_replicator.run_replicate_FPT(driver_rep_model,
				      driver_bt.count_active_forks(),
				      stod(params[0]));
  return 0;
}


int btree_driver::prng_seed(vector<string> &params)
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
