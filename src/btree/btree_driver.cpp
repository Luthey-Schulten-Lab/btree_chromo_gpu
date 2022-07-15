#include <btree/btree_driver.hpp>

// constructor
btree_driver::btree_driver()
{
  reset_root();
}


// destructor
btree_driver::~btree_driver()
{
  // cout << "btree_driver destructor after destroy" << endl;
  destroy_tree();
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

  bool require_topo_update = true;
  bool require_CG_update = true;
  bool regions_present = false;

  cout << "\n---BEGIN EXECUTING DIRECTIVES---\n" << endl;

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
	  if (params.size() != 1)
	    {
	      cout << "ERROR: wrong number of parameters, correct input file" << endl;
	      return 0;
	    }
	  driver_st = read_state(params[0]);
	  prepare_state(driver_st);
	  // require a topology update
	  require_topo_update = true;
	  // require a coarse-graining update
	  require_CG_update = true;
	}

      // apply state transformations from a file
      else if (command == "transforms_file")
	{
	  if (params.size() != 1)
	    {
	      cout << "ERROR: wrong number of parameters, correct input file" << endl;
	      return 0;
	    }
	  driver_tr = read_transforms(params[0]);
	  apply_transforms(driver_tr);
	  // require a topology update
	  require_topo_update = true;
	  // require a coarse-graining update
	  require_CG_update = true;
	}

      // write the state to an output file
      else if (command == "output_state")
	{
	  if (params.size() != 1)
	    {
	      cout << "ERROR: wrong number of parameters, correct input file" << endl;
	      return 0;
	    }
	  write_state(params[0],dump_state());
	}

      // apply state transformations from a file
      else if (command == "regions_file")
	{
	  if (params.size() != 2)
	    {
	      cout << "ERROR: wrong number of parameters, correct input file" << endl;
	      return 0;
	    }
	  driver_rg.clear();
	  driver_rg = read_regions(params[0],stoi(params[1]));
	  // flag for regions being present
	  regions_present = true;
	}

      // apply state transformations from a file
      else if (command == "dump_regions")
	{
	  if (params.size() != 1)
	    {
	      cout << "ERROR: wrong number of parameters, correct input file" << endl;
	      return 0;
	    }
	  if (regions_present == true)
	    {
	      
	      // update topology before updating regions
	      if (require_topo_update == true)
		{
		  solve_topology();
		  require_topo_update = false;
		}
 
	      update_region_counts(driver_rg);
	      dump_regions(params[0],driver_rg);
	    }
	  else
	    {
	      cout << "  missing (read_regions), no regions to dump" << endl;
	    }
	}

      // write the topology to an output file
      else if (command == "dump_topology")
	{
	  if (params.size() != 2)
	    {
	      cout << "ERROR: wrong number of parameters, correct input file" << endl;
	      return 0;
	    }
	  // update topology before dumping
	  if (require_topo_update == true)
	    {
	      solve_topology();
	      require_topo_update = false;
	    }
	  dump_topology(params[0],stoi(params[1]));
	}

      // solve the topology of the current state
      else if (command == "solve_topology")
	{
	  solve_topology();
	  // disable flag for topology updating after an update
	  require_topo_update = false;
	}

      // update the CG map
      else if (command == "update_CG_map")
	{
	  if (params.size() != 1)
	    {
	      cout << "ERROR: wrong number of parameters, correct input file" << endl;
	      return 0;
	    }
	  // update topology before dumping
	  if (require_topo_update == true)
	    {
	      solve_topology();
	      require_topo_update = false;
	    }
	  driver_CG = update_CG_map(stoi(params[0]));
	}

      // update the CG map
      else if (command == "dump_CG_map")
	{
	  if (params.size() != 3)
	    {
	      cout << "ERROR: wrong number of parameters, correct input file" << endl;
	      return 0;
	    }
	  // update topology before updating CG_map
	  if (require_topo_update == true)
	    {
	      solve_topology();
	      require_topo_update = false;
	    }
	  // update CG_map before dumping
	  if (require_CG_update == true)
	    {
	      driver_CG = update_CG_map(stoi(params[1]));
	      require_CG_update = false;
	    }
	  dump_CG_map(params[0],stoi(params[2]),driver_CG);
	}

      // seed the PRNG
      else if (command == "prng_seed")
	{
	  if (params.size() != 1)
	    {
	      cout << "ERROR: wrong number of parameters, correct input file" << endl;
	      return 0;
	    }
	  // seed the PRNG
	  prng_seed(stoi(params[0]));
	  
	}

      // print the current state
      else if (command == "print")
	{
	  // update topology before printing
	  if (require_topo_update == true)
	    {
	      solve_topology();
	      require_topo_update = false;
	    }
	  print_tree();
	}

    }
  cout << "---END EXECUTING DIRECTIVES---\n" << endl;
  
  return 0;
}
