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

  cout << "Reading Directives:" << endl;
  cout << drctvs_filename << endl;

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

  cout << "\n---BEGIN EXECUTING DIRECTIVES---\n" << endl;

  for (string drctv: drctvs)
    {

      delim = drctv.find(cmd_delim);

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
      else
	{
	  command = drctv;

	  cout << "COMMAND: " << command << endl;
	  cout << "\n" << endl;
	}

      if (command == "input_state")
	{
	  driver_st = read_state(params[0]);
	  prepare_state(driver_st);
	}
      else if (command == "transforms_file")
	{
	  driver_tr = read_transforms(params[0]);
	  apply_transforms(driver_tr);
	}
      else if (command == "output_state")
	{
	  write_state(params[0],dump_state());
	}
      else if (command == "dump_topology")
	{
	  dump_topology(params[0],stoi(params[1]));
	}
      else if (command == "solve_topology")
	{
	  solve_topology();
	}
      else if (command == "print")
	{
	  print_tree();
	}

    }
  cout << "\n---END EXECUTING DIRECTIVES---\n" << endl;
  
  return 0;
}
