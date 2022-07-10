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
	      (line.find("#") != 1))
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

  string command, params;
  string cmd_delim = ":";
  int delim;

  for (string drctv: drctvs)
    {

      delim = drctv.find(cmd_delim);

      if (delim != -1)
	{
	  command = drctv.substr(0,delim);
	  params = drctv.substr(delim+1,drctv.length());
	}
      else
	{
	  command = drctv;
	}

      if (command == "input_state")
	{
	  driver_st = read_state(params);
	  prepare_state(driver_st);
	}
      else if (command == "transforms_file")
	{
	  driver_tr = read_transforms(params);
	  apply_transforms(driver_tr);
	}
      else if (command == "output_state")
	{
	  write_state(params,dump_state());
	}
      else if (command == "dump_topology")
	{
	  print_tree();
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
  return 0;
}
