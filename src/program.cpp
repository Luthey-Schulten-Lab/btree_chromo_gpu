#include <btree/btree_driver.hpp>
#include <program.hpp>

int main(int argc, char *argv[])
{

  int error_code = 0;
  btree_driver driver;

  if (argc == 1)
    {
      cout << "ERROR: missing directive file as command-line argument" << endl;
      return 0;
    }

  string drctv_filename(argv[argc-1]);

  // read the directives from the directive file
  driver.read_directives(drctv_filename);

  // print the directives to be executed
  driver.print_directives();

  // parse the directives into commands and parameters
  driver.parse_directives();
  
  // expand the metacommands


  // validate the numbers of parameters
  error_code = driver.validate_command_sequence_parameters();
  if (error_code != 0)
    {
      cout << "error during parameter validation" << endl;
      return 1;
    }
  else
    {
      cout << "\n--- VALID COMMAND PARAMETERS ---\n" << endl;
    }

  // validate the command sequence
  error_code = driver.validate_command_sequence();
  if (error_code != 0)
    {
      cout << "error during sequence validation" << endl;
      return 1;
    }
  else
    {
      cout << "\n--- VALID COMMAND SEQUENCE ---\n" << endl;
    }

  // print the set of commands
  driver.print_commands();
    
  // execute the commands
  error_code = driver.execute_commands();
  if (error_code != 0)
    {
      cout << "error during command execution" << endl;
    } 

  return 0;
}
