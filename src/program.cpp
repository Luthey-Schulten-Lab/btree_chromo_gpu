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

  // validate the numbers of parameters
  cout << "\n--- BEGIN COMMAND PARAMETER VALIDATION ---\n" << endl;
  error_code = driver.validate_command_sequence_parameters();
  if (error_code != 0)
    {
      cout << "\terror in command parameters" << endl;
      cout << "\n--- END COMMAND PARAMETER VALIDATION ---\n" << endl;
      return 1;
    }
  else
    {
      cout << "\tvalid command sequence" << endl;
      cout << "\n--- END COMMAND PARAMETER VALIDATION ---\n" << endl;
    }
  
  // expand the metacommands
  cout << "\n--- BEGIN METACOMMAND EXPANSION ---\n" << endl;
  error_code = driver.expand_metacommands();
  if (error_code != 0)
    {
      cout << "\terror during metacommand expansion" << endl;
      cout << "\n--- END METACOMMAND EXPANSION ---\n" << endl;
      return 1;
    }
  else
    {
      cout << "\tsuccessful metacommand expansion" << endl;
      cout << "\n--- END METACOMMAND EXPANSION ---\n" << endl;
    }

  // validate the command sequence
  cout << "\n--- BEGIN COMMAND SEQUENCE VALIDATION ---\n" << endl;
  error_code = driver.validate_command_sequence();
  if (error_code != 0)
    {
      cout << "\terror in command sequence" << endl;
      cout << "\n--- END COMMAND SEQUENCE VALIDATION ---\n" << endl;
      return 1;
    }
  else
    {
      cout << "\tvalid command sequence" << endl;
      cout << "\n--- END COMMAND SEQUENCE VALIDATION ---\n" << endl;
    }

  // print the set of commands
  driver.print_commands();
    
  // execute the commands
  // error_code = driver.execute_commands();
  // if (error_code != 0)
  //   {
  //     cout << "error during command execution" << endl;
  //   } 

  return 0;
  
}
