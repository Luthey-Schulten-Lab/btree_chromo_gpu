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

  // parse the directives
  
  
  // execute the directives
  error_code = driver.execute_directives();

  if (error_code != 0)
    {
      cout << "error during directive execution" << endl;
    } 

  return 0;
}
