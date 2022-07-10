#include <btree/btree.hpp>
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

  // driver.read_directives("/home/ben/Workspace/btree_chromo/test_case/directives.inp");
  driver.read_directives(drctv_filename);
  driver.print_directives();

  error_code = driver.execute_directives();

  if (error_code != 0)
    {
      cout << "error during directive execution" << endl;
    } 

  
  // program::foo();

  return 0;
}
