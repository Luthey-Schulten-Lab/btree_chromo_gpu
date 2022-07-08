#include <btree/btree.hpp>
#include <btree/btree_driver.hpp>
#include <program.hpp>

int main()
{

  int error_code = 0;
  btree_driver driver;

  driver.read_directives("/home/ben/Workspace/btree_chromo/test_case/directives.inp");
  driver.print_directives();

  error_code = driver.execute_directives();

  if (error_code != 0)
    {
      cout << "error during directive execution" << endl;
    }
  
  // driver.~btree_driver();
  
  // btree btree_test;
  // fork_rho fork_rho_test;
  // btree_state btree_state_test;
  // btree_transforms btree_transforms_test;


  // btree_state_test = btree_test.read_state("/home/ben/Workspace/btree_chromo/test_case/in_state_test.txt");
  // btree_test.prepare_state(btree_state_test);
  // btree_test.print_tree();

  // btree_transforms_test = btree_test.read_transforms("/home/ben/Workspace/btree_chromo/test_case/state_trans_test.txt");

  // btree_test.apply_transforms(btree_transforms_test);
  // btree_test.print_tree();

  // btree_state_test = btree_test.dump_state();

  // btree_test.write_state("/home/ben/Workspace/btree_chromo/test_case/out_state_test.txt"
  // 			 ,btree_state_test);
  
  // btree_test.destroy_tree();

  // btree_test.print_tree();

  // btree_test.initialize_tree(10);

  // btree_test.print_tree();
  
  program::foo();

  return 0;
}
