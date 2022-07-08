#include <btree/btree.hpp>
#include <program.hpp>

int main()
{

  btree btree_test;
  fork_rho fork_rho_test;
  btree_state btree_state_test;
  btree_transforms btree_transforms_test;

  btree_state_test.size = 1000;
  fork_rho_test.fork = "m";
  fork_rho_test.rho = 100;
  btree_state_test.fork_rhos.push_back(fork_rho_test);
  fork_rho_test.fork = "ml";
  fork_rho_test.rho = 150;
  btree_state_test.fork_rhos.push_back(fork_rho_test);
  fork_rho_test.fork = "mlr";
  fork_rho_test.rho = 25;
  btree_state_test.fork_rhos.push_back(fork_rho_test);
  fork_rho_test.fork = "mr";
  fork_rho_test.rho = 15;
  btree_state_test.fork_rhos.push_back(fork_rho_test);

  btree_test.prepare_state(btree_state_test);
  btree_test.print_tree();

  fork_rho_test.fork = "mlr";
  fork_rho_test.rho = 25;
  btree_transforms_test.fork_rhos.push_back(fork_rho_test);
  fork_rho_test.fork = "mrrr";
  fork_rho_test.rho = 15;
  btree_transforms_test.fork_rhos.push_back(fork_rho_test);

  btree_test.apply_transforms(btree_transforms_test);
  btree_test.print_tree();

  btree_state_test = btree_test.dump_state();

  btree_test.destroy_tree();

  cout << "remaking btree" << endl;
  btree_test.prepare_state(btree_state_test);
  btree_test.print_tree();

  btree_test.destroy_tree();

  btree_test.print_tree();

  btree_test.initialize_tree(10);

  btree_test.print_tree();
  
  program::foo();

  return 0;
}
