#include <btree/btree.hpp>
#include <program.hpp>

int main()
{

  btree btree_test;
  fork_rho fork_rho_test;
  btree_state btree_state_test;

  btree_state_test.size = 100;
  fork_rho_test.fork = "m";
  fork_rho_test.rho = 100;
  btree_state_test.fork_rhos.push_back(fork_rho_test);
  fork_rho_test.fork = "ml";
  fork_rho_test.rho = 100;
  btree_state_test.fork_rhos.push_back(fork_rho_test);
  fork_rho_test.fork = "mlr";
  fork_rho_test.rho = 25;
  btree_state_test.fork_rhos.push_back(fork_rho_test);
  fork_rho_test.fork = "mr";
  fork_rho_test.rho = 15;
  btree_state_test.fork_rhos.push_back(fork_rho_test);

  btree_test.prepare_state(btree_state_test);

  // btree_test.initialize_tree(100);

  // btree_test.foo();

  // btree_test.print_tree();

  // btree_test.branch("ml");
  // btree_test.branch("m");
  // btree_test.grow_at_branch("m",10);
  // btree_test.branch("ml");
  // btree_test.grow_at_branch("ml",15);

  // btree_test.print_tree();
  // btree_test.grow_at_branch("m",150);
  // btree_test.branch("mlr");
  // btree_test.grow_at_branch("ml",150);
  // btree_test.branch("mr");

  btree_test.print_tree();

  btree_test.destroy_tree();

  btree_test.print_tree();

  btree_test.initialize_tree(10);

  btree_test.print_tree();
  
  program::foo();

  return 0;
}
