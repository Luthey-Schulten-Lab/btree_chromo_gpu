#ifndef INCLUDE_BTREE_DRIVER_HPP
#define INCLUDE_BTREE_DRIVER_HPP

#include <btree/btree.hpp>

using namespace std;


class btree_driver : protected btree
{
public:

  // constructor and destructor
  btree_driver();
  ~btree_driver();

  // read_directives
  void read_directives(string drctvs_filename);

  // print_directives
  void print_directives();

  // execute directives
  int execute_directives();

private:

  btree_state driver_st; // state structure
  btree_transforms driver_tr; // transform structure
  vector<string> drctvs; // set of directives

};

#endif
