#ifndef INCLUDE_BTREE_DRIVER_HPP
#define INCLUDE_BTREE_DRIVER_HPP

#include <btree/btree.hpp>
#include <btree/gillespie_solver.hpp>

using namespace std;

struct drctv_reqs
{
  bool topo_update;
  bool CG_update;
  bool regions_present;
  bool rep_model_present;
};


class btree_driver //: protected btree
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

  ////////////////
  // directives // 
  ////////////////

  // input-ouput
  int input_state(vector<string> &params, drctv_reqs &reqs);
  int output_state(vector<string> &params);

  // transforms
  int transforms_file(vector<string> &params, drctv_reqs &reqs);
  int transform(vector<string> &params, drctv_reqs &reqs);
  int random_transforms(vector<string> &params, drctv_reqs &reqs);

  // regions
  int regions_file(vector<string> &params, drctv_reqs &reqs);
  int dump_regions(vector<string> &params, drctv_reqs &reqs);

  // topology
  int dump_topology(vector<string> &params, drctv_reqs &reqs);
  int update_topology(drctv_reqs &reqs);

  // coarse-graining
  int update_CG_map(vector<string> &params, drctv_reqs &reqs);
  int dump_CG_map(vector<string> &params, drctv_reqs &reqs);

  // miscellaneous
  int prng_seed(vector<string> &params);
  int print_state(drctv_reqs &reqs);

  // replication model
  int replication_model(vector<string> &params, drctv_reqs &reqs);

  // classes
  gillespie_solver solver;
  btree driver_bt;

  // structs
  rep_model_params driver_rep_model; // replication model
  btree_state driver_st; // state structure
  btree_transforms driver_tr; // transform structure
  vector<chromo_region> driver_rg; // vector of chromo_regions
  vector<string> drctvs; // set of directives
  CG_map driver_CG; // coarse-graining map

};

#endif
