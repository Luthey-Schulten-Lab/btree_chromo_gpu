#ifndef INCLUDE_BINDING_REGION_HPP
#define INCLUDE_BINDING_REGION_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <memory>

#include <btree/btree.hpp>
#include <LAMMPS_sys/vec_quat_manipulator.hpp>

using namespace std;

class binding_region
{
public:

  binding_region(string leaf, int ll, int ul, int size,
		 bool completed,
		 bool ter_crossing, int mid_ll, int mid_ul);
  ~binding_region();

  // getters
  int get_size();
  string get_leaf();

  // randomly select an initial anchor within the region

  // test if queried index is within the region
  bool within_region(int q);

private:

  string leaf; // leaf that binding region belongs to
  int ll, ul, size; // lower limit and upper limit of indices
  bool completed;
  bool ter_crossing;
  int mid_ll, mid_ul;
  
};

#endif
