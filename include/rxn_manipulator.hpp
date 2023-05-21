#ifndef INCLUDE_RXN_MANIPULATOR_HPP
#define INCLUDE_RXN_MANIPULATOR_HPP

#include <vector>

struct species_count
{
  int id;
  int N;
};

struct reaction
{
  double k;
  std::vector<species_count> inputs, outputs;
};


class rxn_manipulator
{
public:

  // constructor and destructor
  rxn_manipulator();
  ~rxn_manipulator();

  void reset_reaction(reaction &r);
  void add_reaction_input(reaction &r, int id, int n);
  void add_reaction_output(reaction &r, int id, int n);
  void add_reaction_rate(reaction &r, double k);
  
private:


};

#endif
