#include <rep_kinetics/rxn_manipulator.hpp>

// constructor
rxn_manipulator::rxn_manipulator()
{
}


// destructor
rxn_manipulator::~rxn_manipulator()
{
}

// reset the reaction
void rxn_manipulator::reset_reaction(reaction &r)
{
  r.inputs.clear();
  r.outputs.clear();
}

// add an input to a reaction
void rxn_manipulator::add_reaction_input(reaction &r, int id, int n)
{
  species_count s_c;
  s_c.id = id;
  s_c.N = n;
  r.inputs.push_back(s_c);
}

// add an output to a reaction
void rxn_manipulator::add_reaction_output(reaction &r, int id, int n)
{
  species_count s_c;
  s_c.id = id;
  s_c.N = n;
  r.outputs.push_back(s_c);
}

