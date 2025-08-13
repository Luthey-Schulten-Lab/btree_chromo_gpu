#include <rxn_manipulator.hpp>

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
  r.k = 0.0;
}

// add an input to a reaction
void rxn_manipulator::add_reaction_input(reaction &r, int id)
{
  r.inputs.push_back(id);
}

// add an output to a reaction
void rxn_manipulator::add_reaction_output(reaction &r, int id)
{
  r.outputs.push_back(id);
}

// add an output to a reaction
void rxn_manipulator::add_reaction_rate(reaction &r, double k)
{
  r.k = k;
}
