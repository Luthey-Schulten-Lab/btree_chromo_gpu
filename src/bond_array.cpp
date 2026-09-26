#include <wcm_fmt.hpp>
#include <bond_array.hpp>

// constructor
bond_array::bond_array()
{
  N = -1;
  bonds = nullptr;
}

// destructor
bond_array::~bond_array()
{
  destroy_bonds();
}

// set the number of elements
void bond_array::set_N(int N)
{
  this->N = N;
  initialize_bonds();
}

// get the number of elements
int bond_array::get_N()
{
  return N;
}

// initialize the array
void bond_array::initialize_bonds()
{
  destroy_bonds();

  if (N > 0)
    {

      bonds = new bond[N];

      for (int i=0; i<N; i++)
	{
	  // bond id
	  bonds[i].id = i + 1;
	  // type
	  bonds[i].type = -1;
	  // bonded atom ids
	  bonds[i].i = -1;
	  bonds[i].j = -1;
	}
      
    }
  
}

// destroy the array
void bond_array::destroy_bonds()
{
  if (bonds != nullptr)
    {
      delete[] bonds;
      bonds = nullptr;
    }
}


// set bond element
void bond_array::set_bond(int i, bond b)
{
  bonds[i] = b;
  bonds[i].id = i + 1;
}


// get bond element
bond bond_array::get_bond(int i)
{
  return bonds[i];
}


// write to stream
// the same text as write()
void bond_array::append(std::string &out) const
{
  if (N > 0)
    {
      out += "\nBonds # bond-ID bond-type i j\n\n";
      for (int i=0; i<N; i++)
        {
      wcm_put_int(out, bonds[i].id); out += '\t';
      wcm_put_int(out, bonds[i].type); out += '\t';
      wcm_put_int(out, bonds[i].i); out += '\t';
      wcm_put_int(out, bonds[i].j); out += '\n';
        }
    }
}

void bond_array::write(std::fstream &data_file)
{

  if (N > 0)
    {
      data_file << "\nBonds # bond-ID bond-type i j\n" << std::endl;

      for (int i=0; i<N; i++)
	{
	  data_file << bonds[i].id << "\t"
		    << bonds[i].type << "\t"
		    << bonds[i].i << "\t"
		    << bonds[i].j << '\n';   // '\n', not std::endl (a flush per line)
	}
    }
  
}
