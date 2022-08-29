#include <LAMMPS_sys/atom_array.hpp>

// constructor
atom_array::atom_array()
{
  N = -1;
  atoms = nullptr;
}

// destructor
atom_array::~atom_array()
{
  destroy_atoms();
}

// set the number of elements
void atom_array::set_N(int N)
{
  this->N = N;
  initialize_atoms();
}

// get the number of elements
int atom_array::get_N()
{
  return N;
}

// initialize the array
void atom_array::initialize_atoms()
{
  destroy_atoms();

  if (N > 0)
    {

      atoms = new atom[N];

      for (int i=0; i<N; i++)
	{
	  // atom id
	  atoms[i].id = i + 1;
	  // type
	  atoms[i].type = -1;
	  // molecule id
	  atoms[i].mol_id = -1;
	  // ellipsoid flag
	  atoms[i].ellipsoid_flag = -1;
	  // density
	  atoms[i].density = 0.0;
	  // position
	  atoms[i].r.x = 0.0;
	  atoms[i].r.y = 0.0;
	  atoms[i].r.z = 0.0;
	}
      
    }
  
}

// destroy the array
void atom_array::destroy_atoms()
{
  if (atoms != nullptr)
    {
      delete[] atoms;
      atoms = nullptr;
    }
}


// setter for elements
void atom_array::set_atom(int i, atom a)
{
  atoms[i] = a;
}


// getter for elements
atom atom_array::get_atom(int i)
{
  return atoms[i];
}


// reset the ids
void atom_array::reset_ids()
{
  for (int i=0; i<N; i++)
    {
      atoms[i].id = i + 1;
    }
}


// set all types
void atom_array::set_type_all(int type)
{
  for (int i=0; i<N; i++)
    {
      atoms[i].type = type;
    }
}


// set all ellipsoid_flags
void atom_array::set_ellipsoid_flag_all(int ellipsoid_flag)
{
  for (int i=0; i<N; i++)
    {
      atoms[i].ellipsoid_flag = ellipsoid_flag;
    }
}


// set all mol_ids
void atom_array::set_mol_id_all(int mol_id)
{
  for (int i=0; i<N; i++)
    {
      atoms[i].mol_id = mol_id;
    }
}


// set all densities
void atom_array::set_density_all(double density)
{
  for (int i=0; i<N; i++)
    {
      atoms[i].density = density;
    }
}


// write to stream
void atom_array::write(fstream &data_file)
{
  data_file << "\nAtoms # atom-ID atom-type x y z molecule-ID ellipsoid-flag density\n" << endl;

  for (int i=0; i<N; i++)
    {
      data_file << atoms[i].id << "\t"
		<< atoms[i].type << "\t"
		<< atoms[i].r.x << "\t"
		<< atoms[i].r.y << "\t"
		<< atoms[i].r.z << "\t"
		<< atoms[i].mol_id << "\t"
		<< atoms[i].ellipsoid_flag << "\t"
		<< atoms[i].density << endl;
    }
  
}
