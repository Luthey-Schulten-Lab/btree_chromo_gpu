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


// set types
void atom_array::set_types(int *&types)
{
  for (int i=0; i<N; i++)
    {
      atoms[i].type = types[i];
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


// set individual coordinate
void atom_array::set_coord(int i, vec r)
{
     atoms[i].r.x = r.x;
     atoms[i].r.y = r.y;
     atoms[i].r.z = r.z;
}


// set element-wise coordinates
void atom_array::set_coords(vector<vec> rs)
{
  for (int i=0; i<N; i++)
    {
      set_coord(i,rs[i]);
    }
}



// set the particle coordinates from a 1D array of doubles
void atom_array::set_coords_arr(double *&x, string order)
{

  vec r;
  
  if (order == "col")
    {
      for (int i=0; i<N; i++)
	{
	  r.x = x[i];
	  r.y = x[N+i];
	  r.z = x[2*N+i];
	  set_coord(i,r);
	}
    }
  else if (order == "row")
    {
      for (int i=0; i<N; i++)
	{
	  r.x = x[3*i];
	  r.y = x[3*i+1];
	  r.z = x[3*i+2];
	  set_coord(i,r);
	}
    }
}


// read coordinates from binary file
int atom_array::read_bin_coords(string data_filename, string order, bool force_resize)
{
  fstream data_file;

  data_file.open(data_filename, ios::in | ios::binary | ios::ate);

  // int data_size = 4;
    
  if (!data_file.is_open())
    {
      cout << "ERROR: file not opened in write_data" << endl;
      return 1;
    }
  else
    {

      // read in the binary file
      streampos size = data_file.tellg();
      char *memblock;
      double *x;

      int data_size = sizeof(double);
      char vals[sizeof(double)];
      int N_data = size/data_size;
      cout << "N_data = " << N_data << endl;
      int N_bin_atoms = N_data/3;
      cout << "N_bin_atoms = " << N_bin_atoms << endl;

      // test for resizing
      if (force_resize)
	{
	  set_N(N_bin_atoms);
	}
      else
	{
	  if (N_bin_atoms != get_N())
	    {
	      cout << "ERROR: incorrect atom_array size" << endl;
	      return 1;
	    }
	}

      memblock = new char[size];
      
      data_file.seekg(0, ios::beg);

      data_file.read(memblock,size);

      data_file.close();

      x = new double[N_data];

      for (int i=0; i<N_data; i++)
	{
	  for (int j=0; j<data_size; j++)
	    {
	      vals[j] = memblock[i*data_size+j];
	    }
	  memcpy(&x[i], vals, 8);
	}

      delete memblock;
      
      set_coords_arr(x,order);

      delete x;
      return 0;
      
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
		<< atoms[i].r.x << "  \t"
		<< atoms[i].r.y << "  \t"
		<< atoms[i].r.z << "  \t"
		<< atoms[i].mol_id << "\t"
		<< atoms[i].ellipsoid_flag << "\t"
		<< atoms[i].density << endl;
    }
  
}


// write the boundary coordinates to an xyz file
void atom_array::write_xyz(string data_filename)
{

  // begin writing data file
  
  fstream data_file;

  data_file.open(data_filename, ios::out);

  if (!data_file)
    {
      cout << "ERROR: file not opened in write_xyz" << endl;
    }
  else
    {

      // write system summary
      data_file << N << "\n" << endl;

      for (int i=0; i<N; i++)
	{
	  data_file << "C\t"
		    << atoms[i].r.x << "\t"
		    << atoms[i].r.y << "\t"
		    << atoms[i].r.z << endl;
	}
      
    }

  data_file.close();
  
}
