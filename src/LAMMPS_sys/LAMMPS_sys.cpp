#include <LAMMPS_sys/LAMMPS_sys.hpp>

// constructor
LAMMPS_sys::LAMMPS_sys()
{
  N_atom_types = 0;
  N_angle_types = 0;
  N_bond_types = 0;
}


// destructor
LAMMPS_sys::~LAMMPS_sys()
{
}


// concatenate atom arrays
void LAMMPS_sys::cat_atom_array(atom_array &in_atoms, atom_array &cat_atoms)
{

  // temporary array used during concatenation
  atom_array temp_atoms;
  
  // get the array sizes
  int N_old = in_atoms.get_N();
  int N_cat = cat_atoms.get_N();
  int N_new = N_old + N_cat;

  // store the current contents of the array
  if (N_old > 0)
    {

      temp_atoms.set_N(N_old);

      for (int i=0; i<N_old; i++)
	{
	  temp_atoms.set_atom(i,in_atoms.get_atom(i));
	}      
      
    }

  // resize the array
  in_atoms.set_N(N_new);

  // copy the prior contents of the array
  if (N_old > 0)
    {
      for (int i=0; i<N_old; i++)
	{
	  in_atoms.set_atom(i,temp_atoms.get_atom(i));
	}
    }

  // concatenate the new contents
  if (N_cat > 0)
    {
      for (int i=0; i<N_cat; i++)
	{ 
	  in_atoms.set_atom(i+N_old,cat_atoms.get_atom(i));
	}
    }

  in_atoms.reset_ids();
  
}


// function to create test data for testing
void LAMMPS_sys::prepare_test_data()
{

  // initialize the mono, ribo, and bdry atoms
  mono_atoms.set_N(10);
  ribo_atoms.set_N(5);
  bdry_atoms.set_N(20);

  mono_ellipsoids.set_N(mono_atoms.get_N());
  ribo_ellipsoids.set_N(ribo_atoms.get_N());

  mono_atoms.set_mol_id_all(3);

  bonds.set_N(10);
  angles.set_N(10);
  
  
}


// merge the system components
void LAMMPS_sys::merge_system_components()
{

  atoms.set_N(0);
  ellipsoids.set_N(0);

  cat_atom_array(atoms,mono_atoms);
  cat_atom_array(atoms,ribo_atoms);
  cat_atom_array(atoms,bdry_atoms);
  
}


// create an appropriately sized bounding box
void LAMMPS_sys::calc_bbox()
{
  
}


// finalize the system
void LAMMPS_sys::finalize_system()
{

  bdry_atoms.set_mol_id_all(1);
  ribo_atoms.set_mol_id_all(2);
  
  bdry_atoms.set_type_all(1);
  ribo_atoms.set_type_all(2);
  
  mono_atoms.set_ellipsoid_flag_all(1);
  ribo_atoms.set_ellipsoid_flag_all(1);
  bdry_atoms.set_ellipsoid_flag_all(0);
  
  merge_system_components();
  calc_bbox();
}


// write the system to a data file
void LAMMPS_sys::write_data(string data_filename)
{

  prepare_test_data();


  // finalize the system before printing
  finalize_system();

  
  // begin writing data file
  
  fstream data_file;

  data_file.open(data_filename, ios::out);

  if (!data_file)
    {
      cout << "ERROR: file not opened in write_data" << endl;
    }
  else
    {

      // write system summary
      data_file << "# LAMMPS data file for replicating chromosomes formed of rigid body monomers\n" << endl;

      data_file << atoms.get_N() << "\t\tatoms" << endl;
      data_file << N_atom_types << "\t\tatom types" << endl;
      data_file << ellipsoids.get_N() << "\t\tatoms" << endl;
      data_file << bonds.get_N() << "\t\tbonds" << endl;
      data_file << N_bond_types << "\t\tbond types" << endl;
      data_file << angles.get_N() << "\t\tangles" << endl;
      data_file << N_angle_types << "\t\tangle types" << endl;

      data_file << "\n\n" << endl;
      
      // write atom information
      atoms.write(data_file);

      // write ellipsoid information

      // write bond information

      // write angle information
      
    }

  data_file.close();
  
}
