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


// set the btree
void LAMMPS_sys::set_btree(btree_state in_state)
{
  internal_btree.prepare_state(in_state);
  internal_btree.solve_topology();
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


// concatenate ellipsoid arrays
void LAMMPS_sys::cat_ellipsoid_array(ellipsoid_array &in_ellipsoids, ellipsoid_array &cat_ellipsoids)
{

  // temporary array used during concatenation
  ellipsoid_array temp_ellipsoids;
  
  // get the array sizes
  int N_old = in_ellipsoids.get_N();
  int N_cat = cat_ellipsoids.get_N();
  int N_new = N_old + N_cat;

  // store the current contents of the array
  if (N_old > 0)
    {

      temp_ellipsoids.set_N(N_old);

      for (int i=0; i<N_old; i++)
	{
	  temp_ellipsoids.set_ellipsoid(i,in_ellipsoids.get_ellipsoid(i));
	}      
      
    }

  // resize the array
  in_ellipsoids.set_N(N_new);

  // copy the prior contents of the array
  if (N_old > 0)
    {
      for (int i=0; i<N_old; i++)
	{
	  in_ellipsoids.set_ellipsoid(i,temp_ellipsoids.get_ellipsoid(i));
	}
    }

  // concatenate the new contents
  if (N_cat > 0)
    {
      for (int i=0; i<N_cat; i++)
	{ 
	  in_ellipsoids.set_ellipsoid(i+N_old,cat_ellipsoids.get_ellipsoid(i));
	}
    }
  
}


// set the bonds in the system
void LAMMPS_sys::set_bonds()
{
  int **c, *t, N;

  c = nullptr;
  t = nullptr;

  internal_btree.prepare_bonds(c,t,N,1);

  // set the bond array size
  bonds.set_N(N);

  // fill the bond array
  bond temp_b;
  temp_b.id = 0;
  for (int i=0; i<N; i++)
    {
      temp_b.type = t[i];
      temp_b.i = c[i][0];
      temp_b.j = c[i][1];
      bonds.set_bond(i,temp_b);
    }
  
  // destroy the arrays
  for (int i=0; i<N; i++)
    {
      delete[] c[i];
    }
  delete[] c;
  delete[] t;
}


// set the angless in the system
void LAMMPS_sys::set_angles()
{
  int **c, *t, N;

  c = nullptr;
  t = nullptr;

  internal_btree.prepare_angles(c,t,N,1);

  // set the angle array size
  angles.set_N(N);

  // fill the angle array
  angle temp_a;
  temp_a.id = 0;
  for (int i=0; i<N; i++)
    {
      temp_a.type = t[i];
      temp_a.i = c[i][0];
      temp_a.j = c[i][1];
      temp_a.k = c[i][2];
      angles.set_angle(i,temp_a);
    }
  
  // destroy the arrays
  for (int i=0; i<N; i++)
    {
      delete[] c[i];
    }
  delete[] c;
  delete[] t;
}


// set the types
void LAMMPS_sys::set_mono_types(int base_type)
{
  int *t, N;

  t = nullptr;

  internal_btree.prepare_types(t,N,base_type);

  // set types in the atom array
  mono_atoms.set_types(t);
  
  delete[] t;
}


// prepare the topology
void LAMMPS_sys::prepare_topology()
{
  set_bonds();
  set_angles();
}


// function to create test data for testing
void LAMMPS_sys::prepare_test_data()
{

  // initialize the mono, ribo, and bdry atoms
  mono_atoms.set_N(internal_btree.total_size());
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

  cat_atom_array(atoms,mono_atoms);
  cat_atom_array(atoms,ribo_atoms);
  cat_atom_array(atoms,bdry_atoms);

  ellipsoids.set_N(0);

  mono_ellipsoids.set_min_id(1);
  cat_ellipsoid_array(ellipsoids,mono_ellipsoids);
  ribo_ellipsoids.set_min_id(mono_ellipsoids.get_N()+1);
  cat_ellipsoid_array(ellipsoids,ribo_ellipsoids);
  ellipsoids.normalize_quats();
  
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
  set_mono_types(3);
  
  mono_atoms.set_ellipsoid_flag_all(1);
  ribo_atoms.set_ellipsoid_flag_all(1);
  bdry_atoms.set_ellipsoid_flag_all(0);

  vec mono_shape;
  mono_shape.x = 34.0;
  mono_shape.y = 34.0;
  mono_shape.z = 34.0;

  vec ribo_shape;
  ribo_shape.x = 200.0;
  ribo_shape.y = 200.0;
  ribo_shape.z = 200.0;

  mono_ellipsoids.set_shape_all(mono_shape);
  ribo_ellipsoids.set_shape_all(ribo_shape);

  prepare_topology();
  
  merge_system_components();
  
  calc_bbox();
}


// write the system to a data file
void LAMMPS_sys::write_data(string data_filename)
{

  internal_btree.print_tree();

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

      data_file << "\n" << endl;
      
      // write atom information
      atoms.write(data_file);

      // write ellipsoid information
      ellipsoids.write(data_file);

      // write bond information
      bonds.write(data_file);

      // write angle information
      angles.write(data_file);
      
    }

  data_file.close();
  
}
