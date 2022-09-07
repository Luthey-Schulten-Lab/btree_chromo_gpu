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
  mono_atoms.set_N(internal_btree.total_size());
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

  mono_ellipsoids.set_N(mono_atoms.get_N());
  ribo_ellipsoids.set_N(ribo_atoms.get_N());

  mono_atoms.set_mol_id_all(3);

  bonds.set_N(10);
  angles.set_N(10);
  
  
}


// merge the system components
void LAMMPS_sys::merge_system_components()
{

  // reinitialize the total atom array
  atoms.set_N(0);

  // concatenate the mono, ribo, and bdry atom arrays
  cat_atom_array(atoms,mono_atoms);
  cat_atom_array(atoms,ribo_atoms);
  cat_atom_array(atoms,bdry_atoms);

  // reinitialize the total ellipsoid array
  ellipsoids.set_N(0);

  // concatenate the mono and ribo ellipsoid arrays
  mono_ellipsoids.set_min_id(1);
  cat_ellipsoid_array(ellipsoids,mono_ellipsoids);
  ribo_ellipsoids.set_min_id(mono_ellipsoids.get_N()+1);
  cat_ellipsoid_array(ellipsoids,ribo_ellipsoids);

  // normalize the quaternions in the total ellipsoid array
  ellipsoids.normalize_quats();
  
}


// create an appropriately sized bounding box
void LAMMPS_sys::calc_bbox()
{
  
  atom a;

  // determine hte maximum and minimum coordinates from the atom array
  for (int i=0; i<atoms.get_N(); i++)
    {
      a = atoms.get_atom(i);

      if (i == 0)
	{

	  // initialize maximums
	  bbox.r_max.x = a.r.x;
	  bbox.r_max.y = a.r.y;
	  bbox.r_max.z = a.r.z;

	  // initialize minimums
	  bbox.r_min.x = a.r.x;
	  bbox.r_min.y = a.r.y;
	  bbox.r_min.z = a.r.z;
	  
	}
      else
	{

	  // update maximums
	  if (a.r.x > bbox.r_max.x) bbox.r_max.x = a.r.x;
	  if (a.r.y > bbox.r_max.y) bbox.r_max.y = a.r.y;
	  if (a.r.z > bbox.r_max.z) bbox.r_max.z = a.r.z;

	  // update minimums
	  if (a.r.x < bbox.r_min.x) bbox.r_min.x = a.r.x;
	  if (a.r.y < bbox.r_min.y) bbox.r_min.y = a.r.y;
	  if (a.r.z < bbox.r_min.z) bbox.r_min.z = a.r.z;
      
	}
    }

  vec r_mid = vqm.v_linterp(0.5,bbox.r_min,bbox.r_max);
  vec dr = vqm.v_xpy(bbox.r_max,vqm.v_inv(bbox.r_min));
  dr = vqm.v_ax(0.5,dr);

  double s = 1.2; // scale factor

  // create bbox that is centered about the system and scaled to enclose all particles
  bbox.r_min = vqm.v_axpy(-s,dr,r_mid);
  bbox.r_max = vqm.v_axpy(s,dr,r_mid);
  
}


// finalize the system
void LAMMPS_sys::finalize_system()
{

  // prepare boundary atoms
  b_surf.generate_sphere(BD_l.r_sphere,BD_l.r_bdry);
  bdry_atoms.set_N(b_surf.get_N_verts());
  bdry_atoms.set_coords(b_surf.get_coords());
  // b_surf.write_xyz("/home/ben/Workspace/btree_chromo/test_case/test_b_surf.xyz");

  // place the ribo and bdry atoms in individual molecules for convenience
  bdry_atoms.set_mol_id_all(1);
  ribo_atoms.set_mol_id_all(2);

  // set the mono, ribo, and bdry types
  bdry_atoms.set_type_all(1);
  ribo_atoms.set_type_all(2);
  set_mono_types(3);

  // set the ellipsoid flags for mono, ribo, and bdry atoms
  mono_atoms.set_ellipsoid_flag_all(1);
  ribo_atoms.set_ellipsoid_flag_all(1);
  bdry_atoms.set_ellipsoid_flag_all(0);

  // set the ellipsoid shapes for mono and ribo atoms
  mono_ellipsoids.set_shape_all(BD_l.mono_shape);
  ribo_ellipsoids.set_shape_all(BD_l.ribo_shape);

  // prepare the system's topology based on the current binary tree
  prepare_topology();

  // merge the individual arrays before writing the system
  merge_system_components();

  // determine the bounding box
  calc_bbox();
  
}


// read the the Brownian dynamics lengths
void LAMMPS_sys::read_BD_lengths(string lengths_filename)
{
  fstream lengths_file;

  string param_delim, param, val;
  int delim;
  
  string line;

  param_delim = "=";

  lengths_file.open(lengths_filename, ios::in);

  if (!lengths_file.is_open())
    {
      cout << "ERROR: file not opened in read_BD_lengths" << endl;
    }
  else
    {
      while (1)
	{
	  lengths_file >> line;
	  if (lengths_file.eof()) break;
	  

	  if ((line.length() > 0) &&
	      (line.find("#") != 0))
	    {

	      delim = line.find(param_delim);

	      if (delim != -1)
		{

		  param = line.substr(0,delim);
		  val = line.substr(delim+1,line.length());

		  // cout << param << "=" << val << endl;

		  if (param == "r_sphere")
		    {
		      BD_l.r_sphere = stod(val);
		    }

		  else if (param == "r_bdry")
		    {
		      BD_l.r_bdry = stod(val);
		    }

		  else if (param == "mono_shape_x")
		    {
		      BD_l.mono_shape.x = stod(val);
		    }

		  else if (param == "mono_shape_y")
		    {
		      BD_l.mono_shape.y = stod(val);
		    }

		  else if (param == "mono_shape_z")
		    {
		      BD_l.mono_shape.z = stod(val);
		    }

		  else if (param == "ribo_shape_x")
		    {
		      BD_l.ribo_shape.x = stod(val);
		    }

		  else if (param == "ribo_shape_y")
		    {
		      BD_l.ribo_shape.y = stod(val);
		    }

		  else if (param == "ribo_shape_z")
		    {
		      BD_l.ribo_shape.z = stod(val);
		    }
		}
	    }
	} // end while loop
      lengths_file.close();
    }
}


// read the monomer coordinates - disallow resizing
int LAMMPS_sys::read_mono_coords(string coords_filename, string order)
{
  return mono_atoms.read_bin_coords(coords_filename,order,true);
}


// read the ribosome coordinates - allow resizing
int LAMMPS_sys::read_ribo_coords(string coords_filename, string order)
{
  return ribo_atoms.read_bin_coords(coords_filename,order,true);
}

// read the boundary coordinates - allow resizing
int LAMMPS_sys::read_bdry_coords(string coords_filename, string order)
{
  return bdry_atoms.read_bin_coords(coords_filename,order,true);
}


// write the system to a data file
void LAMMPS_sys::write_data(string data_filename)
{

  // internal_btree.print_tree();

  prepare_test_data();


  // finalize the system before printing
  finalize_system();

  
  // begin writing data file
  
  fstream data_file;

  data_file.open(data_filename, ios::out);

  if (!data_file.is_open())
    {
      cout << "ERROR: file not opened in write_data" << endl;
    }
  else
    {

      // write system summary
      data_file << "# LAMMPS data file for replicating chromosomes formed of rigid body monomers\n" << endl;

      data_file << atoms.get_N() << "\t\tatoms" << endl;
      data_file << N_atom_types << "\t\tatom types" << endl;
      data_file << ellipsoids.get_N() << "\t\tellipsoids" << endl;
      data_file << bonds.get_N() << "\t\tbonds" << endl;
      data_file << N_bond_types << "\t\tbond types" << endl;
      data_file << angles.get_N() << "\t\tangles" << endl;
      data_file << N_angle_types << "\t\tangle types" << endl;

      data_file << "\n" << endl;

      data_file << bbox.r_min.x << "\t" << bbox.r_max.x << "\txlo xhi" << endl;
      data_file << bbox.r_min.y << "\t" << bbox.r_max.y << "\tylo yhi" << endl;
      data_file << bbox.r_min.z << "\t" << bbox.r_max.z << "\tzlo zhi" << endl;
      
      // write atom information
      atoms.write(data_file);

      // write ellipsoid information
      ellipsoids.write(data_file);

      // write bond information
      bonds.write(data_file);

      // write angle information
      angles.write(data_file);


      data_file.close();
    }
  
}
