#include <LAMMPS_sys.hpp>

// constructor
LAMMPS_sys::LAMMPS_sys()
{
  N_atom_types = 9; // 0 bdry, 1 ribo, 2 mono(m), 3 ori(o), 4 ter(t), 5 fork(f), 6 anchor(a), 7 hinge(h), 8 extruded mono(e)
  N_angle_types = 4; // linear m-m/o/t-m, twist m-m/o/t-m, linear m-f-m, twist m-f-m
  N_bond_types = 3; // m/o/t-m/o/t, loops(a-h)
  
  loop_topo.prng_seed(0);

  t_s.bond = true;
  t_s.bending_angle = true;
  t_s.twisting_angle = true;
  r_s.ellipsoid = true;

  voxel_grid = VoxelGrid();
}


// destructor
LAMMPS_sys::~LAMMPS_sys()
{
}


void LAMMPS_sys::prng_seed(int s)
{
  loop_topo.prng_seed(s);
}


// set the btree
void LAMMPS_sys::set_btree(btree_state in_state)
{
  internal_btree.prepare_state(in_state);
  internal_btree.solve_topology();
  if (mono_atoms.get_N() != internal_btree.total_size())
    {
      mono_atoms.set_N(internal_btree.total_size());
    }
  if (mono_ellipsoids.get_N() != internal_btree.total_size())
    {
      mono_ellipsoids.set_N(internal_btree.total_size());
    }
}


// get the fork partitions
std::vector<fork_partition> LAMMPS_sys::get_all_fork_partitions()
{
  return internal_btree.get_all_fork_partitions();
}


// concatenate atom arrays
void LAMMPS_sys::cat_atom_array(atom_array &in_atoms, atom_array &cat_atoms)
{

  // temporary array used during concatenation
  atom_array temp_atoms;
  
  // get the array sizes
  int N_old = in_atoms.get_N();
  int N_cat = cat_atoms.get_N();
  int N_new = 0;
  if (N_old > 0) N_new += N_old;
  if (N_cat > 0) N_new += N_cat;

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
  int N_new = 0;
  if (N_old > 0) N_new += N_old;
  if (N_cat > 0) N_new += N_cat;

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


// switch the bonds
void LAMMPS_sys::switch_bonds(bool s)
{
  t_s.bond = s;
}


// switch the bonds
void LAMMPS_sys::switch_bending_angles(bool s)
{
  t_s.bending_angle = s;
}


// switch the twisting angles
void LAMMPS_sys::switch_twisting_angles(bool s)
{
  t_s.twisting_angle = s;
}


// switch whether to use ellipsoids
void LAMMPS_sys::switch_ellipsoids(bool s)
{
  r_s.ellipsoid = s;
}


// filter the bonds
void LAMMPS_sys::filter_bonds(int t)
{
  bond_array temp_bonds;
  bond t_b;

  int c;

  // count the number of remaining bonds
  c = 0;
  for (int i=0; i<bonds.get_N(); i++)
    {
      t_b = bonds.get_bond(i);
      if (t_b.type != t)
	{
	  c += 1;
	}
    }

  // size the temporary array
  temp_bonds.set_N(c);

  // fill the temporary array with the remaining bonds
  c = 0;
  for (int i=0; i<bonds.get_N(); i++)
    {
      t_b = bonds.get_bond(i);
      if (t_b.type != t)
	{
	  temp_bonds.set_bond(c,t_b);
	  c += 1;
	}
    }

  // resize the bond array
  bonds.set_N(temp_bonds.get_N());

  // copy the contents of the temporary array to the bond array
  for (int i=0; i<temp_bonds.get_N(); i++)
    {
      bonds.set_bond(i,temp_bonds.get_bond(i));
    }
}


// filter the angless
void LAMMPS_sys::filter_angles(int t)
{
  angle_array temp_angles;
  angle t_a;

  int c;

  // count the number of remaining angles
  c = 0;
  for (int i=0; i<angles.get_N(); i++)
    {
      t_a = angles.get_angle(i);
      if (t_a.type != t)
	{
	  c += 1;
	}
    }

  // size the temporary array
  temp_angles.set_N(c);

  // fill the temporary array with the remaining angles
  c = 0;
  for (int i=0; i<angles.get_N(); i++)
    {
      t_a = angles.get_angle(i);
      if (t_a.type != t)
	{
	  temp_angles.set_angle(c,t_a);
	  c += 1;
	}
    }

  // resize the angle array
  angles.set_N(temp_angles.get_N());

  // copy the contents of the temporary array to the angle array
  for (int i=0; i<temp_angles.get_N(); i++)
    {
      angles.set_angle(i,temp_angles.get_angle(i));
    }
}


// set the types
void LAMMPS_sys::set_mono_types(int base_type)
{
  int *t, N;

  N = get_N_mono();

  t = nullptr;

  internal_btree.prepare_types(t,N,base_type);

  // set types in the atom array
  mono_atoms.set_types(t);
  
  delete[] t;
}


// get the mono types
void LAMMPS_sys::get_types(int *&t)
{
  for (int i=0; i<get_N_total(); i++)
    {
      t[i] = atoms.get_atom(i).type;
    }
}


// prepare the topology
void LAMMPS_sys::prepare_topology()
{
  // set the bond topology
  set_bonds();

  // filter the bonds
  if (t_s.bond == false)
    {
      filter_bonds(1);
    }

  // set the angle topology
  set_angles();

  // filter the bending angles
  if (t_s.bending_angle == false)
    {
      filter_angles(1);
      filter_angles(3);
    }

  // filter the twisting angles
  if (!(t_s.twisting_angle && r_s.ellipsoid))
    {

      std::cout << "Filtering the twisting angles... " << std::endl;
      filter_angles(2);
      filter_angles(4);
    }
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


// generate a spherical boundary
void LAMMPS_sys::generate_spherical_bdry(double r, double x0, double y0, double z0)
{
  // generate a spherical boundary from an interpolated set of triangulated mesh
  b_surf.generate_sphere(r,BD_l.r_bdry);

  // size the bdry_atoms
  bdry_atoms.set_N(b_surf.get_N_verts());

  std::vector<vec> bdry_coords = b_surf.get_coords();

  vec r0 = vqm.v_new(x0,y0,z0);

  // translate the boundary
  for (size_t i=0; i<bdry_coords.size(); i++)
    {
      bdry_coords[i] = vqm.v_xpy(bdry_coords[i],r0);
    }

  // set the bdry_atoms to the coordinates
  bdry_atoms.set_coords(bdry_coords);
}

void LAMMPS_sys::generate_cylindrical_bdry(double l, double r, double x0, double y0, double z0)
{
    // generate a spherical boundary from an interpolated set of triangulated mesh
    b_surf.generate_cylinder(l, r,BD_l.r_bdry);

    // size the bdry_atoms
    bdry_atoms.set_N(b_surf.get_N_verts());

    std::vector<vec> bdry_coords = b_surf.get_coords();

    vec r0 = vqm.v_new(x0,y0,z0);

    // translate the boundary
    for (size_t i=0; i<bdry_coords.size(); i++)
    {
        bdry_coords[i] = vqm.v_xpy(bdry_coords[i],r0);
    }

    // set the bdry_atoms to the coordinates
    bdry_atoms.set_coords(bdry_coords);
}

void LAMMPS_sys::generate_spherocylindrical_bdry(double l, double r, double x0, double y0, double z0)
{
    // generate a spherocylindrical boundary from an interpolated set of triangulated mesh
    b_surf.generate_spherocylinder(l, r,BD_l.r_bdry);

    // size the bdry_atoms
    bdry_atoms.set_N(b_surf.get_N_verts());

    std::vector<vec> bdry_coords = b_surf.get_coords();

    vec r0 = vqm.v_new(x0,y0,z0);

    // translate the boundary
    for (size_t i=0; i<bdry_coords.size(); i++)
    {
        bdry_coords[i] = vqm.v_xpy(bdry_coords[i],r0);
    }

    // set the bdry_atoms to the coordinates
    bdry_atoms.set_coords(bdry_coords);
}

void LAMMPS_sys::generate_overlapping_spheres_bdry(double h, double r, double x0, double y0, double z0, double u, double v, double w)
{
  // generate a spherical boundary from an interpolated set of triangulated mesh
  b_surf.generate_overlapping_spheres(h, r,BD_l.r_bdry, u, v, w);

  // size the bdry_atoms
  bdry_atoms.set_N(b_surf.get_N_verts());

  std::vector<vec> bdry_coords = b_surf.get_coords();

  vec r0 = vqm.v_new(x0,y0,z0);

  // translate the boundary
  for (size_t i=0; i<bdry_coords.size(); i++)
  {
    bdry_coords[i] = vqm.v_xpy(bdry_coords[i],r0);
  }

  // set the bdry_atoms to the coordinates
  bdry_atoms.set_coords(bdry_coords);
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
  // if (bdry_atoms.get_N() == -1)
  //  {
  //    generate_spherical_bdry(BD_l.r_sphere,0.0,0.0,0.0);
  //  }

  // place the ribo and bdry atoms in individual molecules for convenience
  bdry_atoms.set_mol_id_all(1);
  ribo_atoms.set_mol_id_all(2);

  // set the mono, ribo, and bdry types
  bdry_atoms.set_type_all(1);
  ribo_atoms.set_type_all(2);
  set_mono_types(3);

  // set the mono, ribo, and bdry densities
  mono_atoms.set_density_all(1.0);
  ribo_atoms.set_density_all(1.0);
  bdry_atoms.set_density_all(1.0);

  // set the ellipsoid flags for mono, ribo, and bdry atoms
  mono_atoms.set_ellipsoid_flag_all(1);
  ribo_atoms.set_ellipsoid_flag_all(1);
  bdry_atoms.set_ellipsoid_flag_all(0);

  // prepare the ellipsoid arrays
  mono_ellipsoids.set_N(mono_atoms.get_N());
  ribo_ellipsoids.set_N(ribo_atoms.get_N());

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
void LAMMPS_sys::read_BD_lengths(std::string lengths_filename)
{
  std::fstream lengths_file;

  std::string param_delim, param, val;
  int delim;
  
  std::string line;

  param_delim = "=";

  lengths_file.open(lengths_filename, std::ios::in);

  if (!lengths_file.is_open())
    {
      std::cout << "ERROR: file not opened in read_BD_lengths" << std::endl;
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

		  // std::cout << param << "=" << val << std::endl;

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
int LAMMPS_sys::read_mono_coords(std::string coords_filename, std::string order)
{
  return mono_atoms.read_bin_coords(coords_filename,order,false);
}


// read the monomer quaternions - disallow resizing
int LAMMPS_sys::read_mono_quats(std::string quats_filename, std::string order)
{
  return mono_ellipsoids.read_bin_quats(quats_filename,order,false);
}


// read the ribosome coordinates - allow resizing
int LAMMPS_sys::read_ribo_coords(std::string coords_filename, std::string order)
{
  return ribo_atoms.read_bin_coords(coords_filename,order,true);
}


// read the ribosome quaternions - allow resizing
int LAMMPS_sys::read_ribo_quats(std::string quats_filename, std::string order)
{
  return ribo_ellipsoids.read_bin_quats(quats_filename,order,true);
}


// read the boundary coordinates - allow resizing
int LAMMPS_sys::read_bdry_coords(std::string coords_filename, std::string order)
{
  return bdry_atoms.read_bin_coords(coords_filename,order,true);
}


// write the monomer coordinates
int LAMMPS_sys::write_mono_coords(std::string coords_filename, std::string order)
{
  return mono_atoms.write_bin(coords_filename,order);
}


// write the monomer quaternions
int LAMMPS_sys::write_mono_quats(std::string quats_filename, std::string order)
{
  return mono_ellipsoids.write_bin(quats_filename,order);
}


// write the ribosome coordinates
int LAMMPS_sys::write_ribo_coords(std::string coords_filename, std::string order)
{
  return ribo_atoms.write_bin(coords_filename,order);
}


// write the ribosome quaternions
int LAMMPS_sys::write_ribo_quats(std::string quats_filename, std::string order)
{
  return ribo_ellipsoids.write_bin(quats_filename,order);
}


// write the boundary coordinates
int LAMMPS_sys::write_bdry_coords(std::string coords_filename, std::string order)
{
  return bdry_atoms.write_bin(coords_filename,order);
}


// apply mapping to the monomers
void LAMMPS_sys::apply_mono_mapping(std::vector<std::vector<std::array<int,3>>> map)
{
  int N_trans, N_old, N_new;

  atom_array temp_atoms;
  ellipsoid_array temp_ellipsoids;
  atom t_a;
  ellipsoid t_e;
  vec uy = vqm.v_new(0.0,1.0,0.0);
  vec t_r;
  quat t_q;
  
  N_trans = static_cast<int>(map.size());

  for (int i_trans=0; i_trans<N_trans; i_trans++)
    {

      N_old = mono_atoms.get_N();

      // store current atoms in temporary array
      temp_atoms.set_N(N_old);

      // set the elements of the temp array to the current monomers
      for (int j=0; j<mono_atoms.get_N(); j++)
	{
	  temp_atoms.set_atom(j,mono_atoms.get_atom(j));
	}

      // store current ellipsoids in temporary array
      temp_ellipsoids.set_N(N_old);

      // set the elements of the temp array to the current monomers
      for (int j=0; j<mono_ellipsoids.get_N(); j++)
	{
	  temp_ellipsoids.set_ellipsoid(j,mono_ellipsoids.get_ellipsoid(j));
	}


      // resize the atom and ellipsoid arrays
      N_new = static_cast<int>(map[i_trans].size());
      mono_atoms.set_N(N_new);
      mono_ellipsoids.set_N(N_new);

      std::cout << "N_old = " << N_old << std::endl;
      std::cout << "N_new = " << N_new << std::endl;
      

      for (int i=0; i<N_new; i++)
	{

	  // std::cout << map[i_trans][i][0] << "\t"
	  //      << map[i_trans][i][1] << "\t"
	  //      << map[i_trans][i][2] << std::endl;

	  // store the previous atom and ellipsoid
	  t_a = temp_atoms.get_atom(map[i_trans][i][1]);
	  t_e = temp_ellipsoids.get_ellipsoid(map[i_trans][i][1]);


	  if (map[i_trans][i][2] != 0)
	    {
	      
	      if (map[i_trans][i][2] == 1)
		{

		  t_r = vqm.v_ax(1.0*BD_l.mono_shape.y/2.0,uy);

		}
	      else if(map[i_trans][i][2] == -1)
		{

		  t_r = vqm.v_ax(-1.0*BD_l.mono_shape.y/2.0,uy);

		}

	      t_q = vqm.v_to_q(t_r);
	      t_q = vqm.q_mult(t_q,vqm.q_norm(vqm.q_conj(t_e.q)));
	      t_q = vqm.q_mult(vqm.q_norm(t_e.q),t_q);
	      t_r = vqm.q_to_v(t_q);

	      t_a.r = vqm.v_xpy(t_a.r,t_r);
	      
	    }
	  
	  // set the atom array to the mapped elements
	  mono_atoms.set_atom(map[i_trans][i][0],t_a);

	  // set the ellipsoid array to the mapped elements
	  mono_ellipsoids.set_ellipsoid(map[i_trans][i][0],t_e);
	  
	}
      
    }
}

void LAMMPS_sys::apply_RMF()
{
    // if (!r_s.ellipsoid)
    std::cout << "Setting up the loop topology..." << std::endl;
    prepare_loop_topo();
    std::cout << "Calculating the RMF..." << std::endl;
    // loop over binding regions
    std::vector<binding_region> b_rs = loop_topo.get_regions();
    size_t N_reg = b_rs.size();

    for (size_t i_reg=0; i_reg<N_reg; i_reg++)
    {
        std::cout << "Getting size of binding region..." << std::endl;
        int N = b_rs[i_reg].get_size();
        std::vector<std::vector<vec>> R(N, std::vector<vec>(3));
        std::vector<vec> x(N), xm1(N), xm2(N), xp1(N), xp2(N), t(N);
        vec v1, v2, pL, tL, p, s;
        double c1, c2, t1, t2, dot;

        std::cout << "Making atom arrays for p, s, and t..." << std::endl;
        atom_array p_atoms, s_atoms, t_atoms;
        p_atoms.set_N(N);
        s_atoms.set_N(N);
        t_atoms.set_N(N);

        std::cout << "Creating array of monomer positions..." << std::endl;
        for (int i = 0; i < N; i++)
        {
            x[i] = atoms.get_atom(b_rs[i_reg].get_mono_pos(i)).r;
        }

        std::cout << "Creating array of nearest and next nearest neighbors..." << std::endl;
        for (int i = 0; i < N; i++)
        {
            // special case: circular chromosome
            if (b_rs[i_reg].get_completed()) {
                xm1[i] = x[(i - 1 + N) % N];
                xm2[i] = x[(i - 2 + N) % N];
                xp1[i] = x[(i + 1) % N];
                xp2[i] = x[(i + 2) % N];
            }
            else {
                // special case: start of region
                if (i == 0 || i == 1 || i == N-2 || i == N-1) {
                    xm1[i] = vqm.v_null(); // Won't be needing these.
                    xm2[i] = vqm.v_null();
                    xp1[i] = vqm.v_null();
                    xp2[i] = vqm.v_null();
                }
                else {
                    xm1[i] = x[i - 1];
                    xm2[i] = x[i - 2];
                    xp1[i] = x[i + 1];
                    xp2[i] = x[i + 2];
                }
            }
        }

        std::cout << "Calculating the tangent vectors..." << std::endl;
        for (int i = 0; i < N; i++)
        {
            // special case: circular chromosome
            if (b_rs[i_reg].get_completed()) {
                t[i] = vqm.v_axpy(-1, vqm.v_axpy(-8.0, xp1[i], xp2[i]), vqm.v_axpy(-8.0, xm1[i], xm2[i]));
            }
            else {
                // special case: start of region
                if (i == 0) {
                    t[i] = vqm.v_axpy(-25, x[0], vqm.v_axpy(48, x[1], vqm.v_axpy(-36, x[2], vqm.v_axpy(16, x[3], vqm.v_ax(-3, x[4])))));
                }
                // special case: neighbor of start of region
                else if (i == 1) {
                    t[i] = vqm.v_axpy(-3, x[0], vqm.v_axpy(-10, x[1], vqm.v_axpy(18, x[2], vqm.v_axpy(-6, x[3], vqm.v_ax(1, x[4])))));
                }
                // special case: neighbor of end of region
                else if (i == N - 2) {
                    t[i] = vqm.v_axpy(3, x[N-1], vqm.v_axpy(10, x[N-2], vqm.v_axpy(-18, x[N-3], vqm.v_axpy(6, x[N-4], vqm.v_ax(-1, x[N-5])))));
                }
                // special case: end of region
                else if (i == N - 1) {
                    t[i] = vqm.v_axpy(25, x[N-1], vqm.v_axpy(-48, x[N-2], vqm.v_axpy(36, x[N-3], vqm.v_axpy(-16, x[N-4], vqm.v_ax(3, x[N-5])))));
                }
                else {
                    t[i] = vqm.v_axpy(-1, vqm.v_axpy(-8.0, xp1[i], xp2[i]), vqm.v_axpy(-8.0, xm1[i], xm2[i]));
                }
            }
            t[i] = vqm.v_norm(t[i]);
        }

        std::cout << "Calculating the p and s vectors..." << std::endl;
        for (int i = 0; i < N; i++)
        {
            if (i == 0) { // Initialize the first orientation randomly
                dot = 1.0;
                while (dot > 0.95)
                {
                    t1 = rand() * (2 * M_PI / RAND_MAX); // Random number between 0 and 2*pi
                    t2 = rand() * (M_PI / RAND_MAX);     // Random number between 0 and pi

                    p.x = cos(t1) * sin(t2);
                    p.y = sin(t1) * sin(t2);
                    p.z = cos(t2);

                    dot = vqm.v_dot(p, t[i]);
                }

                p = vqm.v_axpy(-1.0 * dot, t[i], p);
                p = vqm.v_norm(p);
            } else {
                // p = R[i - 1][0];
                p = R[(i - 1 + N) % N][0]; // didn't fix it...
                v1 = vqm.v_axpy(-1.0, xm1[i], x[i]); // instead of x and xp1
                c1 = vqm.v_dot(v1, v1);
                pL = vqm.v_axpy(-2.0 / c1 * vqm.v_dot(v1, p), v1, p);
                tL = vqm.v_axpy(-2.0 / c1 * vqm.v_dot(v1, t[i - 1]), v1, t[i - 1]);
                v2 = vqm.v_axpy(-1.0, tL, t[i]);
                c2 = vqm.v_dot(v2, v2);
                p = vqm.v_axpy(-2.0 / c2 * vqm.v_dot(v2, pL), v2, pL);
                p = vqm.v_norm(p);
            }

            s = vqm.v_cross(t[i], p);
            s = vqm.v_norm(s);

            R[i][0] = p;
            R[i][1] = s;
            R[i][2] = t[i];

            // set orientation of DNA monomer
            int reg_pos = b_rs[i_reg].get_mono_pos(i);
            mono_ellipsoids.set_quat(reg_pos, vqm.R_to_q(R[i]));

            vec s2 = vqm.q_to_v(mono_ellipsoids.get_ellipsoid(i).q);
            s2 = vqm.v_norm(s2); // normalize!

            vec p_coord, s_coord, t_coord;
            p_coord = vqm.v_axpy(BD_l.mono_shape.x, p, x[i]);
            s_coord = vqm.v_axpy(BD_l.mono_shape.x, s2, x[i]);
            t_coord = vqm.v_axpy(BD_l.mono_shape.x, t[i], x[i]);
            // std::cout << "Setting orientation atom coordinates..." << std::endl;
            p_atoms.set_coord(i, p_coord);
            s_atoms.set_coord(i, s_coord);
            t_atoms.set_coord(i, t_coord);
        }
        std::cout << "Writing p, s, and t atom coordinates..." << std::endl;
        p_atoms.write_xyz("/home/andrew/Data/btree_chromo/p.xyz");
        s_atoms.write_xyz("/home/andrew/Data/btree_chromo/s.xyz");
        t_atoms.write_xyz("/home/andrew/Data/btree_chromo/t.xyz");
    }
    std::cout << "Done calculating RMF for all binding regions..." << std::endl;
}

// write the system to a data file
void LAMMPS_sys::write_data(std::string data_filename)
{
  // internal_btree.print_tree();

  // prepare_test_data();


  // finalize the system before printing
  finalize_system();
  
  // begin writing data file
  
  std::fstream data_file;

  data_file.open(data_filename, std::ios::out);

  if (!data_file.is_open())
    {
      std::cout << "ERROR: file not opened in write_data" << std::endl;
    }
  else
    {
      if (r_s.ellipsoid) {
          // write system summary
          data_file << "# LAMMPS data file for replicating chromosomes formed of rigid body monomers\n" << std::endl;

          data_file << atoms.get_N() << "\t\tatoms" << std::endl;
          data_file << N_atom_types << "\t\tatom types" << std::endl;
          data_file << ellipsoids.get_N() << "\t\tellipsoids" << std::endl;
          data_file << bonds.get_N() << "\t\tbonds" << std::endl;
          data_file << N_bond_types << "\t\tbond types" << std::endl;
          data_file << angles.get_N() << "\t\tangles" << std::endl;
          data_file << N_angle_types << "\t\tangle types" << std::endl;

          data_file << "\n" << std::endl;

          data_file << bbox.r_min.x << "\t" << bbox.r_max.x << "\txlo xhi" << std::endl;
          data_file << bbox.r_min.y << "\t" << bbox.r_max.y << "\tylo yhi" << std::endl;
          data_file << bbox.r_min.z << "\t" << bbox.r_max.z << "\tzlo zhi" << std::endl;

          // write atom information
          atoms.write(data_file, r_s.ellipsoid);

          // write ellipsoid information
          ellipsoids.write(data_file);

          // write bond information
          bonds.write(data_file);

          // write angle information
          angles.write(data_file);


          data_file.close();
      }
      else {
          // write system summary
          data_file << "# LAMMPS data file for replicating chromosomes formed of rigid body monomers\n" << std::endl;

          data_file << atoms.get_N() << "\t\tatoms" << std::endl;
          data_file << N_atom_types << "\t\tatom types" << std::endl;
          data_file << bonds.get_N() << "\t\tbonds" << std::endl;
          data_file << N_bond_types << "\t\tbond types" << std::endl;
          data_file << angles.get_N() << "\t\tangles" << std::endl;
          data_file << N_angle_types << "\t\tangle types" << std::endl;

          data_file << "\n" << std::endl;

          data_file << bbox.r_min.x << "\t" << bbox.r_max.x << "\txlo xhi" << std::endl;
          data_file << bbox.r_min.y << "\t" << bbox.r_max.y << "\tylo yhi" << std::endl;
          data_file << bbox.r_min.z << "\t" << bbox.r_max.z << "\tzlo zhi" << std::endl;

          // write atom information
          atoms.write(data_file, r_s.ellipsoid);

          // write bond information
          bonds.write(data_file);

          // write angle information
          angles.write(data_file);


          data_file.close();
      }

    }
  
}


// write an xyz file with the monomer atoms
void LAMMPS_sys::write_mono_xyz(std::string data_filename)
{
  mono_atoms.write_xyz(data_filename);
}

// write an xyz file with the monomer atoms
void LAMMPS_sys::write_bdry_xyz(std::string data_filename)
{
  bdry_atoms.write_xyz(data_filename);
}

// write an xyz file with the orientation vector endpoints as atoms
void LAMMPS_sys::write_mono_orientation_xyz(std::string data_filename)
{
    // std::cout << "Setting up the loop topology..." << std::endl;
    prepare_loop_topo();

    std::vector<binding_region> b_rs = loop_topo.get_regions();
    size_t N_reg = b_rs.size();
    int N = 0;

    // create a first loop where you iterate over the binding regions and calculate the total size of p_atoms, N
    for (size_t i_reg=0; i_reg<N_reg; i_reg++)
    {
        N += b_rs[i_reg].get_size();
    }

    // define p_atoms here and then set N
    std::cout << "Making atom arrays for p..." << std::endl;
    atom_array p_atoms;
    p_atoms.set_N(N);

    // create a second loop over the binding regions
    //  within the loop over binding regions, then iterate over the atoms within that binding region
    //      set the coordinates within p_atoms as you do that iteration
    for (size_t i_reg=0; i_reg<N_reg; i_reg++)
    {

        std::cout << "Getting size of binding region..." << std::endl;
        int M = b_rs[i_reg].get_size();

        // Write the RMF vector endpoints to xyz files
        // Here, we create monomers 3.4 nm away from origin of RMF (DNA monomer position) in directions of p, s, and t

        // iterate over the monomers within the ith binding region
        for (int i = 0; i < M; i++)
        {
            // std::cout << "Getting monomer coordinates..." << std::endl;
            vec x = atoms.get_atom(b_rs[i_reg].get_mono_pos(i)).r;
            // std::cout << "Getting orientations..." << std::endl;
            vec p = vqm.q_to_R(mono_ellipsoids.get_ellipsoid(i).q)[0];
            p = vqm.v_norm(p); // normalize!
            // std::cout << "Calculating orientation atom coordinates..." << std::endl;
            vec p_coord;
            p_coord = vqm.v_axpy(BD_l.mono_shape.x, p, x);
            // std::cout << "Setting orientation atom coordinates..." << std::endl;
            p_atoms.set_coord(i, p_coord);
        }

    }

    std::cout << "Writing orientation atom coordinates..." << std::endl;
    p_atoms.write_xyz(data_filename);
}


// get the total number of atoms in the sytem
int LAMMPS_sys::get_N_total()
{
  int N_total = 0;
  int N_mono = get_N_mono();
  int N_ribo = get_N_ribo();
  int N_bdry = get_N_bdry();

  if (N_mono > 0) N_total += N_mono;
  if (N_ribo > 0) N_total += N_ribo;
  if (N_bdry > 0) N_total += N_bdry;
  
  return N_total;
}


// get the number of mono atoms
int LAMMPS_sys::get_N_mono()
{
  return mono_atoms.get_N();
}


// get the number of ribo atoms
int LAMMPS_sys::get_N_ribo()
{
  return ribo_atoms.get_N();
}


// get the number of bdry atoms
int LAMMPS_sys::get_N_bdry()
{
  return bdry_atoms.get_N();
}


// set the total number of atoms
void LAMMPS_sys::set_N_total(int N)
{
  atoms.set_N(N);
}


// set the number of mono atoms
void LAMMPS_sys::set_N_mono(int N_mono)
{
  mono_atoms.set_N(N_mono);
}


// set the number of ribo atoms
void LAMMPS_sys::set_N_ribo(int N_ribo)
{
  ribo_atoms.set_N(N_ribo);
}


// set the number of bdry atoms
void LAMMPS_sys::set_N_bdry(int N_bdry)
{
  bdry_atoms.set_N(N_bdry);
}


// set the total number of ellipsoids
void LAMMPS_sys::set_N_total_ellipsoids(int N)
{
  ellipsoids.set_N(N);
}


// set the number of mono ellipsoids
void LAMMPS_sys::set_N_mono_ellipsoids(int N_mono)
{
  mono_ellipsoids.set_N(N_mono);
}


// set the number of ribo ellipsoids
void LAMMPS_sys::set_N_ribo_ellipsoids(int N_ribo)
{
  ribo_ellipsoids.set_N(N_ribo);
}


// set the total coordinate array for the system
void LAMMPS_sys::set_coords_arr_total(double *&x, std::string order)
{
  atoms.set_coords_arr(x,order);
}


// set the total coordinate array for the system
void LAMMPS_sys::set_quats_arr_total(double *&q, std::string order)
{
  ellipsoids.set_quats_arr(q,order);
}


// sync the subarrays with the total arrays
void LAMMPS_sys::sync_subarrays()
{
  int N_mono, N_ribo, N_bdry;

  N_mono = std::max(0,mono_atoms.get_N());
  N_ribo = std::max(0,ribo_atoms.get_N());
  N_bdry = std::max(0,bdry_atoms.get_N());
  
  // copy to mono_atoms
  for (int i=0; i<N_mono; i++)
    {
      mono_atoms.set_atom(i,atoms.get_atom(i));
    }

  // copy to ribo_atoms
  for (int i=0; i<N_ribo; i++)
    {
      ribo_atoms.set_atom(i,atoms.get_atom(i+N_mono));
    }

  // copy to bdry_atoms
  for (int i=0; i<N_bdry; i++)
    {
      bdry_atoms.set_atom(i,atoms.get_atom(i+N_mono+N_ribo));
    }

  // copy to mono_ellipsoids
  for (int i=0; i<N_mono; i++)
    {
      mono_ellipsoids.set_ellipsoid(i,ellipsoids.get_ellipsoid(i));
    }

  // copy to ribo_ellipsoids
  for (int i=0; i<N_ribo; i++)
    {
      ribo_ellipsoids.set_ellipsoid(i,ellipsoids.get_ellipsoid(i+N_mono));
    }
  
}

void LAMMPS_sys::prepare_loop_topo()
{
    std::vector<std::string> leaves = internal_btree.get_leaves();
    std::vector<theta_topo> leaf_topos;

    for (std::string leaf : leaves)
    {
        leaf_topos.push_back(internal_btree.get_leaf_topo(leaf));
    }

    int *t, N;

    N = get_N_mono();
    t = new int[N];

    for (int i=0; i<N; i++)
    {
        t[i] = mono_atoms.get_atom(i).type;
    }

    // prepare the possible binding regions
    loop_topo.prepare_binding_regions(leaves,leaf_topos,t);

    //std::vector<binding_region> regions = loop_topo.get_regions();

    //std::cout << "BINDING REGIONS" << std::endl;
    //for (binding_region b_r : regions)
    //{
    //std::cout << b_r.get_leaf() << " " << b_r.get_size() << std::endl;
    //}

    delete[] t;
}

// initialize the loop topology
void LAMMPS_sys::initialize_loop_topo(int N_loops, double radius)
{

  prepare_loop_topo();

  loop_topo.set_coords(get_mono_coords());

  // initialize_voxel_grid(250.0);

  // loop_topo.set_anchor_candidates(sample_anchors(N_loops, 50, radius));

  // loop_topo.initialize_loops_on_periphery(0,l_sys_p.min_dist, radius); // set some on periphery gauranteed
  int N_bound = int(N_loops * l_sys_p.p_on /(l_sys_p.p_off + l_sys_p.p_on));
  int N_unbound = N_loops - N_bound;
  loop_topo.initialize_loops_interior(N_bound, N_unbound, l_sys_p.min_dist, radius);
  // loop_topo.initialize_loops_conc(N_loops,l_sys_p.min_dist); // set the rest in the low conc. voxels
}


// update the loop topology
void LAMMPS_sys::update_loop_topo(double threshold)
{

  loop_topo.set_coords(get_mono_coords());

  loop_topo.update_loops(l_sys_p.ext_max,
			 l_sys_p.min_dist,
			 l_sys_p.p_unbinding,
       l_sys_p.p_off,
       l_sys_p.p_on,
			 l_sys_p.r_g,
			 threshold);

  
}


// get the loops and convert them into bonds
std::vector<bond> LAMMPS_sys::get_loop_bonds()
{
  bond loop_bond;
  std::vector<bond> loop_bonds;
  int id = bonds.get_N() + 1;
  
  std::vector<loop> loops = loop_topo.get_loops();

  for (loop l : loops)
  {
      if (l.get_a_bound())
	  {
      loop_bond.id = id;
      loop_bond.type = 2;
      loop_bond.i = l.get_a() + 1;
      loop_bond.j = l.get_h() + 1;
      loop_bonds.push_back(loop_bond);
      id += 1;
      }
  }
  return loop_bonds;
}


// set the loop sys parameters
void LAMMPS_sys::set_loop_sys_params(loop_sys_params &l_sys_p)
{
  this->l_sys_p = l_sys_p;
  // set hinge update parameters
  loop_topo.set_step_dist(l_sys_p.family,l_sys_p.ext_avg,l_sys_p.ext_max);
  // set binding/unbinding parameters
}


// getter for monomer coordinates
std::vector<vec> LAMMPS_sys::get_mono_coords()
{
  return mono_atoms.get_coords();
}

void LAMMPS_sys::initialize_voxel_grid(double voxel_size) {
  // Initialize voxel grid with specific values
  std::cout << "Initializing voxel grid..." << std::endl;
  voxel_grid.initialize_voxels(voxel_size, bbox.r_min, bbox.r_max);
  // Assign each monomer to a voxel
  std::cout << "Starting a loop binning monomers into voxels..." << std::endl;
  for (int i = 0; i < mono_atoms.get_N(); ++i) {
    const vec& position = mono_atoms.get_atom(i).r;
    voxel_grid.add_monomer(position, i);
  }
  std::cout << "Finished binning monomers..." << std::endl;
}

std::vector<std::tuple<int, int, int>> LAMMPS_sys::find_low_concentration_voxels(int threshold) {
  std::vector<std::tuple<int, int, int>> low_voxels;

  // Iterate over the 3D grid dimensions
  for (int iz = 0; iz < voxel_grid.dim_z; iz++) {
    for (int iy = 0; iy < voxel_grid.dim_y; iy++) {
      for (int ix = 0; ix < voxel_grid.dim_x; ix++) {
        int index = voxel_grid.get_index(ix, iy, iz);

        // Calculate "concentration" as the size of the vector at the index
        if (!voxel_grid.voxel_grid[index].empty() && voxel_grid.voxel_grid[index].size() <= static_cast<size_t>(threshold)) {
          low_voxels.emplace_back(ix, iy, iz);
        }
      }
    }
  }
  std::cout << "Found " << low_voxels.size() << " voxels that are below concentration threshold. " << std::endl;
  return low_voxels;
}

std::vector<std::tuple<int, int, int>> LAMMPS_sys::get_top_low_density_voxels(int N) {
    std::vector<std::tuple<int, int, int, size_t>> all_voxels_with_density;

    // Iterate over the 3D grid dimensions
    for (int iz = 0; iz < voxel_grid.dim_z; iz++) {
        for (int iy = 0; iy < voxel_grid.dim_y; iy++) {
            for (int ix = 0; ix < voxel_grid.dim_x; ix++) {
                int index = voxel_grid.get_index(ix, iy, iz);

                // Include only non-empty voxels
                if (!voxel_grid.voxel_grid[index].empty()) {
                    size_t density = voxel_grid.voxel_grid[index].size();
                    all_voxels_with_density.emplace_back(ix, iy, iz, density);
                }
            }
        }
    }

    // Sort voxels by density in ascending order
    std::sort(all_voxels_with_density.begin(), all_voxels_with_density.end(),
              [](const std::tuple<int, int, int, size_t>& a, const std::tuple<int, int, int, size_t>& b) {
                  return std::get<3>(a) < std::get<3>(b);
              });

    // Prepare the result vector with the top N lowest density voxels
    std::vector<std::tuple<int, int, int>> low_density_voxels;
    for (int i = 0; i < std::min(N, static_cast<int>(all_voxels_with_density.size())); ++i) {
        low_density_voxels.emplace_back(std::get<0>(all_voxels_with_density[i]),
                                        std::get<1>(all_voxels_with_density[i]),
                                        std::get<2>(all_voxels_with_density[i]));
    }

    std::cout << "Found " << low_density_voxels.size() << " voxels with the lowest densities." << std::endl;
    return low_density_voxels;
}

std::vector<int> LAMMPS_sys::sample_anchors(int num_anchors, int num_voxels, double radius) {
  // Find low-concentration voxels
  auto low_voxels = find_low_concentration_voxels(num_voxels);

  // Random number generators
  std::mt19937 rand_eng(std::random_device{}());
  std::uniform_int_distribution<int> voxel_dist(0, low_voxels.size() - 1);

  std::vector<int> anchors;

  while (anchors.size() < static_cast<size_t>(num_anchors)) {
    // Randomly pick a low-concentration voxel
    auto voxel_tuple = low_voxels[voxel_dist(rand_eng)];
    int ix = std::get<0>(voxel_tuple);
    int iy = std::get<1>(voxel_tuple);
    int iz = std::get<2>(voxel_tuple);
    int voxel_index = voxel_grid.get_index(ix, iy, iz);
    const auto& candidates = voxel_grid.voxel_grid[voxel_index];

    if (!candidates.empty()) {
      std::uniform_int_distribution<int> candidate_dist(0, candidates.size() - 1);
      // Randomly pick an anchor from the voxel
      int anchor_index = candidates[candidate_dist(rand_eng)];

      // Ensure the anchor is not already in anchor_candidates
      double mono_dist = sqrt(vqm.v_L2(mono_atoms.get_atom(anchor_index).r));
      if (mono_dist < radius) {
        if (std::find(anchors.begin(), anchors.end(), anchor_index) == anchors.end()) {
          std::cout << "Adding to anchor candidates monomer " << anchor_index << std::endl;
          anchors.push_back(anchor_index); // Update anchor_candidates if needed
        } else {
          std::cout << "Anchor " << anchor_index << " is already in anchor candidates, skipping...\n";
        }
      }
    }
  }
  std::cout << "Finished finding all " << num_anchors << " anchors from the " << num_voxels << "lowest density voxels." << std::endl;
  return anchors;
  }


  void LAMMPS_sys::set_boundary_radius(double radius) {
    this->boundary_radius = radius;
  }

  double LAMMPS_sys::get_boundary_radius() {
    return boundary_radius;
  }

  void LAMMPS_sys::update_loop_binding() {
    // update the regions of the existing loops (anchors and hinges. they might lie in different regions, presently...)
    loop_topo.update_anchor_hinge_regions();
  }

  void LAMMPS_sys::add_ribos(int N_ribos, double radius) {
    ribo_atoms.resize(ribo_atoms.get_N()+N_ribos);
    for (int i = 0; i < N_ribos; ++i) {
      double x, y, z;
	  while (true) {
    	x = 2.0 * ((double)std::rand() / RAND_MAX) - 1.0; // Generate x in [-1, 1]
   	 	y = 2.0 * ((double)std::rand() / RAND_MAX) - 1.0; // Generate y in [-1, 1]
    	z = 2.0 * ((double)std::rand() / RAND_MAX) - 1.0; // Generate z in [-1, 1]

    	double magnitude_squared = x * x + y * y + z * z;
    	if (magnitude_squared <= 1.0) {
        	// Scale to the desired radius
        	x *= radius;
        	y *= radius;
        	z *= radius;
        	break;
    	}
	}
	// Output the generated coordinates
	std::cout << "Ribosome placed inside of radius " << radius << " at coordinate: (" << x << ", " << y << ", " << z << ")\n";
    ribo_atoms.set_coord(i, vqm.v_new(x, y, z));
    }
    ribo_atoms.set_type_all(2);
    ribo_atoms.set_mol_id_all(2);
    ribo_atoms.set_density_all(1.0);
    ribo_atoms.set_ellipsoid_flag_all(1);
  }

  void LAMMPS_sys::add_loops(int N_loops, double radius) {
    loop_topo.set_coords(get_mono_coords());
    loop_topo.initialize_loops_interior(N_loops, 0, l_sys_p.min_dist, radius);
  }
