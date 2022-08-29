#include <LAMMPS_sys/ellipsoid_array.hpp>

// constructor
ellipsoid_array::ellipsoid_array()
{
  N = -1;
  ellipsoids = nullptr;
}

// destructor
ellipsoid_array::~ellipsoid_array()
{
  destroy_ellipsoids();
}

// set the number of elements
void ellipsoid_array::set_N(int N)
{
  this->N = N;
  initialize_ellipsoids();
}

// get the number of elements
int ellipsoid_array::get_N()
{
  return N;
}

// initialize the array
void ellipsoid_array::initialize_ellipsoids()
{
  destroy_ellipsoids();

  if (N > 0)
    {

      ellipsoids = new ellipsoid[N];

      for (int i=0; i<N; i++)
	{
	  // atom id
	  ellipsoids[i].id = i + 1;
	  // shape
	  ellipsoids[i].s.x = -1.0;
	  ellipsoids[i].s.y = -1.0;
	  ellipsoids[i].s.z = -1.0;
	  // orientation
	  ellipsoids[i].q.w = -1.0;
	  ellipsoids[i].q.i = -1.0;
	  ellipsoids[i].q.j = -1.0;
	  ellipsoids[i].q.k = -1.0;
	}
      
    }
  
}

// destroy the array
void ellipsoid_array::destroy_ellipsoids()
{
  if (ellipsoids != nullptr)
    {
      delete[] ellipsoids;
      ellipsoids = nullptr;
    }
}


// setter for elements
void ellipsoid_array::set_ellipsoid(int i, ellipsoid e)
{
  ellipsoids[i] = e;
}


// getter for elements
ellipsoid ellipsoid_array::get_ellipsoid(int i)
{
  return ellipsoids[i];
}


// set the minimum atom id and range accordingly [min,min+N)
void ellipsoid_array::set_min_id(int min_id)
{
  for (int i=0; i<N; i++)
    {
      ellipsoids[i].id = min_id + i;
    }
}


// set the shape of all ellipsoids
void ellipsoid_array::set_shape_all(vec shape)
{
  for (int i=0; i<N; i++)
    {
      ellipsoids[i].s = shape;
    }
}


// write to stream
void ellipsoid_array::write(fstream &data_file)
{
  data_file << "\nEllipsoids # atom-ID shapex shapey shapez quatw quati quatj quatk\n" << endl;

  for (int i=0; i<N; i++)
    {        
      data_file << ellipsoids[i].id << "\t"
		<< ellipsoids[i].s.x << "\t"
		<< ellipsoids[i].s.y << "\t"
		<< ellipsoids[i].s.z << "\t"
		<< ellipsoids[i].q.w << "\t"
		<< ellipsoids[i].q.i << "\t"
		<< ellipsoids[i].q.j << "\t"
		<< ellipsoids[i].q.k << endl;
    }
  
}
