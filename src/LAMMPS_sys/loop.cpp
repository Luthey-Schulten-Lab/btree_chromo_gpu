#include <LAMMPS_sys/loop.hpp>

// constructor
loop::loop(int a_region, int h_region)
{
  this->a_region = a_region;
  this->h_region = h_region;
  this->a = -1;
  this->h = -1;
  this->d = 0;
}


// destructor
loop::~loop()
{
}


// set a
void loop::set_a(int a)
{
  this->a = a;
}


// set h
void loop::set_h(int h)
{
  this->h = h;
}


// set d
void loop::set_d(int d)
{
  this->d = d;
}


// set a_region
void loop::set_a_region(int a_region)
{
  this->a_region = a_region;
}


// set h_region
void loop::set_h_region(int h_region)
{
  this->h_region = h_region;
}


// get a
int loop::get_a()
{
  return a;
}


// get h
int loop::get_h()
{
  return h;
}


// get d
int loop::get_d()
{
  return d;
}


// get a_region
int loop::get_a_region()
{
  return a_region;
}


// get h_region
int loop::get_h_region()
{
  return h_region;
}
