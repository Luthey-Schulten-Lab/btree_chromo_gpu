#include <vec_quat_manipulator.hpp>

// constructor
vec_quat_manipulator::vec_quat_manipulator()
{
}

// destructor
vec_quat_manipulator::~vec_quat_manipulator()
{
}

// vector functions

vec vec_quat_manipulator::v_null()
{
  vec v;
  v.x = 0.0;
  v.y = 0.0;
  v.z = 0.0;
  return v;
}

vec vec_quat_manipulator::v_new(double x, double y, double z)
{
  vec v;
  v.x = x;
  v.y = y;
  v.z = z;
  return v;
}

double vec_quat_manipulator::v_dot(vec v, vec w)
{
  return v.x*w.x + v.y*w.y + v.z*w.z;
}

double vec_quat_manipulator::v_L2(vec v)
{
  return v_dot(v,v);
}

vec vec_quat_manipulator::v_norm(vec v)
{
  return v_ax(1.0/sqrt(v_L2(v)),v);
}

vec vec_quat_manipulator::v_inv(vec v)
{
  return v_ax(-1.0,v);
}

vec vec_quat_manipulator::v_ax(double a, vec v)
{
  vec ax;
  ax.x = a*v.x;
  ax.y = a*v.y;
  ax.z = a*v.z;
  return ax;
}

vec vec_quat_manipulator::v_xpy(vec v, vec w)
{
  vec xpy;
  xpy.x = v.x + w.x;
  xpy.y = v.y + w.y;
  xpy.z = v.z + w.z;
  return xpy;
}

vec vec_quat_manipulator::v_axpy(double a, vec v, vec w)
{ 
  return v_xpy(v_ax(a,v),w);
}

vec vec_quat_manipulator::v_cross(vec v, vec w)
{
  vec c;
  c.x = v.y*w.z - v.z*w.y;
  c.y = v.z*w.x - v.x*w.z;
  c.z = v.x*w.y - v.y*w.x;
  return c;
}

vec vec_quat_manipulator::v_linterp(double s, vec v, vec w)
{
  return v_axpy(s,v_xpy(v,v_inv(w)),w);
}

void vec_quat_manipulator::v_array_from_vector(double **&x, int &N, std::vector<vec> &vs)
{

  N = static_cast<int>(vs.size());
  
  x = new double*[N];
  for (int i=0; i<N; i++)
    {
      x[i] = new double[3];
      x[i][0] = vs[i].x;
      x[i][1] = vs[i].y;
      x[i][2] = vs[i].z;
    }
  
}

void vec_quat_manipulator::v_vector_from_array(std::vector<vec> &vs, double **&x, int &N)
{

  vs.clear();
  
  for (int i=0; i<N; i++)
    {
      vs.push_back(v_new(x[i][0],x[i][1],x[i][2]));
    }
  
}

// quaternion functions

quat vec_quat_manipulator::q_null()
{
  quat q;
  q.w = 0.0;
  q.v = v_null();
  return q;
}

double vec_quat_manipulator::q_dot(quat q, quat p)
{
  return q.w*p.w + v_dot(q.v,p.v);
}

double vec_quat_manipulator::q_L2(quat q)
{
  return q_dot(q,q);
}

quat vec_quat_manipulator::q_norm(quat q)
{
  return q_ax(1.0/sqrt(q_L2(q)),q);
}

quat vec_quat_manipulator::q_ax(double a, quat q)
{
  quat ax;
  ax.w = a*q.w;
  ax.v = v_ax(a,q.v);
  return ax;
}

quat vec_quat_manipulator::q_xpy(quat q, quat p)
{
  quat xpy;
  xpy.w = q.w + p.w;
  xpy.v = v_xpy(q.v,p.v);
  return xpy;
}

quat vec_quat_manipulator::q_axpy(double a, quat q, quat p)
{ 
  return q_ax(a,q_xpy(q,p));
}

quat vec_quat_manipulator::q_mult(quat q, quat p)
{
  quat qmp;
  qmp.w = q.w*p.w - v_dot(q.v,p.v);
  qmp.v = v_cross(q.v,p.v);
  qmp.v = v_axpy(p.w,q.v,qmp.v);
  qmp.v = v_axpy(q.w,p.v,qmp.v);
  return qmp;
}

quat vec_quat_manipulator::q_conj(quat q)
{
  quat qc = q;
  qc.v = v_ax(-1.0,qc.v);
  return qc;
}

quat vec_quat_manipulator::q_inv(quat q)
{
  return q_ax(1.0/q_L2(q),q_conj(q));
}


// conversions

// convert a vector to a quaternion
quat vec_quat_manipulator::v_to_q(vec v)
{
  quat q = q_null();
  q.v = v;
  return q;
}

// convert a quaternion to a vector
vec vec_quat_manipulator::q_to_v(quat q)
{
  return q.v;
}

quat vec_quat_manipulator::R_to_q(std::vector<vec> R)
{
    quat q;
    vec p=R[0], s=R[1], t=R[2];
    double trace = p.x + s.y + t.z;
    if (trace > 0) {
        double S = std::sqrt(trace + 1.0) * 2; // S=4*q.w
        q.w = 0.25 * S;
        q.v.x = (t.y - s.z) / S;
        q.v.y = (p.z - t.x) / S;
        q.v.z = (s.x - p.y) / S;
    } else {
        if (p.x > s.y && p.x > t.z) {
            double S = std::sqrt(1.0 + p.x - s.y - t.z) * 2; // S=4*q.v.x
            q.w = (t.y - s.z) / S;
            q.v.x = 0.25 * S;
            q.v.y = (p.y + s.x) / S;
            q.v.z = (p.z + t.x) / S;
        } else if (s.y > t.z) {
            double S = std::sqrt(1.0 + s.y - p.x - t.z) * 2; // S=4*q.v.y
            q.w = (p.z - t.x) / S;
            q.v.x = (p.y + s.x) / S;
            q.v.y = 0.25 * S;
            q.v.z = (s.z + t.y) / S;
        } else {
            double S = std::sqrt(1.0 + t.z - p.x - s.y) * 2; // S=4*q.v.z
            q.w = (s.x - p.y) / S;
            q.v.x = (p.z + t.x) / S;
            q.v.y = (s.z + t.y) / S;
            q.v.z = 0.25 * S;
        }
    }

    return q;
}

std::vector<vec> vec_quat_manipulator::q_to_R(quat q)
{
    double q0 = q.w, q1 = q.v.x, q2 = q.v.y, q3 = q.v.z;
    vec p, s, t; // aka f, v, u
    std::vector<vec> R;

    // fI vector (first column of rotation matrix)
    p.x = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
    p.y = 2.0 * (q1 * q2 + q0 * q3);
    p.z = 2.0 * (q1 * q3 - q0 * q2);

    // vI vector (second column of rotation matrix)
    s.x = 2.0 * (q1 * q2 - q0 * q3);
    s.y = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
    s.z = 2.0 * (q2 * q3 + q0 * q1);

    // uI vector (third column of rotation matrix)
    t.x = 2.0 * (q1 * q3 + q0 * q2);
    t.y = 2.0 * (q2 * q3 - q0 * q1);
    t.z = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;

    R[0] = p;
    R[1] = s;
    R[2] = t;

    return R;
}
