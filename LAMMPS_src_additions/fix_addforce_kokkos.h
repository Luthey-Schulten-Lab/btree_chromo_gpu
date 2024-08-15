/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(addforce/kk,FixAddForceKokkos<LMPDeviceType>);
FixStyle(addforce/kk/device,FixAddForceKokkos<LMPDeviceType>);
FixStyle(addforce/kk/host,FixAddForceKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_ADD_FORCE_KOKKOS_H
#define LMP_FIX_ADD_FORCE_KOKKOS_H

#include "fix_setforce.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct s_double_3a {
  double d0, d1, d2;
  KOKKOS_INLINE_FUNCTION
  s_double_3a() {
    d0 = d1 = d2 = 0.0;
  }
  KOKKOS_INLINE_FUNCTION
  s_double_3a& operator+=(const s_double_3a &rhs) {
    d0 += rhs.d0;
    d1 += rhs.d1;
    d2 += rhs.d2;
    return *this;
  }
};
typedef s_double_3a double_3a;

struct TagFixAddForceConstant{};

struct TagFixAddForceNonConstant{};

template<class DeviceType>
class FixAddForceKokkos : public FixSetForce {
 public:
  typedef DeviceType device_type;
  typedef double_3a value_type;
  typedef ArrayTypes<DeviceType> AT;

  FixAddForceKokkos(class LAMMPS *, int, char **);
  ~FixAddForceKokkos() override;
  void init() override;
  void post_force(int) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAddForceConstant, const int&, double_3a&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAddForceNonConstant, const int&, double_3a&) const;

 private:
  DAT::tdual_ffloat_2d k_sforce;
  typename AT::t_ffloat_2d_randomread d_sforce;
  typename AT::t_int_1d d_match;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread mask;
};

}

#endif
#endif

