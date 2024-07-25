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
FixStyle(brownian/kk,FixBrownianKokkos<LMPDeviceType>);
FixStyle(brownian/kk/device,FixBrownianKokkos<LMPDeviceType>);
FixStyle(brownian/kk/host,FixBrownianKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_BROWNIAN_KOKKOS_H
#define LMP_FIX_BROWNIAN_KOKKOS_H

#include "fix_brownian.h"
#include "kokkos_type.h"
#include "kokkos_base.h"
#include "Kokkos_Random.hpp"

namespace LAMMPS_NS {

    template<class DeviceType>
    class FixBrownianKokkos;

    template <class DeviceType>
    struct FixBrownianKokkosInitialIntegrateFunctor;

    template<class DeviceType>
    class FixBrownianKokkos : public FixBrownian, public KokkosBase {
    public:
        FixBrownianKokkos(class LAMMPS *, int, char **);

        void cleanup_copy();
        void init() override;
        void initial_integrate(int) override;

        KOKKOS_INLINE_FUNCTION
        void initial_integrate_item(int) const;

    private:


        typename ArrayTypes<DeviceType>::t_x_array x;
        typename ArrayTypes<DeviceType>::t_v_array v;
        typename ArrayTypes<DeviceType>::t_f_array_const f;
        typename ArrayTypes<DeviceType>::t_int_1d type;
        typename ArrayTypes<DeviceType>::t_int_1d mask;

        Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
        typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;

    };

    template <class DeviceType>
    struct FixBrownianKokkosInitialIntegrateFunctor  {
        typedef DeviceType  device_type ;
        FixBrownianKokkos<DeviceType> c;

        FixBrownianKokkosInitialIntegrateFunctor(FixBrownianKokkos<DeviceType>* c_ptr):
                c(*c_ptr) {c.cleanup_copy();};
        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const {
            c.initial_integrate_item(i);
        }
    };

}

#endif
#endif

