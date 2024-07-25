// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_brownian_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixBrownianKokkos<DeviceType>::FixBrownianKokkos(LAMMPS *lmp, int narg, char **arg) :
        FixBrownian(lmp, narg, arg),rand_pool(seed + comm->me)
{
    kokkosable = 1;
    fuse_integrate_flag = 0; //???
    atomKK = (AtomKokkos *) atom;
    execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

    datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK | TYPE_MASK;
    datamask_modify = X_MASK | V_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixBrownianKokkos<DeviceType>::init()
{
    FixBrownian::init();
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

template<class DeviceType>
void FixBrownianKokkos<DeviceType>::initial_integrate(int /*vflag*/)
{
    atomKK->sync(execution_space,datamask_read);
    // atomKK->modified(execution_space,datamask_modify);

    x = atomKK->k_x.view<DeviceType>();
    v = atomKK->k_v.view<DeviceType>();
    f = atomKK->k_f.view<DeviceType>();
    type = atomKK->k_type.view<DeviceType>();
    mask = atomKK->k_mask.view<DeviceType>();
    int nlocal = atomKK->nlocal;
    if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

    FixBrownianKokkosInitialIntegrateFunctor<DeviceType> functor(this);
    Kokkos::parallel_for(nlocal, functor);

    atomKK->modified(execution_space,datamask_modify);

}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixBrownianKokkos<DeviceType>::initial_integrate_item(int i) const
{
    double dx, dy, dz;

    if (mask[i] & groupbit) {
        rand_type rand_gen = rand_pool.get_state();
        dx = dt * (g1 * f(i,0) + g2 * (rand_gen.drand() - 0.5));
        dy = dt * (g1 * f(i,1) + g2 * (rand_gen.drand() - 0.5));
        dz = dt * (g1 * f(i,2) + g2 * (rand_gen.drand() - 0.5));

        x(i,0) += dx;
        v(i,0) = dx / dt;

        x(i,1) += dy;
        v(i,1) = dy / dt;

        x(i,2) += dz;
        v(i,2) = dz / dt;
        rand_pool.free_state(rand_gen);
    }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixBrownianKokkos<DeviceType>::cleanup_copy()
{
    gamma_t_inv = nullptr;
    gamma_t_invsqrt = nullptr;
    gamma_r_inv = nullptr;
    gamma_r_invsqrt = nullptr;
    dipole_body = nullptr;
    rng = nullptr;
    id = style = nullptr;
    vatom = nullptr;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
    template class FixBrownianKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
    template class FixBrownianKokkos<LMPHostType>;
#endif
}

