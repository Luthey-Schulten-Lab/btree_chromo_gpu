/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Originally modified from CG-DNA/fix_nve_dotc_langevin.cpp.

   Contributing author: Sam Cameron (University of Bristol)
------------------------------------------------------------------------- */

#include "fix_brownian_omp.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "random_mars.h"

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBrownianOMP::FixBrownianOMP(LAMMPS *lmp, int narg, char **arg) : FixBrownianBaseOMP(lmp, narg, arg)
{
    if (dipole_flag || gamma_t_eigen_flag || gamma_r_eigen_flag || gamma_r_flag) {
        error->all(FLERR, "Illegal fix brownian command.");
    }
    if (!gamma_t_flag) { error->all(FLERR, "Illegal fix brownian command."); }
}

/* ---------------------------------------------------------------------- */

void FixBrownianOMP::init()
{
    FixBrownianBaseOMP::init();
    g1 /= gamma_t;
    g2 /= gamma_t;
}

/* ---------------------------------------------------------------------- */

void FixBrownianOMP::initial_integrate(int /*vflag */)
{
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(noise_flag)
#endif
    { // begin parallel region

        int ifrom, ito, tid;

        loop_setup_thr(ifrom, ito, tid, atom->nlocal, nthreads);
        ThrData *thr = fix->get_thr(tid);
        // thr->timer(Timer::START);
        // ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

        // generate a random number generator instance for
        // all threads != 0. make sure we use unique seeds.
        if ((tid > 0) && (random_thr[tid] == nullptr))
            random_thr[tid] = new RanMars(Fix::lmp, seed + comm->me
                                                    + comm->nprocs * tid);

        if (domain->dimension == 2) {
            if (!noise_flag) {
                initial_integrate_templated<0, 0, 1>(ifrom, ito, thr);
            } else if (gaussian_noise_flag) {
                initial_integrate_templated<0, 1, 1>(ifrom, ito, thr);
            } else {
                initial_integrate_templated<1, 0, 1>(ifrom, ito, thr);
            }
        } else {
            if (!noise_flag) {
                initial_integrate_templated<0, 0, 0>(ifrom, ito, thr);
            } else if (gaussian_noise_flag) {
                initial_integrate_templated<0, 1, 0>(ifrom, ito, thr);
            } else {
                initial_integrate_templated<1, 0, 0>(ifrom, ito, thr);
            }
        }
    } // end parallel region
    return;
}

/* ---------------------------------------------------------------------- */

template <int Tp_UNIFORM, int Tp_GAUSS, int Tp_2D>
void FixBrownianOMP::initial_integrate_templated(int ifrom, int ito, ThrData * const thr)
{
    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    int *mask = atom->mask;
    // int nlocal = atom->nlocal;

    // if (igroup == atom->firstgroup) nlocal = atom->nfirst;

    double dx, dy, dz;

    RanMars &rng = *random_thr[thr->get_tid()];

    for (int i = ifrom; i < ito; i++) {
        if (mask[i] & groupbit) {
            if (Tp_2D) {
                dz = 0;
                if (Tp_UNIFORM) {
                    dx = dt * (g1 * f[i][0] + g2 * (rng.uniform() - 0.5));
                    dy = dt * (g1 * f[i][1] + g2 * (rng.uniform() - 0.5));
                } else if (Tp_GAUSS) {
                    dx = dt * (g1 * f[i][0] + g2 * rng.gaussian());
                    dy = dt * (g1 * f[i][1] + g2 * rng.gaussian());
                } else {
                    dx = dt * g1 * f[i][0];
                    dy = dt * g1 * f[i][1];
                }
            } else {
                if (Tp_UNIFORM) {
                    dx = dt * (g1 * f[i][0] + g2 * (rng.uniform() - 0.5));
                    dy = dt * (g1 * f[i][1] + g2 * (rng.uniform() - 0.5));
                    dz = dt * (g1 * f[i][2] + g2 * (rng.uniform() - 0.5));
                } else if (Tp_GAUSS) {
                    dx = dt * (g1 * f[i][0] + g2 * rng.gaussian());
                    dy = dt * (g1 * f[i][1] + g2 * rng.gaussian());
                    dz = dt * (g1 * f[i][2] + g2 * rng.gaussian());
                } else {
                    dx = dt * g1 * f[i][0];
                    dy = dt * g1 * f[i][1];
                    dz = dt * g1 * f[i][2];
                }
            }

            x[i][0] += dx;
            v[i][0] = dx / dt;

            x[i][1] += dy;
            v[i][1] = dy / dt;

            x[i][2] += dz;
            v[i][2] = dz / dt;
        }
    }

    return;
}
