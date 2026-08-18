#pragma once

// Bridge between btree_chromo data structures and the fused CG minimizer.
// Include this from btree_chromo code to convert atom_array/bond_array/angle_array
// into FusedMinSystem for direct GPU minimization.

#include "fused_cg_minimize.h"
#include <atom_array.hpp>
#include <bond_array.hpp>
#include <angle_array.hpp>

// Populate a FusedMinSystem from btree_chromo arrays.
// Allocates sys fields — caller must call fused_free_system() when done.
// atom_ids are converted from 1-indexed LAMMPS IDs to 0-indexed.
inline void fused_sys_from_btree(
    atom_array &atoms,
    bond_array &bonds,
    angle_array &angles,
    FusedMinSystem *sys)
{
    int N = atoms.get_N();
    int nb = bonds.get_N();
    int na = angles.get_N();

    sys->N       = N;
    sys->nbonds  = nb;
    sys->nangles = na;

    sys->x    = (double*)malloc(N * 3 * sizeof(double));
    sys->xref = (double*)malloc(N * 3 * sizeof(double));
    sys->type = (int*)malloc(N * sizeof(int));

    for (int i = 0; i < N; i++) {
        atom a = atoms.get_atom(i);
        sys->type[i]    = a.type;
        sys->x[3*i]     = a.r.x;
        sys->x[3*i + 1] = a.r.y;
        sys->x[3*i + 2] = a.r.z;
        // Reference positions for spring/self (ribosomes)
        sys->xref[3*i]     = a.r.x;
        sys->xref[3*i + 1] = a.r.y;
        sys->xref[3*i + 2] = a.r.z;
    }

    if (nb > 0) {
        sys->bond_type = (int*)malloc(nb * sizeof(int));
        sys->bond_i    = (int*)malloc(nb * sizeof(int));
        sys->bond_j    = (int*)malloc(nb * sizeof(int));
        for (int b = 0; b < nb; b++) {
            bond bd = bonds.get_bond(b);
            sys->bond_type[b] = bd.type;
            sys->bond_i[b]    = bd.i - 1; // LAMMPS 1-indexed → 0-indexed
            sys->bond_j[b]    = bd.j - 1;
        }
    }

    if (na > 0) {
        sys->angle_type = (int*)malloc(na * sizeof(int));
        sys->angle_i    = (int*)malloc(na * sizeof(int));
        sys->angle_j    = (int*)malloc(na * sizeof(int));
        sys->angle_k    = (int*)malloc(na * sizeof(int));
        for (int a = 0; a < na; a++) {
            angle ag = angles.get_angle(a);
            sys->angle_type[a] = ag.type;
            sys->angle_i[a]    = ag.i - 1;
            sys->angle_j[a]    = ag.j - 1;
            sys->angle_k[a]    = ag.k - 1;
        }
    }
}

// Write minimized positions back to btree_chromo atom_array.
inline void fused_sys_to_btree(
    const FusedMinSystem *sys,
    atom_array &atoms)
{
    int N = atoms.get_N();
    for (int i = 0; i < N; i++) {
        atom a = atoms.get_atom(i);
        a.r.x = sys->x[3*i];
        a.r.y = sys->x[3*i + 1];
        a.r.z = sys->x[3*i + 2];
        atoms.set_atom(i, a);
    }
}
