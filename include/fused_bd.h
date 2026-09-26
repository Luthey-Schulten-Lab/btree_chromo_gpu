// Author(s): Ron Acda (using an iterative LLM-guided workflow, https://github.com/quarkron/iterative-hillclimber/tree/main)
#pragma once
// fused Brownian dynamics for the soft-harmonic DNA run (replaces the LAMMPS Kokkos "run" of run_soft_harmonic).
//
// The model integrated by LAMMPS there (all coefficients are taken from the live LAMMPS instance by the caller):
//   pair   soft:     E = A (1 + cos(pi r / rc)),  r < rc          (per type pair A, rc; A = 0 or rc <= 0: no interaction)
//   bond   harmonic: E = K (r - r0)^2                              (per bond type)
//   angle  cosine:   E = K (1 + cos theta)                         (per angle type)
//   special pairs (LAMMPS special list, lj weights 0) excluded from the pair interaction
//   fix brownian (uniform noise) on the mobile atoms: x += dt * (g1 * f + g2 * (u - 0.5)), u ~ U(0,1) per component
// Neighbour lists with a skin, displacement check every `check_every` steps, rebuild when an atom moved more than skin/2.
// Noise: counter-based Philox4x32-10 keyed by (seed, atom, step): statistically equivalent to LAMMPS's per-atom RanMars /
// Kokkos pool uniforms, not the same stream .

#include <cstddef>

#define FUSED_BD_MAXT 10    // atom types 1..9 (index 0 unused)
#define FUSED_BD_MAXB 8     // bond / angle types 1..7

struct FusedBDSystem {
    int N;
    double *x;              // [3N] row-major, in: start positions, out: end positions
    const int *type;        // [N] 1-based atom types
    const unsigned char *mobile;  // [N] 1 = integrated (member of the Brownian fix's group)
    int nbonds;  const int *bond_type, *bond_i, *bond_j;            // 0-based atom indices
    int nangles; const int *angle_type, *angle_i, *angle_j, *angle_k;
    const int *excl_off, *excl_list;   // CSR [N+1] / [excl_off[N]]: excluded partners per atom, sorted ascending
    double pair_A[FUSED_BD_MAXT][FUSED_BD_MAXT], pair_rc[FUSED_BD_MAXT][FUSED_BD_MAXT];
    double bond_K[FUSED_BD_MAXB], bond_r0[FUSED_BD_MAXB];
    double angle_K[FUSED_BD_MAXB];
    double dt, g1, g2;
    double skin;
    int check_every;
    unsigned long long seed;
    double box_lo[3], box_hi[3];
};

struct FusedBDStats {
    int builds;
    double max_disp_at_build;     // largest displacement since the previous build seen at any rebuild check (A)
    double kernel_ms, build_ms, total_ms;
};

// Runs nsteps steps on GPU 0; returns 0 on success.
int fused_bd_run(FusedBDSystem *s, long nsteps, FusedBDStats *stats);
