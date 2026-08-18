#pragma once

#include <cstddef>

constexpr int NUM_TYPES = 8;  // LAMMPS atom types 1-8

enum StopCode {
    STOP_MAXITER        = 0,
    STOP_ETOL           = 1,
    STOP_FTOL           = 2,
    STOP_MAXEVAL        = 3,
    STOP_NEEDS_REBUILD  = 4,
    STOP_ENERGY_TARGET  = 5,
    STOP_RUNNING        = -1
};

enum FusedPairStyle  { FUSED_PAIR_HARD = 0, FUSED_PAIR_SOFT = 1, FUSED_PAIR_TOPO = 2 };
enum FusedBondStyle  { FUSED_BOND_HARMONIC = 0, FUSED_BOND_FENE = 1 };

// CG convergence budget profile.
// FULL          = original tight defaults (etol=1e-5, ftol=1e-7, maxiter=40000,
//                 dmax=0.1) used for directive-level minimizes (first call,
//                 relax_progressive, equilibration templates). Converges to a
//                 true stationary point.
// REFRESH       = loose, iteration-capped budget used for the FIRST in-loop
//                 minimize of each run_loops iteration (or the post-update_loops
//                 minimize). Handles the high-energy basin after fresh SMC bond
//                 creation. maxiter=8000.
// REFRESH_SHORT = even tighter cap for IN-LOOP refreshes after iter 1 of
//                 simulator_run_loops, where the system is already
//                 well-equilibrated by the preceding minimize+BD step.
//                 In job 1514, hard_harmonic averaged 3000-7000 iters/call
//                 inside the hump; capping at 1500 cuts the in-loop budget
//                 ~3-5x without hitting maxiter on typical iters.
enum FusedMinProfile {
    FUSED_MIN_FULL          = 0,
    FUSED_MIN_REFRESH       = 1,
    FUSED_MIN_REFRESH_SHORT = 2,
};

struct FusedMinParams {
    // Per-type-pair mode: 0=off, 1=WCA, 2=soft
    int pair_mode[9][9];

    // Unified squared cutoff (WCA rc or soft rc, whichever is active)
    double cutoff_sq[9][9];

    // WCA pair parameters (used when pair_mode == 1)
    double sigma_sq[9][9];
    double sigma_6[9][9];
    double epsilon[9][9];

    // Soft pair parameters (used when pair_mode == 2)
    double soft_A[9][9];       // prefactor A
    double soft_rc[9][9];      // cutoff rc

    // Bond coeffs, indexed by LAMMPS bond type (index 0 unused).
    // type 1 = backbone, type 2 = SMC loop, type 3 = replication fork
    // (protein_science: bond_coeff 3 = k_stretch, r0=0). Sized [5] for headroom.
    int    bond_is_fene[5];    // 1 if FENE, 0 if harmonic
    double bond_K[5];
    double bond_r0[5];         // harmonic equilibrium (only when !fene)

    // FENE-specific (only used when bond_is_fene[bt] == 1)
    double fene_R0_sq[5];      // max extension squared
    double fene_eps[5];        // LJ epsilon within FENE
    double fene_sigma_6[5];    // LJ sigma^6
    double fene_inner_sq[5];   // (2^(1/6)*sigma)^2 WCA inner cutoff

    // Angle: types 1-4 (index 0 unused)
    double angle_K[5];

    // Spring/self for ribosomes (type 2)
    double spring_K_ribo;

    // CG convergence
    double etol;
    double ftol;
    int    maxiter;
    int    maxeval;
    double dmax;

    // Neighbor list skin for rebuild detection
    double skin;

    // Adaptive early-stop: exit CG when energy drops below this threshold (0 = disabled)
    double energy_target;

    // Dynamic early-stop: exit when relative energy improvement per iteration
    // drops below this fraction (0 = disabled). Scale-independent.
    // E.g. 1e-4 means stop when (E_prev - E_curr) / E_curr < 0.01% per iteration.
    double early_stop_rtol;
};

struct FusedMinSystem {
    int    N;
    int    nbonds;
    int    nangles;

    double *x;          // [N*3] AOS: x y z x y z ...
    int    *type;       // [N] atom types, 1-indexed
    double *xref;       // [N*3] reference positions for spring/self

    int    *bond_type;  // [nbonds]
    int    *bond_i;     // [nbonds] 0-indexed
    int    *bond_j;     // [nbonds] 0-indexed

    int    *angle_type; // [nangles]
    int    *angle_i;    // [nangles]
    int    *angle_j;    // [nangles] vertex
    int    *angle_k;    // [nangles]
};

struct FusedNeighList {
    int    N;
    int    max_neighs;
    int    *numneigh;   // [N]
    int    *neighlist;  // [N * max_neighs]
};

struct FusedMinResult {
    StopCode stop;
    int      niter;
    int      neval;
    double   energy_initial;
    double   energy_final;
    double   fnorm_final;
    double   elapsed_ms;
};

struct FusedDeviceCache {
    bool initialized;
    int N_cap;              // allocated capacity (atoms)
    int N_cur;              // atom count from last call
    int cached_nbonds;
    int cached_nangles;
    size_t topo_hash;       // fingerprint of bond/angle arrays for cache invalidation
    int gpu_id;
    int max_blocks;
    int block_size;

    double *d_x, *d_f, *d_h, *d_g, *d_x0, *d_x_build, *d_xref;
    int    *d_type;

    int *d_atom_bond_off, *d_atom_bond_data;
    int *d_atom_angle_off, *d_atom_angle_data;
    int bond_data_cap, angle_data_cap;

    int *d_excl_off, *d_excl_list;
    int excl_cap;

    int max_neighs;
    int *d_numneigh, *d_neighlist;

    int *d_cell_of_atom, *d_sorted_atoms, *d_cell_start, *d_cell_count;
    int cell_cap;

    double *d_block_sums_a, *d_block_sums_b, *d_scalars;
    int *d_stop, *d_niter, *d_neval;
    double *d_energy_init, *d_energy_final, *d_fnorm, *d_gg;
};

void fused_cache_init(FusedDeviceCache *cache);
void fused_cache_destroy(FusedDeviceCache *cache);

void fused_min_default_params(FusedMinParams *p);
void fused_min_init_params(FusedMinParams *p, int pair_style, int bond_style,
                           int profile = FUSED_MIN_FULL);

void fused_build_neighbor_list(
    const double *x, const int *type, int N,
    const FusedMinSystem *sys,
    double skin, double max_cutoff,
    const double box_lo[3], const double box_hi[3],
    FusedNeighList *nlist
);

FusedMinResult fused_cg_minimize(
    FusedMinSystem *sys,
    const FusedMinParams *params,
    const double box_lo[3], const double box_hi[3],
    int gpu_id,
    FusedDeviceCache *cache
);

int fused_read_lammps_data(
    const char *filename,
    FusedMinSystem *sys,
    double box_lo[3],
    double box_hi[3]
);

void fused_free_system(FusedMinSystem *sys);
void fused_free_neighlist(FusedNeighList *nlist);
