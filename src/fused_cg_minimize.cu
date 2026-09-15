#include "fused_cg_minimize.h"

#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>
#include <chrono>
#include <thrust/device_ptr.h>
#include <thrust/scan.h>

namespace cg = cooperative_groups;

#define CUDA_CHECK(call) do { \
    cudaError_t err = (call); \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error %s:%d: %s\n", \
                __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(1); \
    } \
} while(0)

static inline bool fused_profile_enabled() {
    const char *e = std::getenv("FUSED_PROFILE");
    return e && e[0] != '\0' && e[0] != '0';
}

static inline double fused_ms_since(std::chrono::steady_clock::time_point t0) {
    using clock = std::chrono::steady_clock;
    return std::chrono::duration<double, std::milli>(clock::now() - t0).count();
}

static constexpr double SIXRT2 = 1.1224620483093730; // 2^(1/6)
static constexpr double FUSED_PI = 3.14159265358979323846;

// ============================================================================
//  Constant memory for force-field parameters
// ============================================================================

__constant__ FusedMinParams d_params;

// ============================================================================
//  Device data passed to the persistent kernel
// ============================================================================

struct KernelData {
    int N;
    double *x, *f, *h, *g, *x0, *x_build;
    int *type;
    double *xref;

    // Per-atom bond topology (CSR): entry = (bond_type, other_atom)
    int *atom_bond_off;   // [N+1]
    int *atom_bond_data;  // [2 * total_bond_entries]

    // Per-atom angle topology (CSR): entry = (angle_type, ai, aj, ak)
    int *atom_angle_off;  // [N+1]
    int *atom_angle_data; // [4 * total_angle_entries]

    int *numneigh, *neighlist;
    int max_neighs;
    int neighlist_stride;

    double *block_sums_a, *block_sums_b;
    double *scalars;

    double skin;
    int start_iter, start_neval;
    int resume;

    int *d_stop;
    int *d_niter, *d_neval;
    double *d_energy_init, *d_energy_final, *d_fnorm, *d_gg;
};

// ============================================================================
//  Device: block and grid reductions
// ============================================================================

__device__ __forceinline__
double block_reduce_sum(double val, double *sdata) {
    sdata[threadIdx.x] = val;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) sdata[threadIdx.x] += sdata[threadIdx.x + s];
        __syncthreads();
    }
    return sdata[0];
}

__device__ __forceinline__
double block_reduce_max(double val, double *sdata) {
    sdata[threadIdx.x] = val;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s)
            sdata[threadIdx.x] = fmax(sdata[threadIdx.x], sdata[threadIdx.x + s]);
        __syncthreads();
    }
    return sdata[0];
}

__device__
double grid_reduce_sum(cg::grid_group &grid, double val,
                       double *sdata, double *bsums, double *scalars, int slot) {
    double bs = block_reduce_sum(val, sdata);
    if (threadIdx.x == 0) bsums[blockIdx.x] = bs;
    grid.sync();
    if (blockIdx.x == 0) {
        double s = 0.0;
        for (int b = threadIdx.x; b < (int)gridDim.x; b += blockDim.x)
            s += bsums[b];
        s = block_reduce_sum(s, sdata);
        if (threadIdx.x == 0) scalars[slot] = s;
    }
    grid.sync();
    return scalars[slot];
}

__device__
double grid_reduce_max(cg::grid_group &grid, double val,
                       double *sdata, double *bsums, double *scalars, int slot) {
    double bm = block_reduce_max(val, sdata);
    if (threadIdx.x == 0) bsums[blockIdx.x] = bm;
    grid.sync();
    if (blockIdx.x == 0) {
        double m = -1e300;
        for (int b = threadIdx.x; b < (int)gridDim.x; b += blockDim.x)
            m = fmax(m, bsums[b]);
        m = block_reduce_max(m, sdata);
        if (threadIdx.x == 0) scalars[slot] = m;
    }
    grid.sync();
    return scalars[slot];
}

// Fused: reduce two sums in parallel using bsums_a (slot 0) and bsums_b (slot 1).
// Saves 2 grid.sync() vs calling grid_reduce_sum twice.
__device__
void grid_reduce_sum_2(cg::grid_group &grid,
                       double val_a, double val_b,
                       double *sdata, double *bsums_a, double *bsums_b,
                       double *scalars,
                       double &out_a, double &out_b) {
    double bs_a = block_reduce_sum(val_a, sdata);
    if (threadIdx.x == 0) bsums_a[blockIdx.x] = bs_a;
    __syncthreads();
    double bs_b = block_reduce_sum(val_b, sdata);
    if (threadIdx.x == 0) bsums_b[blockIdx.x] = bs_b;
    grid.sync();
    if (blockIdx.x == 0) {
        double sa = 0.0;
        for (int b = threadIdx.x; b < (int)gridDim.x; b += blockDim.x)
            sa += bsums_a[b];
        sa = block_reduce_sum(sa, sdata);
        if (threadIdx.x == 0) scalars[0] = sa;
        __syncthreads();
        double sb = 0.0;
        for (int b = threadIdx.x; b < (int)gridDim.x; b += blockDim.x)
            sb += bsums_b[b];
        sb = block_reduce_sum(sb, sdata);
        if (threadIdx.x == 0) scalars[1] = sb;
    }
    grid.sync();
    out_a = scalars[0];
    out_b = scalars[1];
}

// Fused: reduce a sum and a max in parallel.
// Saves 2 grid.sync() vs calling grid_reduce_sum + grid_reduce_max.
__device__
void grid_reduce_sum_max(cg::grid_group &grid,
                         double sum_val, double max_val,
                         double *sdata, double *bsums_a, double *bsums_b,
                         double *scalars,
                         double &out_sum, double &out_max) {
    double bs = block_reduce_sum(sum_val, sdata);
    if (threadIdx.x == 0) bsums_a[blockIdx.x] = bs;
    __syncthreads();
    double bm = block_reduce_max(max_val, sdata);
    if (threadIdx.x == 0) bsums_b[blockIdx.x] = bm;
    grid.sync();
    if (blockIdx.x == 0) {
        double s = 0.0;
        for (int b = threadIdx.x; b < (int)gridDim.x; b += blockDim.x)
            s += bsums_a[b];
        s = block_reduce_sum(s, sdata);
        if (threadIdx.x == 0) scalars[0] = s;
        __syncthreads();
        double m = -1e300;
        for (int b = threadIdx.x; b < (int)gridDim.x; b += blockDim.x)
            m = fmax(m, bsums_b[b]);
        m = block_reduce_max(m, sdata);
        if (threadIdx.x == 0) scalars[1] = m;
    }
    grid.sync();
    out_sum = scalars[0];
    out_max = scalars[1];
}

// ============================================================================
//  Device: atom-centric force computation
// ============================================================================

__device__
double compute_atom_forces(KernelData &kd, int a,
                           double *out_fx, double *out_fy, double *out_fz) {
    int ta = __ldg(&kd.type[a]);

    if (ta == 1) {
        *out_fx = *out_fy = *out_fz = 0.0;
        return 0.0;
    }

    double pe = 0.0;
    double fx = 0.0, fy = 0.0, fz = 0.0;
    double xa = kd.x[3*a], ya = kd.x[3*a+1], za = kd.x[3*a+2];

    // ---- Pair forces (WCA or soft, Newton OFF / full neighbor list) ----
    int nn = __ldg(&kd.numneigh[a]);
    for (int jj = 0; jj < nn; jj++) {
        int j = __ldg(&kd.neighlist[(size_t)jj * kd.neighlist_stride + a]);
        int tj = __ldg(&kd.type[j]);
        double dx = xa - __ldg(&kd.x[3*j]);
        double dy = ya - __ldg(&kd.x[3*j+1]);
        double dz = za - __ldg(&kd.x[3*j+2]);
        double r2 = dx*dx + dy*dy + dz*dz;

        int pm = d_params.pair_mode[ta][tj];
        if (pm == 0) continue;
        double csq = d_params.cutoff_sq[ta][tj];
        if (r2 >= csq || r2 < 1e-20) continue;

        if (pm == 1) { // WCA (shifted LJ)
            double eps = d_params.epsilon[ta][tj];
            if (eps == 0.0) continue;
            double r2inv = 1.0 / r2;
            double r6inv = d_params.sigma_6[ta][tj] * r2inv * r2inv * r2inv;
            double r12inv = r6inv * r6inv;
            double fmag = 48.0 * eps * r2inv * (r12inv - 0.5 * r6inv);
            fx += fmag * dx; fy += fmag * dy; fz += fmag * dz;
            pe += 0.5 * (4.0 * eps * (r12inv - r6inv) + eps);
        } else { // Soft: E = A*(1+cos(pi*r/rc))
            double A  = d_params.soft_A[ta][tj];
            double rc = d_params.soft_rc[ta][tj];
            double r  = sqrt(r2);
            double arg = FUSED_PI * r / rc;
            double fmag = A * FUSED_PI / (rc * r) * sin(arg);
            fx += fmag * dx; fy += fmag * dy; fz += fmag * dz;
            pe += 0.5 * A * (1.0 + cos(arg));
        }
    }

    // ---- Bond forces (harmonic or FENE, atom-centric) ----
    int b_start = __ldg(&kd.atom_bond_off[a]);
    int b_end   = __ldg(&kd.atom_bond_off[a + 1]);
    for (int bi = b_start; bi < b_end; bi++) {
        int bt    = __ldg(&kd.atom_bond_data[2*bi]);
        int other = __ldg(&kd.atom_bond_data[2*bi + 1]);

        double dx = xa - __ldg(&kd.x[3*other]);
        double dy = ya - __ldg(&kd.x[3*other+1]);
        double dz = za - __ldg(&kd.x[3*other+2]);
        double rsq = dx*dx + dy*dy + dz*dz;
        if (rsq < 1e-20) continue;

        if (d_params.bond_is_fene[bt]) {
            double R0sq = d_params.fene_R0_sq[bt];
            double rr = 1.0 - rsq / R0sq;
            if (rr < 0.01) rr = 0.01;

            double fscale = -d_params.bond_K[bt] / rr;
            double e_wca = 0.0;
            if (rsq < d_params.fene_inner_sq[bt]) {
                double r2inv = 1.0 / rsq;
                double r6inv = d_params.fene_sigma_6[bt] * r2inv * r2inv * r2inv;
                fscale += 48.0 * d_params.fene_eps[bt] * r2inv * r6inv * (r6inv - 0.5);
                e_wca = 4.0 * d_params.fene_eps[bt] * r6inv * (r6inv - 1.0)
                        + d_params.fene_eps[bt];
            }
            fx += fscale * dx; fy += fscale * dy; fz += fscale * dz;
            pe += 0.5 * (-0.5 * d_params.bond_K[bt] * R0sq * log(rr) + e_wca);
        } else {
            double r  = sqrt(rsq);
            double dr = r - d_params.bond_r0[bt];
            double fscale = -2.0 * d_params.bond_K[bt] * dr / r;
            fx += fscale * dx; fy += fscale * dy; fz += fscale * dz;
            pe += 0.5 * d_params.bond_K[bt] * dr * dr;
        }
    }

    // ---- Angle forces (cosine, atom-centric) ----
    int ag_start = __ldg(&kd.atom_angle_off[a]);
    int ag_end   = __ldg(&kd.atom_angle_off[a + 1]);
    for (int ai = ag_start; ai < ag_end; ai++) {
        int atype = __ldg(&kd.atom_angle_data[4*ai]);
        int a_i   = __ldg(&kd.atom_angle_data[4*ai + 1]);
        int a_j   = __ldg(&kd.atom_angle_data[4*ai + 2]);
        int a_k   = __ldg(&kd.atom_angle_data[4*ai + 3]);

        double d1x = __ldg(&kd.x[3*a_i])   - __ldg(&kd.x[3*a_j]);
        double d1y = __ldg(&kd.x[3*a_i+1]) - __ldg(&kd.x[3*a_j+1]);
        double d1z = __ldg(&kd.x[3*a_i+2]) - __ldg(&kd.x[3*a_j+2]);
        double d2x = __ldg(&kd.x[3*a_k])   - __ldg(&kd.x[3*a_j]);
        double d2y = __ldg(&kd.x[3*a_k+1]) - __ldg(&kd.x[3*a_j+1]);
        double d2z = __ldg(&kd.x[3*a_k+2]) - __ldg(&kd.x[3*a_j+2]);

        double rsq1 = d1x*d1x + d1y*d1y + d1z*d1z;
        double rsq2 = d2x*d2x + d2y*d2y + d2z*d2z;
        double r1 = sqrt(rsq1), r2 = sqrt(rsq2);
        if (r1 < 1e-20 || r2 < 1e-20) continue;

        double c = (d1x*d2x + d1y*d2y + d1z*d2z) / (r1 * r2);
        c = fmax(-1.0, fmin(1.0, c));

        double K = d_params.angle_K[atype];
        pe += K * (1.0 + c) / 3.0;

        double inv_r1r2 = 1.0 / (r1 * r2);
        double c_rsq1 = c / rsq1;
        double c_rsq2 = c / rsq2;

        double fi_x = -K * (d2x * inv_r1r2 - c_rsq1 * d1x);
        double fi_y = -K * (d2y * inv_r1r2 - c_rsq1 * d1y);
        double fi_z = -K * (d2z * inv_r1r2 - c_rsq1 * d1z);
        double fk_x = -K * (d1x * inv_r1r2 - c_rsq2 * d2x);
        double fk_y = -K * (d1y * inv_r1r2 - c_rsq2 * d2y);
        double fk_z = -K * (d1z * inv_r1r2 - c_rsq2 * d2z);

        if (a == a_i)      { fx += fi_x; fy += fi_y; fz += fi_z; }
        else if (a == a_k) { fx += fk_x; fy += fk_y; fz += fk_z; }
        else               { fx -= (fi_x+fk_x); fy -= (fi_y+fk_y); fz -= (fi_z+fk_z); }
    }

    // ---- Fixes ----
    if (ta == 2) { // ribosome: spring/self
        double dx = xa - __ldg(&kd.xref[3*a]);
        double dy = ya - __ldg(&kd.xref[3*a+1]);
        double dz = za - __ldg(&kd.xref[3*a+2]);
        double Ks = d_params.spring_K_ribo;
        fx -= Ks * dx;
        fy -= Ks * dy;
        fz -= Ks * dz;
        pe += 0.5 * Ks * (dx*dx + dy*dy + dz*dz);
    }

    *out_fx = fx; *out_fy = fy; *out_fz = fz;
    return pe;
}

// Full energy+force evaluation: each thread handles its assigned atoms
__device__
double full_energy_force(KernelData &kd, cg::grid_group &grid,
                         double *sdata, int tid, int nthreads) {
    double pe_local = 0.0;
    for (int a = tid; a < kd.N; a += nthreads) {
        double fx, fy, fz;
        pe_local += compute_atom_forces(kd, a, &fx, &fy, &fz);
        kd.f[3*a]   = fx;
        kd.f[3*a+1] = fy;
        kd.f[3*a+2] = fz;
    }
    grid.sync();
    return grid_reduce_sum(grid, pe_local, sdata, kd.block_sums_a, kd.scalars, 0);
}

// ============================================================================
//  The persistent CG kernel — matches LAMMPS min_cg_kokkos + linemin_quadratic
// ============================================================================

// Constants from LAMMPS min_linesearch_kokkos.cpp
static constexpr double LM_ALPHA_MAX      = 1.0;
static constexpr double LM_ALPHA_REDUCE   = 0.5;
static constexpr double LM_BACKTRACK_SLOPE = 0.4;
static constexpr double LM_QUADRATIC_TOL  = 0.1;
static constexpr double LM_EMACH          = 1.0e-8;
static constexpr double LM_EPS_QUAD       = 1.0e-28;
static constexpr double LM_EPS_ENERGY     = 1.0e-8;

__global__ void fused_cg_kernel(KernelData kd) {
    extern __shared__ double sdata[];
    auto grid = cg::this_grid();
    int tid = grid.thread_rank();
    int nthreads = grid.num_threads();
    int N3 = kd.N * 3;
    double inv_N = 1.0 / (double)kd.N;

    int neval = kd.start_neval;
    int niter = 0;
    int stop = STOP_RUNNING;

    // nlimit for periodic CG restart (LAMMPS: min(MAXSMALLINT, ndoftotal))
    int nlimit = (N3 < 2000000000) ? N3 : 2000000000;

    // ---- Initial energy & forces ----
    double ecurrent = full_energy_force(kd, grid, sdata, tid, nthreads);
    neval++;

    if (tid == 0) *kd.d_energy_init = ecurrent;

    // Initial force norm (LAMMPS: TWO norm => fnorm_sqr = Σf²)
    double ff = 0.0;
    for (int i = tid; i < N3; i += nthreads) ff += kd.f[i] * kd.f[i];
    double gnorm2 = grid_reduce_sum(grid, ff, sdata, kd.block_sums_a, kd.scalars, 0);

    if (gnorm2 < d_params.ftol * d_params.ftol) {
        stop = STOP_FTOL;
    }

    double gg;
    if (!kd.resume) {
        // First launch: init CG with steepest descent
        for (int i = tid; i < N3; i += nthreads) {
            kd.g[i] = kd.f[i];
            kd.h[i] = kd.f[i];
        }
        gg = gnorm2;
    } else {
        // Resume after neighbor rebuild: keep g[] and h[],
        // recover gg from saved scalar
        gg = *kd.d_gg;
    }

    // Always update x_build after a (re)build
    for (int i = tid; i < N3; i += nthreads) kd.x_build[i] = kd.x[i];
    grid.sync();

    // ---- Main CG loop (matches min_cg_kokkos::iterate) ----
    for (niter = 0; niter < d_params.maxiter && stop == STOP_RUNNING; niter++) {
        double eprevious = ecurrent;

        // ================================================================
        //  linemin_quadratic  (matches min_linesearch_kokkos.cpp)
        // ================================================================

        // Fused fdothall (sum) + hmaxall (max) in one pass over h[]/f[]
        double fdoth_l = 0.0, hmax_l = 0.0;
        for (int i = tid; i < N3; i += nthreads) {
            fdoth_l += kd.f[i] * kd.h[i];
            hmax_l = fmax(hmax_l, fabs(kd.h[i]));
        }
        double fdothall, hmaxall;
        grid_reduce_sum_max(grid, fdoth_l, hmax_l, sdata,
                            kd.block_sums_a, kd.block_sums_b, kd.scalars,
                            fdothall, hmaxall);
        fdothall *= inv_N;

        if (fdothall <= 0.0) { stop = STOP_ETOL; break; }
        if (hmaxall == 0.0) { stop = STOP_FTOL; break; }

        double alphamax = fmin(LM_ALPHA_MAX, d_params.dmax / hmaxall);

        // Save x0
        for (int i = tid; i < N3; i += nthreads) kd.x0[i] = kd.x[i];
        grid.sync();

        // Initialize line search state
        double alpha = alphamax;
        double fhprev = fdothall;
        double engprev = eprevious;
        double alphaprev = 0.0;
        int ls_fail = 0;

        // Backtracking loop (matches LAMMPS while(true) in linemin_quadratic)
        for (int ls = 0; ls < 50 && !ls_fail; ls++) {
            // alpha_step: reset to x0, then step to x0 + alpha*h
            for (int i = tid; i < N3; i += nthreads)
                kd.x[i] = kd.x0[i] + alpha * kd.h[i];
            grid.sync();
            ecurrent = full_energy_force(kd, grid, sdata, tid, nthreads);
            neval++;

            // Compute fh = f · h (normalized)
            double fh_l = 0.0;
            for (int i = tid; i < N3; i += nthreads) fh_l += kd.f[i] * kd.h[i];
            double fh = grid_reduce_sum(grid, fh_l, sdata, kd.block_sums_a, kd.scalars, 0);
            fh *= inv_N;

            double delfh = fh - fhprev;

            // ZEROQUAD check
            if (fabs(fh) < LM_EPS_QUAD || fabs(delfh) < LM_EPS_QUAD) {
                for (int i = tid; i < N3; i += nthreads) kd.x[i] = kd.x0[i];
                grid.sync();
                ecurrent = full_energy_force(kd, grid, sdata, tid, nthreads);
                neval++;
                ls_fail = 1;
                break;
            }

            // Quadratic projection via secant method
            double relerr = fabs(1.0 - (0.5 * (alpha - alphaprev) * (fh + fhprev)
                                        + ecurrent) / engprev);
            double alpha0 = alpha - (alpha - alphaprev) * fh / delfh;

            if (relerr <= LM_QUADRATIC_TOL && alpha0 > 0.0 && alpha0 < alphamax) {
                for (int i = tid; i < N3; i += nthreads)
                    kd.x[i] = kd.x0[i] + alpha0 * kd.h[i];
                grid.sync();
                ecurrent = full_energy_force(kd, grid, sdata, tid, nthreads);
                neval++;
                if (ecurrent - eprevious < LM_EMACH)
                    break; // accept secant step
            }

            // Armijo sufficient decrease check
            double de_ideal = -LM_BACKTRACK_SLOPE * alpha * fdothall;
            double de = ecurrent - eprevious;
            if (de <= de_ideal)
                break; // accept

            // Save state and reduce alpha
            fhprev = fh;
            engprev = ecurrent;
            alphaprev = alpha;
            alpha *= LM_ALPHA_REDUCE;

            // ZEROALPHA check
            if (alpha <= 0.0 || de_ideal >= -LM_EMACH) {
                for (int i = tid; i < N3; i += nthreads) kd.x[i] = kd.x0[i];
                grid.sync();
                ecurrent = full_energy_force(kd, grid, sdata, tid, nthreads);
                neval++;
                ls_fail = 1;
                break;
            }
        }

        // LAMMPS aborts on line search failure
        if (ls_fail) { stop = STOP_ETOL; break; }

        // ================================================================
        //  Post-linesearch convergence checks (matches min_cg_kokkos)
        // ================================================================

        if (neval >= d_params.maxeval) { stop = STOP_MAXEVAL; break; }

        // Energy tolerance (matches LAMMPS EPS_ENERGY = 1e-8)
        double ediff = fabs(ecurrent - eprevious);
        double eref  = 0.5 * (fabs(ecurrent) + fabs(eprevious) + LM_EPS_ENERGY);
        if (ediff < d_params.etol * eref) { stop = STOP_ETOL; break; }

        // Adaptive early-stop: exit when energy drops below target threshold
        if (d_params.energy_target > 0.0 && ecurrent < d_params.energy_target) {
            stop = STOP_ENERGY_TARGET; break;
        }

        // Dynamic early-stop: stop when relative improvement per iteration is
        // negligible. Scale-independent — works for any starting energy.
        if (d_params.early_stop_rtol > 0.0 && ecurrent > 0.0) {
            double rel_drop = (eprevious - ecurrent) / ecurrent;
            if (rel_drop < d_params.early_stop_rtol) {
                stop = STOP_ENERGY_TARGET; break;
            }
        }

        // Force tolerance: LAMMPS TWO norm => Σf² < ftol²
        // Fused fdotf + fdotg reduction
        double fdotf_l = 0.0, fdotg_l = 0.0;
        for (int i = tid; i < N3; i += nthreads) {
            double fi = kd.f[i];
            fdotf_l += fi * fi;
            fdotg_l += fi * kd.g[i];
        }
        double fdotg;
        grid_reduce_sum_2(grid, fdotf_l, fdotg_l, sdata,
                          kd.block_sums_a, kd.block_sums_b, kd.scalars,
                          gnorm2, fdotg);

        if (gnorm2 < d_params.ftol * d_params.ftol) { stop = STOP_FTOL; break; }

        // ================================================================
        //  CG direction update (matches min_cg_kokkos Polak-Ribiere)
        // ================================================================

        double beta = (gg > 1e-30) ? fmax(0.0, (gnorm2 - fdotg) / gg) : 0.0;

        // Periodic CG restart every nlimit iterations
        if (((kd.start_iter + niter + 1) % nlimit) == 0) beta = 0.0;

        gg = gnorm2;

        for (int i = tid; i < N3; i += nthreads) {
            kd.g[i] = kd.f[i];
            kd.h[i] = kd.g[i] + beta * kd.h[i];
        }
        grid.sync();

        // Check if new direction is downhill; reset to steepest descent if not
        double gdoth_l = 0.0;
        for (int i = tid; i < N3; i += nthreads) gdoth_l += kd.g[i] * kd.h[i];
        double gdoth = grid_reduce_sum(grid, gdoth_l, sdata, kd.block_sums_a, kd.scalars, 0);
        if (gdoth <= 0.0) {
            for (int i = tid; i < N3; i += nthreads) kd.h[i] = kd.g[i];
            grid.sync();
        }

        // Displacement check for neighbor rebuild
        if ((niter & 3) == 3) {
            double dmax_l = 0.0;
            for (int i = tid; i < kd.N; i += nthreads) {
                double dx = kd.x[3*i]   - kd.x_build[3*i];
                double dy = kd.x[3*i+1] - kd.x_build[3*i+1];
                double dz = kd.x[3*i+2] - kd.x_build[3*i+2];
                dmax_l = fmax(dmax_l, dx*dx + dy*dy + dz*dz);
            }
            double max_disp2 = grid_reduce_max(grid, dmax_l, sdata,
                                               kd.block_sums_b, kd.scalars, 1);
            if (max_disp2 > kd.skin * kd.skin * 0.25) {
                stop = STOP_NEEDS_REBUILD;
                break;
            }
        }
    }

    if (stop == STOP_RUNNING) stop = STOP_MAXITER;
    if (tid == 0) {
        *kd.d_stop         = stop;
        *kd.d_niter        = kd.start_iter + niter;
        *kd.d_neval        = neval;
        *kd.d_energy_final = ecurrent;
        *kd.d_fnorm        = sqrt(gnorm2 * inv_N);
        *kd.d_gg           = gg;
    }
}

// ============================================================================
//  Host: build per-atom topology tables (CSR format)
// ============================================================================

struct AtomTopology {
    std::vector<int> bond_off, bond_data;
    std::vector<int> angle_off, angle_data;
};

static AtomTopology build_atom_topology(const FusedMinSystem *sys) {
    AtomTopology topo;
    int N = sys->N;

    // Count bonds per atom
    std::vector<int> bond_count(N, 0);
    for (int b = 0; b < sys->nbonds; b++) {
        bond_count[sys->bond_i[b]]++;
        bond_count[sys->bond_j[b]]++;
    }

    topo.bond_off.resize(N + 1);
    topo.bond_off[0] = 0;
    for (int i = 0; i < N; i++) topo.bond_off[i+1] = topo.bond_off[i] + bond_count[i];
    int total_bond_entries = topo.bond_off[N];
    topo.bond_data.resize(2 * total_bond_entries);

    std::vector<int> bpos(N, 0);
    for (int b = 0; b < sys->nbonds; b++) {
        int ai = sys->bond_i[b], aj = sys->bond_j[b], bt = sys->bond_type[b];
        int idx;
        idx = topo.bond_off[ai] + bpos[ai]++;
        topo.bond_data[2*idx]     = bt;
        topo.bond_data[2*idx + 1] = aj;
        idx = topo.bond_off[aj] + bpos[aj]++;
        topo.bond_data[2*idx]     = bt;
        topo.bond_data[2*idx + 1] = ai;
    }

    // Count angles per atom
    std::vector<int> angle_count(N, 0);
    for (int a = 0; a < sys->nangles; a++) {
        angle_count[sys->angle_i[a]]++;
        angle_count[sys->angle_j[a]]++;
        angle_count[sys->angle_k[a]]++;
    }

    topo.angle_off.resize(N + 1);
    topo.angle_off[0] = 0;
    for (int i = 0; i < N; i++) topo.angle_off[i+1] = topo.angle_off[i] + angle_count[i];
    int total_angle_entries = topo.angle_off[N];
    topo.angle_data.resize(4 * total_angle_entries);

    std::vector<int> apos(N, 0);
    for (int a = 0; a < sys->nangles; a++) {
        int at = sys->angle_type[a];
        int ai = sys->angle_i[a], aj = sys->angle_j[a], ak = sys->angle_k[a];
        int atoms[3] = {ai, aj, ak};
        for (int w = 0; w < 3; w++) {
            int idx = topo.angle_off[atoms[w]] + apos[atoms[w]]++;
            topo.angle_data[4*idx]     = at;
            topo.angle_data[4*idx + 1] = ai;
            topo.angle_data[4*idx + 2] = aj;
            topo.angle_data[4*idx + 3] = ak;
        }
    }

    return topo;
}

// ============================================================================
//  Host: build exclusion list in CSR format (computed once per minimize call)
// ============================================================================

struct ExclusionCSR {
    std::vector<int> offsets;  // [N+1]
    std::vector<int> indices;  // sorted excluded atom indices
};

static ExclusionCSR build_exclusion_csr(const FusedMinSystem *sys) {
    int N = sys->N;

    std::vector<std::vector<int>> adj(N);
    for (int b = 0; b < sys->nbonds; b++) {
        adj[sys->bond_i[b]].push_back(sys->bond_j[b]);
        adj[sys->bond_j[b]].push_back(sys->bond_i[b]);
    }

    // Match LAMMPS "special_bonds angle yes" with default LJ weights (0 0 0):
    //   1-2: always excluded
    //   1-3: excluded ONLY if the pair is the first/last atoms of an angle
    //   1-4: always excluded
    // Pre-replication all 1-3 pairs are in angles (linear chain), so behavior
    // matches the old code. Post-replication, fork atoms have 3 bonded neighbors
    // creating 1-3 pairs NOT in any angle — those must get pair interactions.
    std::set<std::pair<int,int>> angle_13;
    for (int a = 0; a < sys->nangles; a++) {
        angle_13.insert({sys->angle_i[a], sys->angle_k[a]});
        angle_13.insert({sys->angle_k[a], sys->angle_i[a]});
    }

    std::vector<std::set<int>> excl(N);
    for (int i = 0; i < N; i++) {
        for (int j : adj[i]) {
            excl[i].insert(j);                               // 1-2: always excluded
            for (int k : adj[j]) {
                if (k == i) continue;
                if (angle_13.count({i, k}))
                    excl[i].insert(k);                       // 1-3: only if in an angle
                for (int l : adj[k]) {
                    if (l == j) continue;
                    excl[i].insert(l);                       // 1-4: always excluded
                }
            }
        }
    }

    ExclusionCSR csr;
    csr.offsets.resize(N + 1);
    csr.offsets[0] = 0;
    for (int i = 0; i < N; i++)
        csr.offsets[i + 1] = csr.offsets[i] + (int)excl[i].size();
    csr.indices.resize(csr.offsets[N]);
    for (int i = 0; i < N; i++) {
        int pos = csr.offsets[i];
        for (int j : excl[i])   // std::set iterates in sorted order
            csr.indices[pos++] = j;
    }

    return csr;
}

// ============================================================================
//  Host: default parameters for DNA hard-harmonic potential
// ============================================================================

void fused_min_init_params(FusedMinParams *p, int pair_style, int bond_style,
                           int profile) {
    memset(p, 0, sizeof(FusedMinParams));

    double kBT    = 0.616032;
    // protein_science lmp.DNA_physical_params: r_bdry = 200 (NOT 85; the 85 here
    // was tuned for the btree_chromo_optmiz 4DWCM model). For the TOPO pair model
    // below r_bdry is unused (boundary uses ribo-sized soft cutoffs), but keep it
    // faithful for the HARD/SOFT-derived WCA sigma table.
    double r_mono = 17.0, r_ribo = 100.0, r_bdry = 200.0;
    double eps_LJ = kBT;
    double eps_soft = kBT;
    double eps_soft_smc = 10.0 * kBT;

    double s_mm = 2*r_mono, s_mr = r_mono+r_ribo, s_rr = 2*r_ribo;
    double s_mb = r_mono+r_bdry, s_rb = r_ribo+r_bdry, s_bb = 2*r_bdry;

    // Build sigma table (1-indexed)
    double sig[9][9];
    memset(sig, 0, sizeof(sig));
    for (int i = 1; i <= 8; i++)
        for (int j = 1; j <= 8; j++) sig[i][j] = s_mm;
    for (int j = 3; j <= 8; j++) { sig[2][j] = s_mr; sig[j][2] = s_mr; }
    sig[2][2] = s_rr;
    for (int j = 3; j <= 8; j++) { sig[1][j] = s_mb; sig[j][1] = s_mb; }
    sig[1][2] = s_rb; sig[2][1] = s_rb;
    sig[1][1] = s_bb;
    for (int i = 3; i <= 6; i++) { sig[i][7] = s_rr; sig[7][i] = s_rr; }

    // Precompute WCA derived quantities (always available)
    for (int i = 1; i <= 8; i++)
        for (int j = i; j <= 8; j++) {
            double s = sig[i][j], s2 = s*s, s6 = s2*s2*s2;
            p->sigma_sq[i][j] = p->sigma_sq[j][i] = s2;
            p->sigma_6[i][j]  = p->sigma_6[j][i]  = s6;
            p->epsilon[i][j]  = p->epsilon[j][i]   = eps_LJ;
        }
    p->epsilon[1][1] = 0.0;

    // ---- Pair mode & cutoffs ----
    if (pair_style == FUSED_PAIR_HARD) {
        for (int i = 1; i <= 8; i++)
            for (int j = 1; j <= 8; j++) {
                p->pair_mode[i][j] = 1; // WCA
                double rc = sig[i][j] * SIXRT2;
                p->cutoff_sq[i][j] = rc * rc;
            }
        p->pair_mode[1][1] = 0;  // bdry-bdry off
    }
    else if (pair_style == FUSED_PAIR_SOFT) {
        for (int i = 1; i <= 8; i++)
            for (int j = 1; j <= 8; j++) {
                if (i == 1 && j == 1) {
                    p->pair_mode[i][j] = 0; // bdry-bdry off
                } else if (i == 1 || j == 1) {
                    p->pair_mode[i][j] = 1; // boundary: WCA
                    double rc = sig[i][j] * SIXRT2;
                    p->cutoff_sq[i][j] = rc * rc;
                } else {
                    p->pair_mode[i][j] = 2; // non-boundary: soft
                    double rc = sig[i][j];
                    p->cutoff_sq[i][j] = rc * rc;
                    p->soft_A[i][j]  = eps_soft;
                    p->soft_rc[i][j] = rc;
                }
            }
        // SMC anchor override: 3-6 vs 7 use 10x soft with ribo-ribo sigma
        for (int i = 3; i <= 6; i++) {
            p->soft_A[i][7] = p->soft_A[7][i] = eps_soft_smc;
            p->soft_rc[i][7] = p->soft_rc[7][i] = s_rr;
            p->cutoff_sq[i][7] = p->cutoff_sq[7][i] = s_rr * s_rr;
        }
    }
    else if (pair_style == FUSED_PAIR_TOPO) {
        // protein_science topoisomerase pair model (lmp.DNA_pair_topo_kk).
        // EVERY pair is `soft` with A = epsilon_soft (= kBT), with ONE exception:
        //   - bdry-bdry (1,1): off (pair_coeff 1 1 lj/cut 0.0 ...).
        // Soft cutoff = the protein_science sigma for that pair class:
        //   DNA-DNA          (3-8 x 3-8): rc = sigma_mono_mono (= 2*r_mono   = 34)
        //   bdry-DNA         (1   x 3-8): rc = sigma_mono_bdry (= r_mono+r_bdry = 217)
        //   ribo-DNA         (2   x 3-8): rc = sigma_mono_ribo (= r_mono+r_ribo = 117)
        //   {bdry|ribo} pair (excl 1,1) : rc = sigma_ribo_ribo (= 2*r_ribo   = 200)
        //
        // bdry-DNA must use r_bdry, not r_ribo. soft/kk is finite at contact,
        // so the undersized 117 cutoff left a cheap path through the envelope
        // and DNA leaked out during replication.
        //
        // CRITICAL: boundary (type 1) is SOFT here, NOT WCA. The previous WCA
        // boundary diverges as r->0 and minimizes to a config that detonates
        // protein_science's soft-only run_soft_harmonic BD (cudaErrorIllegalAddress).
        // Ribosomes (type 2) are also soft here, not off.
        // type-9 atoms remain off by default (pair_mode memset 0; not reached
        // while NUM_TYPES <= 8), matching `pair_coeff 1*9 9 soft 0.0`.
        for (int i = 1; i <= 8; i++)
            for (int j = 1; j <= 8; j++) {
                if (i == 1 && j == 1) {
                    p->pair_mode[i][j] = 0;  // bdry-bdry off
                    continue;
                }
                double rc;
                bool bdry_dna = (i == 1 && j >= 3) || (j == 1 && i >= 3);
                if (i >= 3 && j >= 3)      rc = s_mm;  // DNA-DNA
                else if (bdry_dna)         rc = s_mb;  // bdry-DNA
                else if (i < 3 && j < 3)   rc = s_rr;  // (1,2)/(2,1)/(2,2)
                else                       rc = s_mr;  // ribo-DNA
                p->pair_mode[i][j]  = 2;   // soft
                p->cutoff_sq[i][j]  = rc * rc;
                p->soft_A[i][j]     = eps_soft;
                p->soft_rc[i][j]    = rc;
            }
    }

    // ---- Bonds ----
    double sigma_mm = s_mm;

    if (bond_style == FUSED_BOND_HARMONIC) {
        p->bond_is_fene[1] = 0;
        p->bond_K[1]  = 1000.0 * kBT / (sigma_mm * sigma_mm);
        p->bond_r0[1] = 2.0 * r_mono - 6.0;
    } else { // FUSED_BOND_FENE
        p->bond_is_fene[1] = 1;
        p->bond_K[1] = 100.0 * kBT / (sigma_mm * sigma_mm);
        double R0 = 1.5 * sigma_mm;
        p->fene_R0_sq[1]   = R0 * R0;
        p->fene_eps[1]     = kBT;
        double fs = sigma_mm;
        p->fene_sigma_6[1] = fs*fs*fs*fs*fs*fs;
        p->fene_inner_sq[1] = fs*fs * SIXRT2*SIXRT2;
    }
    // Loops (type 2) are always harmonic.
    // protein_science lmp.DNA_loops + lmp.loop_bond: k_loop = kBT,
    // r_loop = 2*r_mono - 6 (NOT the 4DWCM 4*kBT/(37.5-sigma)^2, 4*r_mono).
    p->bond_is_fene[2] = 0;
    p->bond_K[2]  = kBT;
    p->bond_r0[2] = 2.0 * r_mono - 6.0;

    // Replication fork (type 3): protein_science lmp.loop_bond `bond_coeff 3
    // k_stretch 0.0` -- harmonic with the backbone stiffness and zero rest
    // length (pulls newly-replicated daughter monomers onto their template).
    p->bond_is_fene[3] = 0;
    p->bond_K[3]  = 1000.0 * kBT / (sigma_mm * sigma_mm);
    p->bond_r0[3] = 0.0;

    // ---- Angles (cosine, same for all variants) ----
    double K_bend = kBT * 450.0 / sigma_mm;
    for (int t = 1; t <= 4; t++) p->angle_K[t] = K_bend;

    // ---- Spring/self for ribosomes ----
    p->spring_K_ribo = 1.0 * kBT;

    // ---- CG parameters ----
    // FULL profile: tight, matches protein_science LAMMPS minimize 1e-5 1e-7 40000 400000.
    p->etol    = 1.0e-5;
    p->ftol    = 1.0e-7;
    p->maxiter = 40000;
    p->maxeval = 400000;
    p->dmax    = 0.1;
    p->skin    = 3.0;
    p->energy_target = 0.0;
    p->early_stop_rtol = 0.0;

    // REFRESH profile: loose, iteration-capped budget for in-loop minimizes
    // inside simulator_run_loops where the following BD run re-thermalizes.
    //
    // IMPORTANT: early_stop_rtol is intentionally disabled. An earlier version
    // (job 1513) used early_stop_rtol=0.02 + maxiter=2000 and crashed with
    // "Bond atoms missing at step 5" because the per-iteration hard_harmonic
    // call can start at E~1e10 right after update_loop_bonds() creates
    // freshly-stretched SMC loop bonds. A ~91% relative drop per check still
    // leaves absolute E at ~1e9, i.e. bonds are still catastrophically
    // stretched, but early-stop triggers and hands BD an unstable config.
    //
    // Safer design: pure iteration cap + modestly looser tolerances.
    // FULL maxiter is 40000; log histograms (job 1301) showed hard_harmonic
    // averaging ~3068 iters and topoDNA_harmonic ~6981 iters. Capping at 8000
    // covers the typical case with headroom while still bounding worst-case
    // 36009-iter calls to ~4.5x fewer iters.
    if (profile == FUSED_MIN_REFRESH) {
        p->etol    = 1.0e-4;
        p->ftol    = 1.0e-6;
        p->maxiter = 8000;
        p->maxeval = 80000;
        p->dmax    = 0.1;
        p->early_stop_rtol = 0.0;   // disabled: wrong metric for steep basins
    }

    // REFRESH_SHORT profile: shorter-budget per-iteration in-loop minimize
    // for `minimize_hard_harmonic` ONLY (the dominant cost in run_loops).
    //
    // PATH 2 RETUNE (after job 1601 crashed at maxiter=1500):
    //   The original 1500 cap was unsafe -- under-converged harmonic CG left
    //   FENE bonds slightly past R0 and the LAMMPS FENE polish (hardcoded
    //   maxiter=500) couldn't recover before BD broke the bond.
    //
    //   Bumping maxiter to 4000 (etol/ftol matched to REFRESH so we don't
    //   stop early on a partial basin). 4000 is well above any tail seen in
    //   1514's hard_harmonic histogram, so the FENE polish always gets a
    //   well-converged input.
    //
    //   Scope is also narrowed (in LAMMPS_simulator.cpp): REFRESH_SHORT is
    //   only used for minimize_hard_harmonic on iter 2+ of run_loops. The
    //   topo block (topoDNA_*, soft_*) always uses REFRESH (8000) because
    //   that's where the steepest basins live (post-topoisomerase action).
    if (profile == FUSED_MIN_REFRESH_SHORT) {
        p->etol    = 1.0e-4;     // matched to REFRESH (was 1e-3)
        p->ftol    = 1.0e-6;     // matched to REFRESH (was 1e-5)
        p->maxiter = 4000;       // PATH 2: was 1500 (caused job 1601 crash)
        p->maxeval = 40000;
        p->dmax    = 0.1;
        p->early_stop_rtol = 0.0;
    }
}

void fused_min_default_params(FusedMinParams *p) {
    fused_min_init_params(p, FUSED_PAIR_HARD, FUSED_BOND_HARMONIC);
}

// ============================================================================
//  Host: topology fingerprint (FNV-1a hash of bond/angle connectivity)
// ============================================================================

static size_t compute_topo_hash(const FusedMinSystem *sys) {
    size_t h = 14695981039346656037ULL; // FNV offset basis
    auto mix = [&](const void *data, int count) {
        const unsigned char *p = (const unsigned char *)data;
        for (int i = 0; i < count; i++) {
            h ^= p[i];
            h *= 1099511628211ULL; // FNV prime
        }
    };
    mix(&sys->N,       sizeof(int));
    mix(&sys->nbonds,  sizeof(int));
    mix(&sys->nangles, sizeof(int));
    if (sys->nbonds > 0) {
        mix(sys->bond_type, sys->nbonds * sizeof(int));
        mix(sys->bond_i,    sys->nbonds * sizeof(int));
        mix(sys->bond_j,    sys->nbonds * sizeof(int));
    }
    if (sys->nangles > 0) {
        mix(sys->angle_type, sys->nangles * sizeof(int));
        mix(sys->angle_i,    sys->nangles * sizeof(int));
        mix(sys->angle_j,    sys->nangles * sizeof(int));
        mix(sys->angle_k,    sys->nangles * sizeof(int));
    }
    return h;
}

// ============================================================================
//  GPU neighbor list build kernels
// ============================================================================

__global__ void gpu_compute_cell_idx_kernel(
    const double *x, int N,
    int ncx, int ncy, int ncz,
    double cx, double cy, double cz,
    double blo_x, double blo_y, double blo_z,
    int *cell_of_atom)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    int ix = min(max((int)((x[3*i]   - blo_x) / cx), 0), ncx - 1);
    int iy = min(max((int)((x[3*i+1] - blo_y) / cy), 0), ncy - 1);
    int iz = min(max((int)((x[3*i+2] - blo_z) / cz), 0), ncz - 1);
    cell_of_atom[i] = ix * ncy * ncz + iy * ncz + iz;
}

__global__ void gpu_count_cell_atoms_kernel(
    const int *cell_of_atom, int N, int *cell_count)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    atomicAdd(&cell_count[cell_of_atom[i]], 1);
}

__global__ void gpu_scatter_atoms_kernel(
    const int *cell_of_atom, int N,
    const int *cell_start, int *cell_fill,
    int *sorted_atoms)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    int cell = cell_of_atom[i];
    int pos = cell_start[cell] + atomicAdd(&cell_fill[cell], 1);
    sorted_atoms[pos] = i;
}

__device__ __forceinline__
bool is_excluded(const int *excl_list, int start, int end, int j) {
    while (start < end) {
        int mid = (start + end) >> 1;
        int val = excl_list[mid];
        if (val == j) return true;
        if (val < j) start = mid + 1;
        else end = mid;
    }
    return false;
}

__global__ void gpu_build_nlist_kernel(
    const double *x, const int *type, int N,
    const int *sorted_atoms, const int *cell_start,
    int ncx, int ncy, int ncz,
    double cx, double cy, double cz,
    double blo_x, double blo_y, double blo_z,
    double rcut2,
    const int *excl_off, const int *excl_list,
    int *numneigh, int *neighlist, int max_neighs,
    int nlist_stride)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    int ti = type[i];
    double xi = x[3*i], yi = x[3*i+1], zi = x[3*i+2];
    int nn = 0;

    int ix = min(max((int)((xi - blo_x) / cx), 0), ncx - 1);
    int iy = min(max((int)((yi - blo_y) / cy), 0), ncy - 1);
    int iz = min(max((int)((zi - blo_z) / cz), 0), ncz - 1);

    int e_start = excl_off[i];
    int e_end   = excl_off[i + 1];

    for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dz = -1; dz <= 1; dz++) {
        int jx = ix + dx, jy = iy + dy, jz = iz + dz;
        if (jx < 0 || jx >= ncx || jy < 0 || jy >= ncy || jz < 0 || jz >= ncz)
            continue;
        int cidx = jx * ncy * ncz + jy * ncz + jz;
        int cs = cell_start[cidx];
        int ce = cell_start[cidx + 1];

        for (int si = cs; si < ce; si++) {
            int j = sorted_atoms[si];
            if (j == i) continue;
            if (ti == 1 && type[j] == 1) continue;
            if (is_excluded(excl_list, e_start, e_end, j)) continue;

            double ddx = xi - x[3*j], ddy = yi - x[3*j+1], ddz = zi - x[3*j+2];
            if (ddx*ddx + ddy*ddy + ddz*ddz < rcut2) {
                if (nn < max_neighs)
                    neighlist[(size_t)nn * nlist_stride + i] = j;
                nn++;
            }
        }
    }
    numneigh[i] = (nn <= max_neighs) ? nn : max_neighs;
}

static void gpu_rebuild_neighbor_list(
    double *d_x, int *d_type, int N,
    int *d_excl_off, int *d_excl_list,
    double rcut2,
    const double box_lo[3], const double box_hi[3],
    int ncx, int ncy, int ncz,
    double cx, double cy, double cz,
    int ncells,
    int *d_cell_of_atom, int *d_sorted_atoms,
    int *d_cell_start, int *d_cell_count,
    int *d_numneigh, int *d_neighlist, int max_neighs,
    int nlist_stride)
{
    int threads = 256;
    int blocks_N = (N + threads - 1) / threads;

    gpu_compute_cell_idx_kernel<<<blocks_N, threads>>>(
        d_x, N, ncx, ncy, ncz, cx, cy, cz,
        box_lo[0], box_lo[1], box_lo[2], d_cell_of_atom);

    CUDA_CHECK(cudaMemset(d_cell_count, 0, (ncells + 1) * sizeof(int)));
    gpu_count_cell_atoms_kernel<<<blocks_N, threads>>>(d_cell_of_atom, N, d_cell_count);

    thrust::device_ptr<int> d_count_ptr(d_cell_count);
    thrust::device_ptr<int> d_start_ptr(d_cell_start);
    thrust::exclusive_scan(d_count_ptr, d_count_ptr + ncells + 1, d_start_ptr);

    CUDA_CHECK(cudaMemset(d_cell_count, 0, (ncells + 1) * sizeof(int)));
    gpu_scatter_atoms_kernel<<<blocks_N, threads>>>(
        d_cell_of_atom, N, d_cell_start, d_cell_count, d_sorted_atoms);

    gpu_build_nlist_kernel<<<blocks_N, threads>>>(
        d_x, d_type, N, d_sorted_atoms, d_cell_start,
        ncx, ncy, ncz, cx, cy, cz,
        box_lo[0], box_lo[1], box_lo[2], rcut2,
        d_excl_off, d_excl_list,
        d_numneigh, d_neighlist, max_neighs,
        nlist_stride);

    CUDA_CHECK(cudaDeviceSynchronize());
}

// ============================================================================
//  Host: CPU neighbor list builder (kept for reference, no longer called)
// ============================================================================

void fused_build_neighbor_list(
    const double *x, const int *type, int N,
    const FusedMinSystem *sys,
    double skin, double max_cutoff,
    const double box_lo[3], const double box_hi[3],
    FusedNeighList *nlist)
{
    double rcut = max_cutoff + skin;
    double rcut2 = rcut * rcut;

    // Build adjacency list from bond topology
    std::vector<std::vector<int>> adj(N);
    for (int b = 0; b < sys->nbonds; b++) {
        adj[sys->bond_i[b]].push_back(sys->bond_j[b]);
        adj[sys->bond_j[b]].push_back(sys->bond_i[b]);
    }

    // DEAD CODE (kept for reference) — see build_exclusion_csr for the live version.
    // Exclusion list matching LAMMPS "special_bonds angle yes":
    //   1-2: always excluded, 1-3: only if in angle, 1-4: always excluded
    std::set<std::pair<int,int>> angle_13;
    for (int a = 0; a < sys->nangles; a++) {
        angle_13.insert({sys->angle_i[a], sys->angle_k[a]});
        angle_13.insert({sys->angle_k[a], sys->angle_i[a]});
    }
    std::vector<std::set<int>> excl(N);
    for (int i = 0; i < N; i++) {
        for (int j : adj[i]) {
            excl[i].insert(j);                       // 1-2
            for (int k : adj[j]) {
                if (k == i) continue;
                if (angle_13.count({i, k}))
                    excl[i].insert(k);               // 1-3 (only if in angle)
                for (int l : adj[k]) {
                    if (l == j) continue;
                    excl[i].insert(l);               // 1-4
                }
            }
        }
    }

    // Cell list
    double Lx = box_hi[0] - box_lo[0];
    double Ly = box_hi[1] - box_lo[1];
    double Lz = box_hi[2] - box_lo[2];
    int ncx = std::max(1, (int)(Lx / rcut));
    int ncy = std::max(1, (int)(Ly / rcut));
    int ncz = std::max(1, (int)(Lz / rcut));
    double cx = Lx / ncx, cy = Ly / ncy, cz = Lz / ncz;

    std::vector<std::vector<int>> cells(ncx * ncy * ncz);
    for (int i = 0; i < N; i++) {
        int ix = std::min((int)((x[3*i]   - box_lo[0]) / cx), ncx - 1);
        int iy = std::min((int)((x[3*i+1] - box_lo[1]) / cy), ncy - 1);
        int iz = std::min((int)((x[3*i+2] - box_lo[2]) / cz), ncz - 1);
        ix = std::max(0, ix); iy = std::max(0, iy); iz = std::max(0, iz);
        cells[ix * ncy * ncz + iy * ncz + iz].push_back(i);
    }

    nlist->N = N;
    nlist->max_neighs = 1024;
    nlist->numneigh = (int*)malloc(N * sizeof(int));
    nlist->neighlist = (int*)malloc((size_t)N * nlist->max_neighs * sizeof(int));

    int overflow_count = 0;

    for (int i = 0; i < N; i++) {
        int ti = type[i];
        double xi = x[3*i], yi = x[3*i+1], zi = x[3*i+2];
        int nn = 0;

        int ix = std::min((int)((xi - box_lo[0]) / cx), ncx - 1);
        int iy = std::min((int)((yi - box_lo[1]) / cy), ncy - 1);
        int iz = std::min((int)((zi - box_lo[2]) / cz), ncz - 1);
        ix = std::max(0, ix); iy = std::max(0, iy); iz = std::max(0, iz);

        for (int dx = -1; dx <= 1; dx++)
        for (int dy = -1; dy <= 1; dy++)
        for (int dz = -1; dz <= 1; dz++) {
            int jx = ix + dx, jy = iy + dy, jz = iz + dz;
            if (jx < 0 || jx >= ncx || jy < 0 || jy >= ncy || jz < 0 || jz >= ncz)
                continue;
            int cidx = jx * ncy * ncz + jy * ncz + jz;
            for (int j : cells[cidx]) {
                if (j == i) continue;
                if (ti == 1 && type[j] == 1) continue;
                if (excl[i].count(j)) continue;

                double ddx = xi - x[3*j], ddy = yi - x[3*j+1], ddz = zi - x[3*j+2];
                if (ddx*ddx + ddy*ddy + ddz*ddz < rcut2) {
                    if (nn < nlist->max_neighs) {
                        nlist->neighlist[i * nlist->max_neighs + nn] = j;
                        nn++;
                    } else {
                        overflow_count++;
                    }
                }
            }
        }
        nlist->numneigh[i] = nn;
    }

    if (overflow_count > 0)
        fprintf(stderr, "WARNING: neighbor list overflow for %d entries (max_neighs=%d)\n",
                overflow_count, nlist->max_neighs);
}

// ============================================================================
//  Host: persistent cache init / destroy
// ============================================================================

void fused_cache_init(FusedDeviceCache *c) {
    memset(c, 0, sizeof(FusedDeviceCache));
}

static void cache_free_n_arrays(FusedDeviceCache *c) {
    if (!c->d_x) return;
    cudaFree(c->d_x);       cudaFree(c->d_f);
    cudaFree(c->d_h);       cudaFree(c->d_g);
    cudaFree(c->d_x0);      cudaFree(c->d_x_build);
    cudaFree(c->d_xref);    cudaFree(c->d_type);
    cudaFree(c->d_numneigh); cudaFree(c->d_neighlist);
    cudaFree(c->d_cell_of_atom); cudaFree(c->d_sorted_atoms);
    c->d_x = nullptr;
}

static void cache_free_topo(FusedDeviceCache *c) {
    if (!c->d_atom_bond_off) return;
    cudaFree(c->d_atom_bond_off);  cudaFree(c->d_atom_bond_data);
    cudaFree(c->d_atom_angle_off); cudaFree(c->d_atom_angle_data);
    cudaFree(c->d_excl_off);       cudaFree(c->d_excl_list);
    c->d_atom_bond_off = nullptr;
    c->bond_data_cap = c->angle_data_cap = c->excl_cap = 0;
}

void fused_cache_destroy(FusedDeviceCache *c) {
    if (!c->initialized) return;
    cache_free_n_arrays(c);
    cache_free_topo(c);
    if (c->d_cell_start) { cudaFree(c->d_cell_start); cudaFree(c->d_cell_count); }
    if (c->d_block_sums_a) {
        cudaFree(c->d_block_sums_a); cudaFree(c->d_block_sums_b);
        cudaFree(c->d_scalars);
        cudaFree(c->d_stop);  cudaFree(c->d_niter);  cudaFree(c->d_neval);
        cudaFree(c->d_energy_init); cudaFree(c->d_energy_final); cudaFree(c->d_fnorm);
        cudaFree(c->d_gg);
    }
    memset(c, 0, sizeof(FusedDeviceCache));
}

// ============================================================================
//  Host: main CG minimizer wrapper (with persistent cache)
// ============================================================================

FusedMinResult fused_cg_minimize(
    FusedMinSystem *sys,
    const FusedMinParams *params,
    const double box_lo[3], const double box_hi[3],
    int gpu_id,
    FusedDeviceCache *cache)
{
    CUDA_CHECK(cudaSetDevice(gpu_id));

    int N = sys->N;
    int N3 = N * 3;
    FusedMinResult result;
    memset(&result, 0, sizeof(result));

    // ================================================================
    //  Phase A: ensure device memory is allocated (only when N grows)
    // ================================================================
    bool first_call = !cache->initialized;
    bool n_grew     = (N > cache->N_cap);

    if (first_call || n_grew) {
        if (n_grew && cache->initialized) {
            cache_free_n_arrays(cache);
            cache_free_topo(cache);
        }

        int N_alloc = N * 3 / 2;  // 1.5x headroom
        size_t n3a = (size_t)N_alloc * 3 * sizeof(double);
        size_t na  = (size_t)N_alloc * sizeof(int);

        CUDA_CHECK(cudaMalloc(&cache->d_x,       n3a));
        CUDA_CHECK(cudaMalloc(&cache->d_f,       n3a));
        CUDA_CHECK(cudaMalloc(&cache->d_h,       n3a));
        CUDA_CHECK(cudaMalloc(&cache->d_g,       n3a));
        CUDA_CHECK(cudaMalloc(&cache->d_x0,      n3a));
        CUDA_CHECK(cudaMalloc(&cache->d_x_build, n3a));
        CUDA_CHECK(cudaMalloc(&cache->d_xref,    n3a));
        CUDA_CHECK(cudaMalloc(&cache->d_type,    na));

        cache->max_neighs = 1024;
        CUDA_CHECK(cudaMalloc(&cache->d_numneigh,  na));
        CUDA_CHECK(cudaMalloc(&cache->d_neighlist,
                              (size_t)N_alloc * cache->max_neighs * sizeof(int)));

        CUDA_CHECK(cudaMalloc(&cache->d_cell_of_atom, na));
        CUDA_CHECK(cudaMalloc(&cache->d_sorted_atoms, na));

        cache->N_cap = N_alloc;

        if (first_call) {
            int block_size = 256;
            int num_blocks_per_sm = 0;
            CUDA_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                &num_blocks_per_sm, fused_cg_kernel, block_size,
                block_size * sizeof(double)));
            int device;
            CUDA_CHECK(cudaGetDevice(&device));
            cudaDeviceProp prop;
            CUDA_CHECK(cudaGetDeviceProperties(&prop, device));
            cache->max_blocks = num_blocks_per_sm * prop.multiProcessorCount;
            cache->block_size = block_size;
            cache->gpu_id = gpu_id;

            CUDA_CHECK(cudaMalloc(&cache->d_block_sums_a,
                                  cache->max_blocks * sizeof(double)));
            CUDA_CHECK(cudaMalloc(&cache->d_block_sums_b,
                                  cache->max_blocks * sizeof(double)));
            CUDA_CHECK(cudaMalloc(&cache->d_scalars, 16 * sizeof(double)));

            CUDA_CHECK(cudaMalloc(&cache->d_stop,         sizeof(int)));
            CUDA_CHECK(cudaMalloc(&cache->d_niter,        sizeof(int)));
            CUDA_CHECK(cudaMalloc(&cache->d_neval,        sizeof(int)));
            CUDA_CHECK(cudaMalloc(&cache->d_energy_init,  sizeof(double)));
            CUDA_CHECK(cudaMalloc(&cache->d_energy_final, sizeof(double)));
            CUDA_CHECK(cudaMalloc(&cache->d_fnorm,        sizeof(double)));
            CUDA_CHECK(cudaMalloc(&cache->d_gg,           sizeof(double)));

            cache->cell_cap = 0;
        }

        cache->initialized = true;
        cache->cached_nbonds  = -1;
        cache->cached_nangles = -1;
        cache->topo_hash      = 0;

        fprintf(stderr, "Fused CG: allocated for N_cap=%d (N=%d), blocks=%d, threads/block=%d\n",
                cache->N_cap, N, cache->max_blocks, cache->block_size);
    }

    int max_blocks = cache->max_blocks;
    int block_size = cache->block_size;

    // ================================================================
    //  Phase B: rebuild topology/exclusion CSR only when topology changes
    // ================================================================
    using clock = std::chrono::steady_clock;
    double topo_host_ms = 0.0;
    int topo_rebuilt = 0;

    clock::time_point t_topo0 = clock::now();
    size_t new_hash = compute_topo_hash(sys);
    bool topo_changed = (new_hash != cache->topo_hash) ||
                        !cache->d_atom_bond_off;
    if (topo_changed) {
        clock::time_point t_build0 = clock::now();
        AtomTopology topo = build_atom_topology(sys);
        ExclusionCSR excl_csr = build_exclusion_csr(sys);

        int bond_entries  = topo.bond_off[N];
        int angle_entries = topo.angle_off[N];
        int total_excl    = excl_csr.offsets[N];

        // Bond CSR: realloc if capacity exceeded
        int need_bond = 2 * bond_entries;
        if (need_bond > cache->bond_data_cap || !cache->d_atom_bond_off) {
            if (cache->d_atom_bond_off) {
                cudaFree(cache->d_atom_bond_off);
                cudaFree(cache->d_atom_bond_data);
            }
            int alloc_off = cache->N_cap + 1;
            int alloc_data = std::max(need_bond, need_bond * 3 / 2);
            CUDA_CHECK(cudaMalloc(&cache->d_atom_bond_off,  alloc_off * sizeof(int)));
            CUDA_CHECK(cudaMalloc(&cache->d_atom_bond_data, alloc_data * sizeof(int)));
            cache->bond_data_cap = alloc_data;
        }
        CUDA_CHECK(cudaMemcpy(cache->d_atom_bond_off,  topo.bond_off.data(),
                              (N+1)*sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(cache->d_atom_bond_data, topo.bond_data.data(),
                              need_bond*sizeof(int), cudaMemcpyHostToDevice));

        // Angle CSR: realloc if capacity exceeded
        int need_angle = 4 * angle_entries;
        if (need_angle > cache->angle_data_cap || !cache->d_atom_angle_off) {
            if (cache->d_atom_angle_off) {
                cudaFree(cache->d_atom_angle_off);
                cudaFree(cache->d_atom_angle_data);
            }
            int alloc_off = cache->N_cap + 1;
            int alloc_data = std::max(need_angle, need_angle * 3 / 2);
            CUDA_CHECK(cudaMalloc(&cache->d_atom_angle_off,  alloc_off * sizeof(int)));
            CUDA_CHECK(cudaMalloc(&cache->d_atom_angle_data, alloc_data * sizeof(int)));
            cache->angle_data_cap = alloc_data;
        }
        CUDA_CHECK(cudaMemcpy(cache->d_atom_angle_off,  topo.angle_off.data(),
                              (N+1)*sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(cache->d_atom_angle_data, topo.angle_data.data(),
                              need_angle*sizeof(int), cudaMemcpyHostToDevice));

        // Exclusion CSR: realloc if capacity exceeded
        int need_excl = std::max(1, total_excl);
        if (need_excl > cache->excl_cap || !cache->d_excl_off) {
            if (cache->d_excl_off) {
                cudaFree(cache->d_excl_off);
                cudaFree(cache->d_excl_list);
            }
            int alloc_off = cache->N_cap + 1;
            int alloc_list = std::max(need_excl, need_excl * 3 / 2);
            CUDA_CHECK(cudaMalloc(&cache->d_excl_off,  alloc_off * sizeof(int)));
            CUDA_CHECK(cudaMalloc(&cache->d_excl_list, alloc_list * sizeof(int)));
            cache->excl_cap = alloc_list;
        }
        CUDA_CHECK(cudaMemcpy(cache->d_excl_off, excl_csr.offsets.data(),
                              (N+1)*sizeof(int), cudaMemcpyHostToDevice));
        if (total_excl > 0)
            CUDA_CHECK(cudaMemcpy(cache->d_excl_list, excl_csr.indices.data(),
                                  total_excl*sizeof(int), cudaMemcpyHostToDevice));

        cache->cached_nbonds  = sys->nbonds;
        cache->cached_nangles = sys->nangles;
        cache->N_cur          = N;
        cache->topo_hash      = new_hash;
        topo_host_ms = fused_ms_since(t_build0);
        topo_rebuilt = 1;
        fprintf(stderr, "  Topology rebuilt (hash=%zu)\n", new_hash);
    } // end if (topo_changed)
    double phaseB_ms = fused_ms_since(t_topo0);

    // ================================================================
    //  Phase C: per-call work (always)
    // ================================================================

    size_t n3_bytes = N3 * sizeof(double);
    size_t n_bytes  = N * sizeof(int);

    clock::time_point t_h2d0 = clock::now();
    CUDA_CHECK(cudaMemcpyToSymbol(d_params, params, sizeof(FusedMinParams)));
    CUDA_CHECK(cudaMemcpy(cache->d_x,    sys->x,    n3_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(cache->d_xref, sys->xref, n3_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(cache->d_type, sys->type, n_bytes,  cudaMemcpyHostToDevice));
    double h2d_ms = fused_ms_since(t_h2d0);

    // Build KernelData from cache pointers
    KernelData kd;
    kd.N = N;
    kd.max_neighs      = cache->max_neighs;
    kd.neighlist_stride = cache->N_cap;
    kd.skin         = params->skin;
    kd.start_iter   = 0;
    kd.start_neval  = 0;
    kd.x = cache->d_x;  kd.f = cache->d_f;  kd.h = cache->d_h;  kd.g = cache->d_g;
    kd.x0 = cache->d_x0;  kd.x_build = cache->d_x_build;
    kd.xref = cache->d_xref;  kd.type = cache->d_type;
    kd.atom_bond_off  = cache->d_atom_bond_off;
    kd.atom_bond_data = cache->d_atom_bond_data;
    kd.atom_angle_off  = cache->d_atom_angle_off;
    kd.atom_angle_data = cache->d_atom_angle_data;
    kd.numneigh  = cache->d_numneigh;
    kd.neighlist = cache->d_neighlist;
    kd.block_sums_a = cache->d_block_sums_a;
    kd.block_sums_b = cache->d_block_sums_b;
    kd.scalars      = cache->d_scalars;
    kd.d_stop         = cache->d_stop;
    kd.d_niter        = cache->d_niter;
    kd.d_neval        = cache->d_neval;
    kd.d_energy_init  = cache->d_energy_init;
    kd.d_energy_final = cache->d_energy_final;
    kd.d_fnorm        = cache->d_fnorm;
    kd.d_gg           = cache->d_gg;

    // Cell grid parameters (rcut depends on pair_style)
    double max_cut = 0.0;
    for (int i = 1; i <= 8; i++)
        for (int j = 1; j <= 8; j++)
            if (params->pair_mode[i][j] != 0)
                max_cut = std::max(max_cut, sqrt(params->cutoff_sq[i][j]));

    double rcut = max_cut + params->skin;
    double rcut2 = rcut * rcut;
    double Lx = box_hi[0] - box_lo[0];
    double Ly = box_hi[1] - box_lo[1];
    double Lz = box_hi[2] - box_lo[2];
    int ncx = std::max(1, (int)(Lx / rcut));
    int ncy = std::max(1, (int)(Ly / rcut));
    int ncz = std::max(1, (int)(Lz / rcut));
    int ncells = ncx * ncy * ncz;
    double cx = Lx / ncx, cy = Ly / ncy, cz = Lz / ncz;

    // Ensure cell list arrays are large enough
    if (ncells + 1 > cache->cell_cap) {
        if (cache->d_cell_start) {
            cudaFree(cache->d_cell_start);
            cudaFree(cache->d_cell_count);
        }
        int alloc_cells = std::max(ncells + 1, (ncells + 1) * 3 / 2);
        CUDA_CHECK(cudaMalloc(&cache->d_cell_start, alloc_cells * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&cache->d_cell_count, alloc_cells * sizeof(int)));
        cache->cell_cap = alloc_cells;
    }

    // Timing
    cudaEvent_t t_start, t_stop;
    CUDA_CHECK(cudaEventCreate(&t_start));
    CUDA_CHECK(cudaEventCreate(&t_stop));
    CUDA_CHECK(cudaEventRecord(t_start));

    // Multi-launch loop (exits and re-enters on neighbor rebuild)
    int total_iter = 0, total_neval = 0;
    int h_stop = STOP_NEEDS_REBUILD;
    int launch_count = 0;
    double nlist_time_ms = 0;

    while (h_stop == STOP_NEEDS_REBUILD) {
        cudaEvent_t nl_start, nl_stop;
        CUDA_CHECK(cudaEventCreate(&nl_start));
        CUDA_CHECK(cudaEventCreate(&nl_stop));
        CUDA_CHECK(cudaEventRecord(nl_start));

        gpu_rebuild_neighbor_list(
            kd.x, kd.type, N,
            cache->d_excl_off, cache->d_excl_list,
            rcut2, box_lo, box_hi,
            ncx, ncy, ncz, cx, cy, cz, ncells,
            cache->d_cell_of_atom, cache->d_sorted_atoms,
            cache->d_cell_start, cache->d_cell_count,
            kd.numneigh, kd.neighlist, kd.max_neighs,
            cache->N_cap);

        CUDA_CHECK(cudaEventRecord(nl_stop));
        CUDA_CHECK(cudaEventSynchronize(nl_stop));
        float nl_ms = 0;
        CUDA_CHECK(cudaEventElapsedTime(&nl_ms, nl_start, nl_stop));
        nlist_time_ms += nl_ms;
        CUDA_CHECK(cudaEventDestroy(nl_start));
        CUDA_CHECK(cudaEventDestroy(nl_stop));

        kd.start_iter  = total_iter;
        kd.start_neval = total_neval;
        kd.resume      = (launch_count > 0) ? 1 : 0;

        dim3 grid_dim(max_blocks), block_dim(block_size);
        void *args[] = { &kd };
        CUDA_CHECK(cudaLaunchCooperativeKernel(
            (void*)fused_cg_kernel, grid_dim, block_dim,
            args, block_size * sizeof(double), 0));
        CUDA_CHECK(cudaDeviceSynchronize());

        if (launch_count == 0)
            CUDA_CHECK(cudaMemcpy(&result.energy_initial, kd.d_energy_init,
                                  sizeof(double), cudaMemcpyDeviceToHost));

        CUDA_CHECK(cudaMemcpy(&h_stop,       kd.d_stop,         sizeof(int),    cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(&total_iter,   kd.d_niter,        sizeof(int),    cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(&total_neval,  kd.d_neval,        sizeof(int),    cudaMemcpyDeviceToHost));
        launch_count++;

        if (h_stop == STOP_NEEDS_REBUILD)
            fprintf(stderr, "  Neighbor rebuild at iter %d (neval=%d)\n", total_iter, total_neval);
    }

    CUDA_CHECK(cudaEventRecord(t_stop));
    CUDA_CHECK(cudaEventSynchronize(t_stop));
    float ms = 0;
    CUDA_CHECK(cudaEventElapsedTime(&ms, t_start, t_stop));

    fprintf(stderr, "  %d kernel launches, %.1f ms neighbor list builds (%.1f%% of total)\n",
            launch_count, nlist_time_ms, 100.0 * nlist_time_ms / ms);

    if (fused_profile_enabled()) {
        double cg_ms = ms - nlist_time_ms;
        fprintf(stderr,
                "[FUSED_PROFILE] N=%d topo_rebuilt=%d phaseB_ms=%.2f topo_build_cpu_ms=%.2f "
                "h2d_ms=%.2f gpu_total_ms=%.2f nlist_ms=%.2f cg_kernel_ms=%.2f launches=%d "
                "niter=%d neval=%d\n",
                N, topo_rebuilt, phaseB_ms, topo_host_ms, h2d_ms, ms, nlist_time_ms, cg_ms,
                launch_count, total_iter, total_neval);
    }

    // Copy final positions back
    CUDA_CHECK(cudaMemcpy(sys->x, kd.x, n3_bytes, cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaMemcpy(&result.energy_final,   kd.d_energy_final, sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(&result.fnorm_final,    kd.d_fnorm,        sizeof(double), cudaMemcpyDeviceToHost));
    result.stop       = (StopCode)h_stop;
    result.niter      = total_iter;
    result.neval      = total_neval;
    result.elapsed_ms = ms;

    CUDA_CHECK(cudaEventDestroy(t_start));
    CUDA_CHECK(cudaEventDestroy(t_stop));

    return result;
}

// ============================================================================
//  Host: LAMMPS data file reader (atom_style full)
// ============================================================================

int fused_read_lammps_data(const char *filename, FusedMinSystem *sys,
                           double box_lo[3], double box_hi[3]) {
    FILE *fp = fopen(filename, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", filename); return -1; }

    memset(sys, 0, sizeof(FusedMinSystem));
    char line[1024];

    // Parse header — use strstr to verify keywords before sscanf
    while (fgets(line, sizeof(line), fp)) {
        int val;
        double lo, hi;
        if (strstr(line, "atoms") && !strstr(line, "atom types") &&
            sscanf(line, "%d", &val) == 1)
            sys->N = val;
        else if (strstr(line, "bonds") && !strstr(line, "bond types") &&
                 sscanf(line, "%d", &val) == 1)
            sys->nbonds = val;
        else if (strstr(line, "angles") && !strstr(line, "angle types") &&
                 sscanf(line, "%d", &val) == 1)
            sys->nangles = val;
        else if (strstr(line, "xlo xhi") && sscanf(line, "%lf %lf", &lo, &hi) == 2)
            { box_lo[0] = lo; box_hi[0] = hi; }
        else if (strstr(line, "ylo yhi") && sscanf(line, "%lf %lf", &lo, &hi) == 2)
            { box_lo[1] = lo; box_hi[1] = hi; }
        else if (strstr(line, "zlo zhi") && sscanf(line, "%lf %lf", &lo, &hi) == 2)
            { box_lo[2] = lo; box_hi[2] = hi; }
        if (strstr(line, "Atoms") || strstr(line, "Masses")) break;
    }

    int N = sys->N;
    if (N == 0) { fclose(fp); return -1; }

    sys->x    = (double*)calloc(N * 3, sizeof(double));
    sys->xref = (double*)calloc(N * 3, sizeof(double));
    sys->type = (int*)calloc(N, sizeof(int));

    if (sys->nbonds > 0) {
        sys->bond_type = (int*)malloc(sys->nbonds * sizeof(int));
        sys->bond_i    = (int*)malloc(sys->nbonds * sizeof(int));
        sys->bond_j    = (int*)malloc(sys->nbonds * sizeof(int));
    }
    if (sys->nangles > 0) {
        sys->angle_type = (int*)malloc(sys->nangles * sizeof(int));
        sys->angle_i    = (int*)malloc(sys->nangles * sizeof(int));
        sys->angle_j    = (int*)malloc(sys->nangles * sizeof(int));
        sys->angle_k    = (int*)malloc(sys->nangles * sizeof(int));
    }

    // Seek to sections and parse
    rewind(fp);
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "Atoms", 5) == 0) {
            (void)fgets(line, sizeof(line), fp); // blank line
            for (int n = 0; n < N; n++) {
                if (!fgets(line, sizeof(line), fp)) break;
                int id = 0, mol = 0, tp = 0;
                double xx = 0, yy = 0, zz = 0, q = 0;
                int nc = sscanf(line, "%d %d %d %lf %lf %lf %lf",
                                &id, &mol, &tp, &q, &xx, &yy, &zz);
                if (nc == 6) // atom_style molecular: id mol type x y z
                    sscanf(line, "%d %d %d %lf %lf %lf",
                           &id, &mol, &tp, &xx, &yy, &zz);
                if (nc >= 6) {
                    int idx = id - 1;
                    if (idx >= 0 && idx < N) {
                        sys->type[idx] = tp;
                        sys->x[3*idx]   = xx;
                        sys->x[3*idx+1] = yy;
                        sys->x[3*idx+2] = zz;
                    }
                }
            }
        } else if (strncmp(line, "Bonds", 5) == 0 && sys->nbonds > 0) {
            (void)fgets(line, sizeof(line), fp);
            for (int n = 0; n < sys->nbonds; n++) {
                if (!fgets(line, sizeof(line), fp)) break;
                int id, tp, ai, aj;
                if (sscanf(line, "%d %d %d %d", &id, &tp, &ai, &aj) == 4) {
                    sys->bond_type[n] = tp;
                    sys->bond_i[n]    = ai - 1;
                    sys->bond_j[n]    = aj - 1;
                }
            }
        } else if (strncmp(line, "Angles", 6) == 0 && sys->nangles > 0) {
            (void)fgets(line, sizeof(line), fp);
            for (int n = 0; n < sys->nangles; n++) {
                if (!fgets(line, sizeof(line), fp)) break;
                int id, tp, ai, aj, ak;
                if (sscanf(line, "%d %d %d %d %d", &id, &tp, &ai, &aj, &ak) == 5) {
                    sys->angle_type[n] = tp;
                    sys->angle_i[n]    = ai - 1;
                    sys->angle_j[n]    = aj - 1;
                    sys->angle_k[n]    = ak - 1;
                }
            }
        }
    }

    // Set xref = initial positions (reference for spring/self)
    memcpy(sys->xref, sys->x, N * 3 * sizeof(double));

    fclose(fp);
    fprintf(stderr, "Read %s: %d atoms, %d bonds, %d angles\n",
            filename, N, sys->nbonds, sys->nangles);
    fprintf(stderr, "  Box: [%.1f,%.1f] x [%.1f,%.1f] x [%.1f,%.1f]\n",
            box_lo[0], box_hi[0], box_lo[1], box_hi[1], box_lo[2], box_hi[2]);
    return 0;
}

// ============================================================================
//  Host: cleanup
// ============================================================================

void fused_free_system(FusedMinSystem *sys) {
    free(sys->x);    free(sys->xref);  free(sys->type);
    free(sys->bond_type); free(sys->bond_i); free(sys->bond_j);
    free(sys->angle_type); free(sys->angle_i); free(sys->angle_j); free(sys->angle_k);
    memset(sys, 0, sizeof(FusedMinSystem));
}

void fused_free_neighlist(FusedNeighList *nlist) {
    free(nlist->numneigh);
    free(nlist->neighlist);
    memset(nlist, 0, sizeof(FusedNeighList));
}
