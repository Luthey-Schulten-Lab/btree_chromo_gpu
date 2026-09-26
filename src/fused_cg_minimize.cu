#include "wcm_arch.cuh"
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
#include <thread>
#include <cub/device/device_scan.cuh>
#include <cub/device/device_radix_sort.cuh>

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
// FP32 copies of the pair tables; the pair term runs in FP32 on the FP64 coordinate
// difference and its force/energy are accumulated in FP64. WCM_FP64=1 keeps the FP64 pair term (d_pair_fp64 = 1).
__constant__ float d_pf_csq[9][9], d_pf_eps[9][9], d_pf_sig6[9][9], d_pf_A[9][9], d_pf_invrc[9][9], d_pf_fpre[9][9];
__constant__ int d_pair_fp64;
// neighbour-list cutoff per type pair, (sqrt(cutoff_sq[i][j]) + skin)^2, or -1 for pairs that never interact
// (pair_mode 0). The list used one cutoff for every pair (the largest, bdry-DNA 217 A + skin) although DNA-DNA pairs interact
// only within 34 A. A pair outside its own cutoff + skin cannot come within its cutoff before the next rebuild (same skin/2
// displacement criterion), and compute_atom_forces skips it anyway, so the contributing pairs are the same.
__constant__ double d_nl_cut2[9][9];

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
    const int *order;     // atom handled at grid-stride slot k (balanced permutation, rebuilt with the list)
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

// One pair's contribution (the former loop body, arithmetic unchanged; dx = xa - xj etc.)
__device__ __forceinline__
void pair_term(int ta, int tj, double dx, double dy, double dz,
               double &fx, double &fy, double &fz, double &pe) {
    if (!d_pair_fp64) {
        int pm = d_params.pair_mode[ta][tj];
        if (pm == 0) return;
        const float x = (float)dx, y = (float)dy, z = (float)dz;
        const float r2 = x*x + y*y + z*z;
        if (r2 >= d_pf_csq[ta][tj] || r2 < 1e-20f) return;
        float fmag, e;
        if (pm == 1) {
            const float eps = d_pf_eps[ta][tj];
            if (eps == 0.0f) return;
            const float r2inv = 1.0f / r2;
            const float r6inv = d_pf_sig6[ta][tj] * r2inv * r2inv * r2inv;
            fmag = 48.0f * eps * r2inv * (r6inv * r6inv - 0.5f * r6inv);
            e = 0.5f * (4.0f * eps * (r6inv * r6inv - r6inv) + eps);
        } else {   // E = A (1 + cos(pi r / rc)); F = A pi / rc sin(pi r / rc) / r
            const float r = sqrtf(r2);
            float sn, cs;
            sincospif(r * d_pf_invrc[ta][tj], &sn, &cs);
            fmag = d_pf_fpre[ta][tj] * sn / r;
            e = 0.5f * d_pf_A[ta][tj] * (1.0f + cs);
        }
        fx += (double)(fmag * x); fy += (double)(fmag * y); fz += (double)(fmag * z);
        pe += (double)e;
        return;
    }
    double r2 = dx*dx + dy*dy + dz*dz;

    int pm = d_params.pair_mode[ta][tj];
    if (pm == 0) return;
    double csq = d_params.cutoff_sq[ta][tj];
    if (r2 >= csq || r2 < 1e-20) return;

    if (pm == 1) { // WCA (shifted LJ)
        double eps = d_params.epsilon[ta][tj];
        if (eps == 0.0) return;
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
    // the loop is latency-bound on the busiest thread (a DNA atom near the membrane has ~100 boundary neighbours and
    // every neighbour cost two dependent global loads plus FP64 sqrt/sin/cos in a serial chain). Neighbours are now loaded four
    // at a time (index, type, position: independent loads in flight together) and then accumulated one by one in list order
    // with the unchanged arithmetic (pair_term), so every atom's force and energy are bit-identical to the one-at-a-time loop.
    int nn = __ldg(&kd.numneigh[a]);
    const int *nl = kd.neighlist + a;
    const size_t nls = (size_t)kd.neighlist_stride;
    int jj = 0;
    for (; jj + 4 <= nn; jj += 4) {
        int j0 = __ldg(nl + (size_t)(jj + 0) * nls);
        int j1 = __ldg(nl + (size_t)(jj + 1) * nls);
        int j2 = __ldg(nl + (size_t)(jj + 2) * nls);
        int j3 = __ldg(nl + (size_t)(jj + 3) * nls);
        int t0 = __ldg(&kd.type[j0]), t1 = __ldg(&kd.type[j1]), t2 = __ldg(&kd.type[j2]), t3 = __ldg(&kd.type[j3]);
        double x0 = __ldg(&kd.x[3*j0]), y0 = __ldg(&kd.x[3*j0+1]), z0 = __ldg(&kd.x[3*j0+2]);
        double x1 = __ldg(&kd.x[3*j1]), y1 = __ldg(&kd.x[3*j1+1]), z1 = __ldg(&kd.x[3*j1+2]);
        double x2 = __ldg(&kd.x[3*j2]), y2 = __ldg(&kd.x[3*j2+1]), z2 = __ldg(&kd.x[3*j2+2]);
        double x3 = __ldg(&kd.x[3*j3]), y3 = __ldg(&kd.x[3*j3+1]), z3 = __ldg(&kd.x[3*j3+2]);
        pair_term(ta, t0, xa - x0, ya - y0, za - z0, fx, fy, fz, pe);
        pair_term(ta, t1, xa - x1, ya - y1, za - z1, fx, fy, fz, pe);
        pair_term(ta, t2, xa - x2, ya - y2, za - z2, fx, fy, fz, pe);
        pair_term(ta, t3, xa - x3, ya - y3, za - z3, fx, fy, fz, pe);
    }
    for (; jj < nn; jj++) {
        int j = __ldg(nl + (size_t)jj * nls);
        int tj = __ldg(&kd.type[j]);
        pair_term(ta, tj, xa - __ldg(&kd.x[3*j]), ya - __ldg(&kd.x[3*j+1]), za - __ldg(&kd.x[3*j+2]), fx, fy, fz, pe);
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
        } else if (!d_pair_fp64) {   // harmonic bond in FP32, FP64 accumulation
            const float x = (float)dx, y = (float)dy, z = (float)dz;
            const float r = sqrtf(x*x + y*y + z*z), K = (float)d_params.bond_K[bt], dr = r - (float)d_params.bond_r0[bt];
            const float fscale = -2.0f * K * dr / r;
            fx += (double)(fscale * x); fy += (double)(fscale * y); fz += (double)(fscale * z);
            pe += (double)(0.5f * K * dr * dr);
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

        if (!d_pair_fp64) {   // cosine angle in FP32 on FP64 differences, FP64 accumulation
            const float e1x = (float)(__ldg(&kd.x[3*a_i])   - __ldg(&kd.x[3*a_j]));
            const float e1y = (float)(__ldg(&kd.x[3*a_i+1]) - __ldg(&kd.x[3*a_j+1]));
            const float e1z = (float)(__ldg(&kd.x[3*a_i+2]) - __ldg(&kd.x[3*a_j+2]));
            const float e2x = (float)(__ldg(&kd.x[3*a_k])   - __ldg(&kd.x[3*a_j]));
            const float e2y = (float)(__ldg(&kd.x[3*a_k+1]) - __ldg(&kd.x[3*a_j+1]));
            const float e2z = (float)(__ldg(&kd.x[3*a_k+2]) - __ldg(&kd.x[3*a_j+2]));
            const float q1 = e1x*e1x + e1y*e1y + e1z*e1z, q2 = e2x*e2x + e2y*e2y + e2z*e2z;
            const float s1 = sqrtf(q1), s2 = sqrtf(q2);
            if (s1 < 1e-20f || s2 < 1e-20f) continue;
            const float inv12 = 1.0f / (s1 * s2);
            float cf = (e1x*e2x + e1y*e2y + e1z*e2z) * inv12;
            cf = fmaxf(-1.0f, fminf(1.0f, cf));
            const float Kf = (float)d_params.angle_K[atype];
            pe += (double)(Kf * (1.0f + cf) / 3.0f);
            const float cq1 = cf / q1, cq2 = cf / q2;
            const float gix = -Kf * (e2x * inv12 - cq1 * e1x), giy = -Kf * (e2y * inv12 - cq1 * e1y), giz = -Kf * (e2z * inv12 - cq1 * e1z);
            const float gkx = -Kf * (e1x * inv12 - cq2 * e2x), gky = -Kf * (e1y * inv12 - cq2 * e2y), gkz = -Kf * (e1z * inv12 - cq2 * e2z);
            if (a == a_i)      { fx += (double)gix; fy += (double)giy; fz += (double)giz; }
            else if (a == a_k) { fx += (double)gkx; fy += (double)gky; fz += (double)gkz; }
            else               { fx -= (double)(gix+gkx); fy -= (double)(giy+gky); fz -= (double)(giz+gkz); }
            continue;
        }
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

// K-value grid reduction behind ONE grid barrier.
// Each block writes its K block partials (block_reduce_sum / block_reduce_max, as before) to buf[k*gridDim.x + blockIdx.x], then
// one grid.sync, then EVERY block sums all partials itself in the order block 0 used before (thread t takes b = t, t+blockDim, ...,
// then block_reduce), so every block holds the value the two-barrier version broadcast through `scalars`, bit for bit. The caller
// alternates two buffers (parity rp): a block may rewrite buffer p only after the next reduction's barrier, which every block
// reaches only after it has finished reading p.
template <int K>
__device__ __forceinline__
void grid_reduce_k(cg::grid_group &grid, const double (&v)[K], const bool (&is_max)[K],
                   double *sdata, KernelData &kd, int &rp, double (&out)[K]) {
    double *buf = rp ? kd.block_sums_b : kd.block_sums_a;
    rp ^= 1;
    for (int k = 0; k < K; k++) {
        double b = is_max[k] ? block_reduce_max(v[k], sdata) : block_reduce_sum(v[k], sdata);
        if (threadIdx.x == 0) buf[k * gridDim.x + blockIdx.x] = b;
        __syncthreads();
    }
    grid.sync();
    for (int k = 0; k < K; k++) {
        double s;
        if (is_max[k]) {
            s = -1e300;
            for (int b = threadIdx.x; b < (int)gridDim.x; b += blockDim.x) s = fmax(s, buf[k * gridDim.x + b]);
            s = block_reduce_max(s, sdata);
        } else {
            s = 0.0;
            for (int b = threadIdx.x; b < (int)gridDim.x; b += blockDim.x) s += buf[k * gridDim.x + b];
            s = block_reduce_sum(s, sdata);
        }
        out[k] = s;
        __syncthreads();
    }
}

// Full energy+force evaluation plus the dot products the caller needs from the new forces, in one grid reduction:
//   out[0] = potential energy, out[1] = f.h, out[2] = f.f, out[3] = f.g
// The element sums run in the same element order (i = tid, tid + nthreads, ...) as the separate loops they replace, so each value
// is identical to the one the separate fh / fdotf+fdotg / initial-norm reductions produced. The barrier after the force loop stays:
// the element loop reads forces written by other threads.
__device__
void energy_force_dots(KernelData &kd, cg::grid_group &grid, double *sdata, int tid, int nthreads, int N3,
                       int &rp, double (&out)[4]) {
    double pe_local = 0.0;
    for (int k = tid; k < kd.N; k += nthreads) {
        int a = kd.order[k];
        double fx, fy, fz;
        pe_local += compute_atom_forces(kd, a, &fx, &fy, &fz);
        kd.f[3*a]   = fx;
        kd.f[3*a+1] = fy;
        kd.f[3*a+2] = fz;
    }
    grid.sync();
    double fh_l = 0.0, ff_l = 0.0, fg_l = 0.0;
    for (int i = tid; i < N3; i += nthreads) {
        double fi = kd.f[i];
        fh_l += fi * kd.h[i];
        ff_l += fi * fi;
        fg_l += fi * kd.g[i];
    }
    const double v[4] = {pe_local, fh_l, ff_l, fg_l};
    const bool mx[4] = {false, false, false, false};
    grid_reduce_k<4>(grid, v, mx, sdata, kd, rp, out);
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

// same iteration as before (LAMMPS min_cg_kokkos + linemin_quadratic), with fewer grid barriers per iteration
// (typically 4 instead of ~14): every reduction is one barrier (grid_reduce_k); pe, f.h, f.f and f.g come out of each force
// evaluation together; the next iteration's f.h and max|h| are taken from the direction update (after g = f, f.h_new == g.h_new
// with the same products in the same order, and after a steepest-descent reset f.h == f.f); the displacement check rides on the
// direction-update reduction; the barrier after saving x0 is dropped (x0 is read back only by the thread that wrote it).
// 2 blocks/SM as before (grid of 296 on B200 unchanged), so the unrolled pair loop cannot lower occupancy
__global__ void __launch_bounds__(256, 2) fused_cg_kernel(KernelData kd) {
    extern __shared__ double sdata[];
    auto grid = cg::this_grid();
    int tid = grid.thread_rank();
    int nthreads = grid.num_threads();
    int N3 = kd.N * 3;
    double inv_N = 1.0 / (double)kd.N;
    int rp = 0;              // reduction buffer parity (grid-uniform)
    double ev[4];

    int neval = kd.start_neval;
    int niter = 0;
    int stop = STOP_RUNNING;

    // nlimit for periodic CG restart (LAMMPS: min(MAXSMALLINT, ndoftotal))
    int nlimit = (N3 < 2000000000) ? N3 : 2000000000;

    // ---- Initial energy & forces (g, h are only read for the unused f.h / f.g on a first launch) ----
    energy_force_dots(kd, grid, sdata, tid, nthreads, N3, rp, ev);
    double ecurrent = ev[0];
    neval++;

    if (tid == 0) *kd.d_energy_init = ecurrent;

    // Initial force norm (LAMMPS: TWO norm => fnorm_sqr = Σf²)
    double gnorm2 = ev[2];

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

    bool carry = false;               // next iteration's f.h and max|h| already known from the direction update
    double carry_fdoth = 0.0, carry_hmax = 0.0;

    // ---- Main CG loop (matches min_cg_kokkos::iterate) ----
    for (niter = 0; niter < d_params.maxiter && stop == STOP_RUNNING; niter++) {
        double eprevious = ecurrent;

        // ================================================================
        //  linemin_quadratic  (matches min_linesearch_kokkos.cpp)
        // ================================================================

        double fdothall, hmaxall;
        if (carry) {
            fdothall = carry_fdoth;
            hmaxall = carry_hmax;
        } else {
            // Fused fdothall (sum) + hmaxall (max) in one pass over h[]/f[]
            double fdoth_l = 0.0, hmax_l = 0.0;
            for (int i = tid; i < N3; i += nthreads) {
                fdoth_l += kd.f[i] * kd.h[i];
                hmax_l = fmax(hmax_l, fabs(kd.h[i]));
            }
            const double v[2] = {fdoth_l, hmax_l};
            const bool mx[2] = {false, true};
            double o[2];
            grid_reduce_k<2>(grid, v, mx, sdata, kd, rp, o);
            fdothall = o[0]; hmaxall = o[1];
        }
        fdothall *= inv_N;

        if (fdothall <= 0.0) { stop = STOP_ETOL; break; }
        if (hmaxall == 0.0) { stop = STOP_FTOL; break; }

        double alphamax = fmin(LM_ALPHA_MAX, d_params.dmax / hmaxall);

        // Save x0 (read back only by the writing thread: no barrier)
        for (int i = tid; i < N3; i += nthreads) kd.x0[i] = kd.x[i];

        // Initialize line search state
        double alpha = alphamax;
        double fhprev = fdothall;
        double engprev = eprevious;
        double alphaprev = 0.0;
        int ls_fail = 0;
        double ff_last = 0.0, fg_last = 0.0;   // f.f and f.g of the last force evaluation

        // Backtracking loop (matches LAMMPS while(true) in linemin_quadratic)
        for (int ls = 0; ls < 50 && !ls_fail; ls++) {
            // alpha_step: reset to x0, then step to x0 + alpha*h
            for (int i = tid; i < N3; i += nthreads)
                kd.x[i] = kd.x0[i] + alpha * kd.h[i];
            grid.sync();
            energy_force_dots(kd, grid, sdata, tid, nthreads, N3, rp, ev);
            ecurrent = ev[0]; ff_last = ev[2]; fg_last = ev[3];
            neval++;

            // fh = f · h (normalized)
            double fh = ev[1];
            fh *= inv_N;

            double delfh = fh - fhprev;

            // ZEROQUAD check
            if (fabs(fh) < LM_EPS_QUAD || fabs(delfh) < LM_EPS_QUAD) {
                for (int i = tid; i < N3; i += nthreads) kd.x[i] = kd.x0[i];
                grid.sync();
                energy_force_dots(kd, grid, sdata, tid, nthreads, N3, rp, ev);
                ecurrent = ev[0]; ff_last = ev[2]; fg_last = ev[3];
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
                energy_force_dots(kd, grid, sdata, tid, nthreads, N3, rp, ev);
                ecurrent = ev[0]; ff_last = ev[2]; fg_last = ev[3];
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
                energy_force_dots(kd, grid, sdata, tid, nthreads, N3, rp, ev);
                ecurrent = ev[0]; ff_last = ev[2]; fg_last = ev[3];
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

        // Force tolerance: LAMMPS TWO norm => Σf² < ftol²  (f.f and f.g of the last evaluation)
        gnorm2 = ff_last;
        double fdotg = fg_last;

        if (gnorm2 < d_params.ftol * d_params.ftol) { stop = STOP_FTOL; break; }

        // ================================================================
        //  CG direction update (matches min_cg_kokkos Polak-Ribiere)
        // ================================================================

        double beta = (gg > 1e-30) ? fmax(0.0, (gnorm2 - fdotg) / gg) : 0.0;

        // Periodic CG restart every nlimit iterations
        if (((kd.start_iter + niter + 1) % nlimit) == 0) beta = 0.0;

        gg = gnorm2;

        // Direction update + g.h (downhill check, and next iteration's f.h since g == f) + max|h| + max|g| in one pass;
        // every 4th iteration also the displacement since the last neighbour build.
        double gdoth_l = 0.0, hmh_l = 0.0, hmg_l = 0.0, dmax_l = 0.0;
        for (int i = tid; i < N3; i += nthreads) {
            kd.g[i] = kd.f[i];
            kd.h[i] = kd.g[i] + beta * kd.h[i];
            gdoth_l += kd.g[i] * kd.h[i];
            hmh_l = fmax(hmh_l, fabs(kd.h[i]));
            hmg_l = fmax(hmg_l, fabs(kd.g[i]));
        }
        bool disp_check = ((niter & 3) == 3);
        if (disp_check) {
            for (int i = tid; i < kd.N; i += nthreads) {
                double dx = kd.x[3*i]   - kd.x_build[3*i];
                double dy = kd.x[3*i+1] - kd.x_build[3*i+1];
                double dz = kd.x[3*i+2] - kd.x_build[3*i+2];
                dmax_l = fmax(dmax_l, dx*dx + dy*dy + dz*dz);
            }
        }
        const double v[4] = {gdoth_l, hmh_l, hmg_l, dmax_l};
        const bool mx[4] = {false, true, true, true};
        double o[4];
        grid_reduce_k<4>(grid, v, mx, sdata, kd, rp, o);
        double gdoth = o[0];

        // Check if new direction is downhill; reset to steepest descent if not
        if (gdoth <= 0.0) {
            for (int i = tid; i < N3; i += nthreads) kd.h[i] = kd.g[i];   // element-owned: no barrier
            carry_fdoth = gnorm2;     // f.h with h = g = f is f.f, same products and order
            carry_hmax = o[2];
        } else {
            carry_fdoth = gdoth;
            carry_hmax = o[1];
        }
        carry = true;

        // Displacement check for neighbor rebuild
        if (disp_check) {
            double max_disp2 = o[3];
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
//  split CG, the same iteration as fused_cg_kernel run as a chain of ordinary kernels
// ============================================================================
// The persistent cooperative kernel is register-bound (128 regs, 2 blocks/SM): its force loop, latency-bound on dependent
// gathers, took ~86 us per evaluation with 16 warps per SM (probe). Here each "tick" is three launches, enqueued back to back
// without host synchronisation:
//   w_k_vec   (elements, 256 threads/block): the vector work of the current mode, plus block partials of its reductions
//   w_k_force (atoms, one per thread, high occupancy): compute_atom_forces (unchanged) and the atom's own f.h, f.f, f.g terms
//   w_k_ctrl  (one block): reduces the partials in a fixed order and runs the line search / CG logic of fused_cg_kernel
//             (same tests, same constants, same order), which sets the next tick's mode. Inactive kernels return at once.
// The host enqueues ticks in small batches and polls a host-mapped done flag; at the end the controller writes the same outputs
// as fused_cg_kernel (stop code, iterations, evaluations, energies, |f|, gg), so the rebuild/relaunch loop is unchanged.
// Reductions are partial sums in a different order than grid_reduce_k . WCM_CG_PERSISTENT=1 keeps the old kernel.
enum { WM_DONE = 0, WM_EVAL = 1, WM_FH = 2, WM_UPD = 3, WM_RESET = 4 };
enum { WP_INIT = 0, WP_FH = 1, WP_LS_ALPHA = 2, WP_LS_SECANT = 3, WP_FAIL = 4, WP_UPD = 5, WP_RESET_STOP = 6 };

struct WCGState {
    int mode, phase;
    int init_gh, set_xbuild, reset_h, save_x0, copy_x0, disp_check;
    double a, beta;
    double ecurrent, eprevious, gnorm2, gg, fdothall, alpha, alphamax, fhprev, engprev, alphaprev, fh, ff_last, fg_last;
    double carry_fdoth, carry_hmax;
    int ls, niter, neval, stop;
};

__device__ __forceinline__ double w_warp_sum(double v) { for (int o = 16; o > 0; o >>= 1) v += __shfl_down_sync(0xffffffffu, v, o); return v; }
__device__ __forceinline__ double w_warp_max(double v) { for (int o = 16; o > 0; o >>= 1) v = fmax(v, __shfl_down_sync(0xffffffffu, v, o)); return v; }

// block reduction of up to 4 values (sum or max per value), result in thread 0; deterministic order
template <int K>
__device__ __forceinline__ void w_block_reduce(double (&v)[K], const bool (&is_max)[K], double *sh) {
    const int lane = threadIdx.x & 31, wid = threadIdx.x >> 5, nw = (blockDim.x + 31) >> 5;
    #pragma unroll
    for (int k = 0; k < K; k++) v[k] = is_max[k] ? w_warp_max(v[k]) : w_warp_sum(v[k]);
    __syncthreads();
    if (lane == 0) {
        #pragma unroll
        for (int k = 0; k < K; k++) sh[k * 32 + wid] = v[k];
    }
    __syncthreads();
    if (wid == 0) {
        #pragma unroll
        for (int k = 0; k < K; k++) {
            double x = (lane < nw) ? sh[k * 32 + lane] : (is_max[k] ? 0.0 : 0.0);
            v[k] = is_max[k] ? w_warp_max(x) : w_warp_sum(x);
        }
    }
}

__global__ void __launch_bounds__(256) w_k_vec(KernelData kd, const WCGState *st, double *part, int P) {
    wcm_grid_dependency_sync();
    __shared__ double sh[4 * 32];
    const int mode = st->mode;
    if (mode == WM_DONE) return;
    const int N3 = 3 * kd.N;
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    double v[4] = {0.0, 0.0, 0.0, 0.0};
    if (i < N3) {
        if (st->init_gh) { kd.g[i] = kd.f[i]; kd.h[i] = kd.f[i]; }
        if (st->set_xbuild) kd.x_build[i] = kd.x[i];
        if (mode == WM_FH) {
            v[0] = kd.f[i] * kd.h[i];
            v[1] = fabs(kd.h[i]);
        } else if (mode == WM_UPD) {
            kd.g[i] = kd.f[i];
            kd.h[i] = kd.g[i] + st->beta * kd.h[i];
            v[0] = kd.g[i] * kd.h[i];
            v[1] = fabs(kd.h[i]);
            v[2] = fabs(kd.g[i]);
            if (st->disp_check && (i % 3) == 0) {
                double dx = kd.x[i] - kd.x_build[i], dy = kd.x[i+1] - kd.x_build[i+1], dz = kd.x[i+2] - kd.x_build[i+2];
                v[3] = dx*dx + dy*dy + dz*dz;
            }
        } else if (mode == WM_EVAL) {
            if (st->reset_h) kd.h[i] = kd.g[i];
            if (st->save_x0) kd.x0[i] = kd.x[i];
            if (st->copy_x0) kd.x[i] = kd.x0[i];
            else if (st->a != 0.0 || st->save_x0) kd.x[i] = kd.x0[i] + st->a * kd.h[i];
        } else if (mode == WM_RESET) {
            kd.h[i] = kd.g[i];
        }
    }
    if (mode == WM_FH || mode == WM_UPD) {
        const bool mx[4] = {false, true, true, true};
        w_block_reduce<4>(v, mx, sh);
        if (threadIdx.x == 0) for (int k = 0; k < 4; k++) part[k * P + blockIdx.x] = v[k];
    }
}

// w_k_vec with one ATOM (three elements) per thread: a third of the blocks and of the partials the controller sums;
// the per-element arithmetic is unchanged, the partial sums run in a different fixed order. WCM_CG_ATOM_VECTOR_OFF=1 keeps w_k_vec.
__global__ void __launch_bounds__(256) w_k_vec3(KernelData kd, const WCGState *st, double *part, int P) {
    wcm_grid_dependency_sync();
    __shared__ double sh[4 * 32];
    const int mode = st->mode;
    if (mode == WM_DONE) return;
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    double v[4] = {0.0, 0.0, 0.0, 0.0};
    if (a < kd.N) {
        const int init_gh = st->init_gh, set_xbuild = st->set_xbuild;
        if (mode == WM_FH) {
            #pragma unroll
            for (int c = 0; c < 3; c++) {
                const int i = 3 * a + c;
                if (init_gh) { kd.g[i] = kd.f[i]; kd.h[i] = kd.f[i]; }
                if (set_xbuild) kd.x_build[i] = kd.x[i];
                v[0] += kd.f[i] * kd.h[i];
                v[1] = fmax(v[1], fabs(kd.h[i]));
            }
        } else if (mode == WM_UPD) {
            const double beta = st->beta;
            #pragma unroll
            for (int c = 0; c < 3; c++) {
                const int i = 3 * a + c;
                if (init_gh) { kd.g[i] = kd.f[i]; kd.h[i] = kd.f[i]; }
                if (set_xbuild) kd.x_build[i] = kd.x[i];
                const double g = kd.f[i], h = g + beta * kd.h[i];
                kd.g[i] = g; kd.h[i] = h;
                v[0] += g * h;
                v[1] = fmax(v[1], fabs(h));
                v[2] = fmax(v[2], fabs(g));
            }
            if (st->disp_check) {
                const int i = 3 * a;
                double dx = kd.x[i] - kd.x_build[i], dy = kd.x[i+1] - kd.x_build[i+1], dz = kd.x[i+2] - kd.x_build[i+2];
                v[3] = dx*dx + dy*dy + dz*dz;
            }
        } else {
            const int reset_h = st->reset_h, save_x0 = st->save_x0, copy_x0 = st->copy_x0;
            const double sa = st->a;
            #pragma unroll
            for (int c = 0; c < 3; c++) {
                const int i = 3 * a + c;
                if (init_gh) { kd.g[i] = kd.f[i]; kd.h[i] = kd.f[i]; }
                if (set_xbuild) kd.x_build[i] = kd.x[i];
                if (mode == WM_EVAL) {
                    if (reset_h) kd.h[i] = kd.g[i];
                    if (save_x0) kd.x0[i] = kd.x[i];
                    if (copy_x0) kd.x[i] = kd.x0[i];
                    else if (sa != 0.0 || save_x0) kd.x[i] = kd.x0[i] + sa * kd.h[i];
                } else if (mode == WM_RESET) {
                    kd.h[i] = kd.g[i];
                }
            }
        }
    }
    if (mode == WM_FH || mode == WM_UPD) {
        const bool mx[4] = {false, true, true, true};
        w_block_reduce<4>(v, mx, sh);
        if (threadIdx.x == 0) for (int k = 0; k < 4; k++) part[k * P + blockIdx.x] = v[k];
    }
}

__global__ void __launch_bounds__(256) w_k_force(KernelData kd, const WCGState *st, double *part, int P) {
    wcm_grid_dependency_sync();
    __shared__ double sh[4 * 32];
    if (st->mode != WM_EVAL) return;
    // atoms in the work-ranked chunk order of 045 (kd.order: the heaviest 32-atom chunks in the first slots), so the
    // warps with the longest neighbour lists (DNA near the membrane) start in the first wave instead of setting a tail
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    const int a = (k < kd.N) ? kd.order[k] : 0;
    double v[4] = {0.0, 0.0, 0.0, 0.0};
    if (k < kd.N) {
        double fx, fy, fz;
        v[0] = compute_atom_forces(kd, a, &fx, &fy, &fz);
        kd.f[3*a] = fx; kd.f[3*a+1] = fy; kd.f[3*a+2] = fz;
        const double hx = kd.h[3*a], hy = kd.h[3*a+1], hz = kd.h[3*a+2];
        const double gx = kd.g[3*a], gy = kd.g[3*a+1], gz = kd.g[3*a+2];
        v[1] = fx*hx + fy*hy + fz*hz;
        v[2] = fx*fx + fy*fy + fz*fz;
        v[3] = fx*gx + fy*gy + fz*gz;
    }
    const bool mx[4] = {false, false, false, false};
    w_block_reduce<4>(v, mx, sh);
    if (threadIdx.x == 0) for (int k = 0; k < 4; k++) part[k * P + blockIdx.x] = v[k];
}

// ============================================================================================================================
// w_k_force with K lanes per atom. The pair terms of an atom's neighbours are computed K at a time: lane l of the atom's
// group takes neighbours jj = l, l+K, l+2K, ... (four per lane in flight, as compute_atom_forces loads them), each lane computes
// the FP32 force/energy terms of its pairs (the mixed-fp32 change arithmetic, same expressions), and lane 0 receives them by
// shuffle and adds them to its FP64 accumulators IN LIST ORDER, skipping pairs pair_term would skip. Bonds and then angles are
// split the same way (entry b_start + l, ...; FP32 terms of the fp32-bonded change branches, added by lane 0 in entry order), and lane 0
// adds the FP64 ribosome spring last. So every atom's force and energy are the same FP64 sums in the same order as w_k_force's
// (the FP64 sum of FP32 terms is NOT order-independent: a lane-partial version changed 6 bytes of the t2996 monomers). A block
// still owns the same 256 slots of kd.order (in R = 256K/T rounds of T/K atoms); each atom's pe, f.h, f.f, f.g go to shared
// memory by slot and the first 256 threads reduce them with exactly the tree of w_block_reduce<4> over 256 threads; the three
// dot products are written with the contraction w_k_force's build has (compiler choice differed otherwise), so the block
// partials are bit for bit those of w_k_force. Used only in the mixed-precision mode with harmonic bonds (FP64 pair mode or FENE
// bonds: w_k_force). WCM_CG_FORCE_LANES_OFF=1: w_k_force.
// ============================================================================================================================
// the FP32 terms pair_term (the mixed-fp32 change branch) adds for one pair; false when it adds nothing
__device__ __forceinline__ bool w_cg_force_lanes_pair_f32(int ta, int tj, double dx, double dy, double dz, float &tx, float &ty, float &tz, float &te) {
    int pm = d_params.pair_mode[ta][tj];
    if (pm == 0) return false;
    const float x = (float)dx, y = (float)dy, z = (float)dz;
    const float r2 = x*x + y*y + z*z;
    if (r2 >= d_pf_csq[ta][tj] || r2 < 1e-20f) return false;
    float fmag, e;
    if (pm == 1) {
        const float eps = d_pf_eps[ta][tj];
        if (eps == 0.0f) return false;
        const float r2inv = 1.0f / r2;
        const float r6inv = d_pf_sig6[ta][tj] * r2inv * r2inv * r2inv;
        fmag = 48.0f * eps * r2inv * (r6inv * r6inv - 0.5f * r6inv);
        e = 0.5f * (4.0f * eps * (r6inv * r6inv - r6inv) + eps);
    } else {
        const float r = sqrtf(r2);
        float sn, cs;
        sincospif(r * d_pf_invrc[ta][tj], &sn, &cs);
        fmag = d_pf_fpre[ta][tj] * sn / r;
        e = 0.5f * d_pf_A[ta][tj] * (1.0f + cs);
    }
    tx = fmag * x; ty = fmag * y; tz = fmag * z; te = e;
    return true;
}

// the ribosome spring of compute_atom_forces (FP64), same text
__device__ __forceinline__ void w_cg_force_lanes_spring(const KernelData &kd, int a, double xa, double ya, double za,
                                            double &fx, double &fy, double &fz, double &pe) {
        double dx = xa - __ldg(&kd.xref[3*a]);
        double dy = ya - __ldg(&kd.xref[3*a+1]);
        double dz = za - __ldg(&kd.xref[3*a+2]);
        double Ks = d_params.spring_K_ribo;
        fx -= Ks * dx;
        fy -= Ks * dy;
        fz -= Ks * dz;
        pe += 0.5 * Ks * (dx*dx + dy*dy + dz*dz);
}

template <int K, int T, int U = 4, int MINB = 1>
__global__ void __launch_bounds__(T, MINB) w_k_force_lanes(KernelData kd, const WCGState *st, double *part, int P) {
    wcm_grid_dependency_sync();
    __shared__ double sv[4][256];
    __shared__ double sh[4 * 32];
    if (st->mode != WM_EVAL) return;
    constexpr int APR = T / K, R = 256 / APR;
    static_assert(32 % K == 0 && 256 % APR == 0, "K lanes must tile a warp, rounds must tile 256 slots");
    const int lane = threadIdx.x % K, sub = threadIdx.x / K;
    const int gbase = (threadIdx.x & 31) & ~(K - 1);
    for (int r = 0; r < R; r++) {
        const int slot = r * APR + sub;
        const int k = blockIdx.x * 256 + slot;
        double fx = 0.0, fy = 0.0, fz = 0.0, pe = 0.0;
        int a = 0, ta = 1;
        if (k < kd.N) { a = kd.order[k]; ta = __ldg(&kd.type[a]); }
        // the neighbour loop runs warp-uniformly (to the longest list of the warp) so that every shuffle is a full-warp one;
        // lanes past their atom's list (or of boundary / empty slots) contribute nothing
        const bool live = (ta != 1);
        double xa = 0.0, ya = 0.0, za = 0.0;
        int nn = 0;
        if (live) { xa = kd.x[3*a]; ya = kd.x[3*a+1]; za = kd.x[3*a+2]; nn = __ldg(&kd.numneigh[a]); }
        const int nn_w = wcm_warp_max(nn);
        const int *nl = kd.neighlist + a;
        const size_t nls = (size_t)kd.neighlist_stride;
        for (int base = 0; base < nn_w; base += U * K) {
            float tx[U], ty[U], tz[U], te[U]; bool ok[U];
            int j[U], tj[U]; double xj[U], yj[U], zj[U];
            #pragma unroll
            for (int u = 0; u < U; u++) { const int jj = base + u * K + lane; j[u] = (jj < nn) ? __ldg(nl + (size_t)jj * nls) : a; }
            #pragma unroll
            for (int u = 0; u < U; u++) { tj[u] = __ldg(&kd.type[j[u]]); xj[u] = __ldg(&kd.x[3*j[u]]); yj[u] = __ldg(&kd.x[3*j[u]+1]); zj[u] = __ldg(&kd.x[3*j[u]+2]); }
            #pragma unroll
            for (int u = 0; u < U; u++) ok[u] = (base + u * K + lane < nn) && w_cg_force_lanes_pair_f32(ta, tj[u], xa - xj[u], ya - yj[u], za - zj[u], tx[u], ty[u], tz[u], te[u]);
            #pragma unroll
            for (int u = 0; u < U; u++) {
                if (K == 1) {
                    if (ok[u]) { fx += (double)tx[u]; fy += (double)ty[u]; fz += (double)tz[u]; pe += (double)te[u]; }
                } else {
                    const unsigned okb = __ballot_sync(0xffffffffu, ok[u]);
                    if (lane == 0 && ok[u]) { fx += (double)tx[u]; fy += (double)ty[u]; fz += (double)tz[u]; pe += (double)te[u]; }
                    #pragma unroll
                    for (int q = 1; q < K; q++) {
                        const float sx = __shfl_sync(0xffffffffu, tx[u], gbase + q), sy = __shfl_sync(0xffffffffu, ty[u], gbase + q);
                        const float sz = __shfl_sync(0xffffffffu, tz[u], gbase + q), se = __shfl_sync(0xffffffffu, te[u], gbase + q);
                        if (lane == 0 && ((okb >> (gbase + q)) & 1u)) { fx += (double)sx; fy += (double)sy; fz += (double)sz; pe += (double)se; }
                    }
                }
            }
        }
        {
            // bonds, then angles, K at a time: lane l computes entry base + l, lane 0 adds the terms in entry order
            int b_start = 0, b_end = 0, g_start = 0, g_end = 0;
            if (live) { b_start = __ldg(&kd.atom_bond_off[a]); b_end = __ldg(&kd.atom_bond_off[a + 1]);
                        g_start = __ldg(&kd.atom_angle_off[a]); g_end = __ldg(&kd.atom_angle_off[a + 1]); }
            const int nb_w = wcm_warp_max(b_end - b_start), ng_w = wcm_warp_max(g_end - g_start);
            for (int base = 0; base < nb_w; base += K) {
                const int bi = b_start + base + lane;
                bool ok = false; float tx = 0.0f, ty = 0.0f, tz = 0.0f, te = 0.0f;
                if (bi < b_end) {
                    const int bt = __ldg(&kd.atom_bond_data[2*bi]), other = __ldg(&kd.atom_bond_data[2*bi + 1]);
                    const double dx = xa - __ldg(&kd.x[3*other]), dy = ya - __ldg(&kd.x[3*other+1]), dz = za - __ldg(&kd.x[3*other+2]);
                    const double rsq = dx*dx + dy*dy + dz*dz;
                    if (rsq >= 1e-20) {
                        const float x = (float)dx, y = (float)dy, z = (float)dz;
                        const float r = sqrtf(x*x + y*y + z*z), Kb = (float)d_params.bond_K[bt], dr = r - (float)d_params.bond_r0[bt];
                        const float fscale = -2.0f * Kb * dr / r;
                        tx = fscale * x; ty = fscale * y; tz = fscale * z; te = 0.5f * Kb * dr * dr; ok = true;
                    }
                }
                const unsigned okb = __ballot_sync(0xffffffffu, ok);
                if (lane == 0 && ok) { fx += (double)tx; fy += (double)ty; fz += (double)tz; pe += (double)te; }
                #pragma unroll
                for (int q = 1; q < K; q++) {
                    const float sx = __shfl_sync(0xffffffffu, tx, gbase + q), sy = __shfl_sync(0xffffffffu, ty, gbase + q);
                    const float sz = __shfl_sync(0xffffffffu, tz, gbase + q), se = __shfl_sync(0xffffffffu, te, gbase + q);
                    if (lane == 0 && ((okb >> (gbase + q)) & 1u)) { fx += (double)sx; fy += (double)sy; fz += (double)sz; pe += (double)se; }
                }
            }
            for (int base = 0; base < ng_w; base += K) {
                const int ai = g_start + base + lane;
                bool ok = false; float tx = 0.0f, ty = 0.0f, tz = 0.0f, te = 0.0f;
                int kind = 0;   // 0: a is a_i (+gi), 1: a is a_k (+gk), 2: vertex (-(gi+gk))
                if (ai < g_end) {
                    const int atype = __ldg(&kd.atom_angle_data[4*ai]), a_i = __ldg(&kd.atom_angle_data[4*ai + 1]);
                    const int a_j = __ldg(&kd.atom_angle_data[4*ai + 2]), a_k = __ldg(&kd.atom_angle_data[4*ai + 3]);
                    const float e1x = (float)(__ldg(&kd.x[3*a_i])   - __ldg(&kd.x[3*a_j]));
                    const float e1y = (float)(__ldg(&kd.x[3*a_i+1]) - __ldg(&kd.x[3*a_j+1]));
                    const float e1z = (float)(__ldg(&kd.x[3*a_i+2]) - __ldg(&kd.x[3*a_j+2]));
                    const float e2x = (float)(__ldg(&kd.x[3*a_k])   - __ldg(&kd.x[3*a_j]));
                    const float e2y = (float)(__ldg(&kd.x[3*a_k+1]) - __ldg(&kd.x[3*a_j+1]));
                    const float e2z = (float)(__ldg(&kd.x[3*a_k+2]) - __ldg(&kd.x[3*a_j+2]));
                    const float q1 = e1x*e1x + e1y*e1y + e1z*e1z, q2 = e2x*e2x + e2y*e2y + e2z*e2z;
                    const float s1 = sqrtf(q1), s2 = sqrtf(q2);
                    if (!(s1 < 1e-20f || s2 < 1e-20f)) {
                        const float inv12 = 1.0f / (s1 * s2);
                        float cf = (e1x*e2x + e1y*e2y + e1z*e2z) * inv12;
                        cf = fmaxf(-1.0f, fminf(1.0f, cf));
                        const float Kf = (float)d_params.angle_K[atype];
                        te = Kf * (1.0f + cf) / 3.0f;
                        const float cq1 = cf / q1, cq2 = cf / q2;
                        const float gix = -Kf * (e2x * inv12 - cq1 * e1x), giy = -Kf * (e2y * inv12 - cq1 * e1y), giz = -Kf * (e2z * inv12 - cq1 * e1z);
                        const float gkx = -Kf * (e1x * inv12 - cq2 * e2x), gky = -Kf * (e1y * inv12 - cq2 * e2y), gkz = -Kf * (e1z * inv12 - cq2 * e2z);
                        if (a == a_i)      { tx = gix; ty = giy; tz = giz; kind = 0; }
                        else if (a == a_k) { tx = gkx; ty = gky; tz = gkz; kind = 1; }
                        else               { tx = gix+gkx; ty = giy+gky; tz = giz+gkz; kind = 2; }
                        ok = true;
                    }
                }
                const unsigned okb = __ballot_sync(0xffffffffu, ok), subb = __ballot_sync(0xffffffffu, kind == 2);
                if (lane == 0 && ok) {
                    pe += (double)te;
                    if (kind != 2) { fx += (double)tx; fy += (double)ty; fz += (double)tz; }
                    else           { fx -= (double)tx; fy -= (double)ty; fz -= (double)tz; }
                }
                #pragma unroll
                for (int q = 1; q < K; q++) {
                    const float sx = __shfl_sync(0xffffffffu, tx, gbase + q), sy = __shfl_sync(0xffffffffu, ty, gbase + q);
                    const float sz = __shfl_sync(0xffffffffu, tz, gbase + q), se = __shfl_sync(0xffffffffu, te, gbase + q);
                    if (lane == 0 && ((okb >> (gbase + q)) & 1u)) {
                        pe += (double)se;
                        if (!((subb >> (gbase + q)) & 1u)) { fx += (double)sx; fy += (double)sy; fz += (double)sz; }
                        else                               { fx -= (double)sx; fy -= (double)sy; fz -= (double)sz; }
                    }
                }
            }
            if (lane == 0 && ta == 2) w_cg_force_lanes_spring(kd, a, xa, ya, za, fx, fy, fz, pe);   // ribosome spring (FP64), last
        }
        if (lane == 0) {
            double v0 = 0.0, v1 = 0.0, v2 = 0.0, v3 = 0.0;
            if (k < kd.N) {
                kd.f[3*a] = fx; kd.f[3*a+1] = fy; kd.f[3*a+2] = fz;
                const double hx = kd.h[3*a], hy = kd.h[3*a+1], hz = kd.h[3*a+2];
                const double gx = kd.g[3*a], gy = kd.g[3*a+1], gz = kd.g[3*a+2];
                v0 = pe;
                // the contraction w_k_force's build gives these three dot products (SASS: y product, x fused onto it, then z);
                // written out so that the compiler's contraction choice here cannot differ
                v1 = __fma_rn(fz, hz, __fma_rn(fx, hx, __dmul_rn(fy, hy)));
                v2 = __fma_rn(fz, fz, __fma_rn(fx, fx, __dmul_rn(fy, fy)));
                v3 = __fma_rn(fz, gz, __fma_rn(fx, gx, __dmul_rn(fy, gy)));
            }
            sv[0][slot] = v0; sv[1][slot] = v1; sv[2][slot] = v2; sv[3][slot] = v3;
        }
    }
    __syncthreads();
    // w_block_reduce<4> (all sums) as a 256-thread block computes it, on the slot values
    const int t = threadIdx.x, wl = t & 31, wid = t >> 5;
    double v[4];
    #pragma unroll
    for (int q = 0; q < 4; q++) v[q] = w_warp_sum(t < 256 ? sv[q][t] : 0.0);
    if (wl == 0 && wid < 8) {
        #pragma unroll
        for (int q = 0; q < 4; q++) sh[q * 32 + wid] = v[q];
    }
    __syncthreads();
    if (wid == 0) {
        #pragma unroll
        for (int q = 0; q < 4; q++) v[q] = w_warp_sum(wl < 8 ? sh[q * 32 + wl] : 0.0);
        if (wl == 0) for (int q = 0; q < 4; q++) part[q * P + blockIdx.x] = v[q];
    }
}

static bool force_lanes_ok = false;

   // set with the parameters: mixed-precision pair mode and harmonic bonds only

__device__ void w_finish(WCGState &s, const KernelData &kd, double inv_N, volatile int *done) {
    if (s.stop == STOP_RUNNING) s.stop = STOP_MAXITER;
    s.mode = WM_DONE;
    *kd.d_stop = s.stop;
    *kd.d_niter = kd.start_iter + s.niter;
    *kd.d_neval = s.neval;
    *kd.d_energy_final = s.ecurrent;
    *kd.d_fnorm = sqrt(s.gnorm2 * inv_N);
    *kd.d_gg = s.gg;
    __threadfence_system();
    *done = 1;
}

__device__ void w_begin_iter(WCGState &s, const KernelData &kd, double fd_raw, double hm, double inv_N, volatile int *done) {
    s.eprevious = s.ecurrent;
    s.fdothall = fd_raw * inv_N;
    if (s.fdothall <= 0.0) { s.stop = STOP_ETOL; w_finish(s, kd, inv_N, done); return; }
    if (hm == 0.0) { s.stop = STOP_FTOL; w_finish(s, kd, inv_N, done); return; }
    s.alphamax = fmin(LM_ALPHA_MAX, d_params.dmax / hm);
    s.alpha = s.alphamax; s.fhprev = s.fdothall; s.engprev = s.eprevious; s.alphaprev = 0.0; s.ls = 0;
    s.mode = WM_EVAL; s.save_x0 = 1; s.copy_x0 = 0; s.a = s.alpha; s.phase = WP_LS_ALPHA;
}

__device__ void w_post_ls(WCGState &s, const KernelData &kd, double inv_N, int nlimit, volatile int *done) {
    if (s.neval >= d_params.maxeval) { s.stop = STOP_MAXEVAL; w_finish(s, kd, inv_N, done); return; }
    double ediff = fabs(s.ecurrent - s.eprevious);
    double eref  = 0.5 * (fabs(s.ecurrent) + fabs(s.eprevious) + LM_EPS_ENERGY);
    if (ediff < d_params.etol * eref) { s.stop = STOP_ETOL; w_finish(s, kd, inv_N, done); return; }
    if (d_params.energy_target > 0.0 && s.ecurrent < d_params.energy_target) { s.stop = STOP_ENERGY_TARGET; w_finish(s, kd, inv_N, done); return; }
    if (d_params.early_stop_rtol > 0.0 && s.ecurrent > 0.0) {
        double rel_drop = (s.eprevious - s.ecurrent) / s.ecurrent;
        if (rel_drop < d_params.early_stop_rtol) { s.stop = STOP_ENERGY_TARGET; w_finish(s, kd, inv_N, done); return; }
    }
    s.gnorm2 = s.ff_last;
    double fdotg = s.fg_last;
    if (s.gnorm2 < d_params.ftol * d_params.ftol) { s.stop = STOP_FTOL; w_finish(s, kd, inv_N, done); return; }
    double beta = (s.gg > 1e-30) ? fmax(0.0, (s.gnorm2 - fdotg) / s.gg) : 0.0;
    if (((kd.start_iter + s.niter + 1) % nlimit) == 0) beta = 0.0;
    s.gg = s.gnorm2;
    s.beta = beta;
    s.disp_check = ((s.niter & 3) == 3) ? 1 : 0;
    s.mode = WM_UPD; s.phase = WP_UPD;
}

__device__ void w_armijo(WCGState &s, const KernelData &kd, double inv_N, int nlimit, volatile int *done) {
    double de_ideal = -LM_BACKTRACK_SLOPE * s.alpha * s.fdothall;
    double de = s.ecurrent - s.eprevious;
    if (de <= de_ideal) { w_post_ls(s, kd, inv_N, nlimit, done); return; }
    s.fhprev = s.fh; s.engprev = s.ecurrent; s.alphaprev = s.alpha; s.alpha *= LM_ALPHA_REDUCE; s.ls++;
    if (s.alpha <= 0.0 || de_ideal >= -LM_EMACH) { s.mode = WM_EVAL; s.copy_x0 = 1; s.phase = WP_FAIL; return; }
    if (s.ls < 50) { s.mode = WM_EVAL; s.copy_x0 = 0; s.a = s.alpha; s.phase = WP_LS_ALPHA; return; }
    w_post_ls(s, kd, inv_N, nlimit, done);    // line search exhausted without failing (as the 50-pass loop)
}

__global__ void __launch_bounds__(1024) w_k_ctrl(KernelData kd, WCGState *st, const double *part, int P, int nforce, int nvec, volatile int *done) {
    wcm_grid_dependency_sync();
    __shared__ double sh[4 * 32];
    const int mode = st->mode;
    if (mode == WM_DONE) return;
    double v[4] = {0.0, 0.0, 0.0, 0.0};
    bool mx[4] = {false, false, false, false};
    int nb = 0;
    if (mode == WM_EVAL) nb = nforce;
    else if (mode == WM_FH || mode == WM_UPD) { nb = nvec; mx[1] = mx[2] = mx[3] = true; }
    for (int b = threadIdx.x; b < nb; b += blockDim.x)
        for (int k = 0; k < 4; k++) v[k] = mx[k] ? fmax(v[k], part[k * P + b]) : v[k] + part[k * P + b];
    {
        const bool m4[4] = {mx[0], mx[1], mx[2], mx[3]};
        w_block_reduce<4>(v, m4, sh);
    }
    if (threadIdx.x != 0) return;
    WCGState s = *st;
    const double inv_N = 1.0 / (double)kd.N;
    const int N3 = 3 * kd.N;
    const int nlimit = (N3 < 2000000000) ? N3 : 2000000000;
    s.init_gh = 0; s.set_xbuild = 0; s.save_x0 = 0; s.copy_x0 = 0;
    const int was_reset = s.reset_h; s.reset_h = 0;
    (void)was_reset;
    switch (s.phase) {
    case WP_INIT: {
        s.ecurrent = v[0]; s.neval++;
        *kd.d_energy_init = s.ecurrent;
        s.gnorm2 = v[2];
        if (s.gnorm2 < d_params.ftol * d_params.ftol) s.stop = STOP_FTOL;
        if (!kd.resume) { s.init_gh = 1; s.gg = s.gnorm2; } else { s.gg = *kd.d_gg; }
        s.set_xbuild = 1;
        if (s.stop != STOP_RUNNING || d_params.maxiter <= 0) {
            // the original applied g = h = f and x_build = x before its loop; with no further launch they are not used
            w_finish(s, kd, inv_N, done);
        } else {
            s.mode = WM_FH; s.phase = WP_FH;
        }
        break;
    }
    case WP_FH:
        w_begin_iter(s, kd, v[0], v[1], inv_N, done);
        break;
    case WP_LS_ALPHA: {
        s.ecurrent = v[0]; s.ff_last = v[2]; s.fg_last = v[3]; s.neval++;
        s.fh = v[1] * inv_N;
        double delfh = s.fh - s.fhprev;
        if (fabs(s.fh) < LM_EPS_QUAD || fabs(delfh) < LM_EPS_QUAD) { s.mode = WM_EVAL; s.copy_x0 = 1; s.phase = WP_FAIL; break; }
        double relerr = fabs(1.0 - (0.5 * (s.alpha - s.alphaprev) * (s.fh + s.fhprev) + s.ecurrent) / s.engprev);
        double alpha0 = s.alpha - (s.alpha - s.alphaprev) * s.fh / delfh;
        if (relerr <= LM_QUADRATIC_TOL && alpha0 > 0.0 && alpha0 < s.alphamax) {
            s.mode = WM_EVAL; s.a = alpha0; s.phase = WP_LS_SECANT; break;
        }
        w_armijo(s, kd, inv_N, nlimit, done);
        break;
    }
    case WP_LS_SECANT:
        s.ecurrent = v[0]; s.ff_last = v[2]; s.fg_last = v[3]; s.neval++;
        if (s.ecurrent - s.eprevious < LM_EMACH) w_post_ls(s, kd, inv_N, nlimit, done);
        else w_armijo(s, kd, inv_N, nlimit, done);
        break;
    case WP_FAIL:
        s.ecurrent = v[0]; s.ff_last = v[2]; s.fg_last = v[3]; s.neval++;
        s.stop = STOP_ETOL;
        w_finish(s, kd, inv_N, done);
        break;
    case WP_UPD: {
        const double gdoth = v[0], hmh = v[1], hmg = v[2], dmax2 = v[3];
        int reset = 0;
        if (gdoth <= 0.0) { reset = 1; s.carry_fdoth = s.gnorm2; s.carry_hmax = hmg; }
        else              { s.carry_fdoth = gdoth; s.carry_hmax = hmh; }
        if (s.disp_check && dmax2 > kd.skin * kd.skin * 0.25) {
            s.stop = STOP_NEEDS_REBUILD;
            if (reset) { s.mode = WM_RESET; s.phase = WP_RESET_STOP; }   // h = g must be in place for the resumed launch
            else w_finish(s, kd, inv_N, done);
            break;
        }
        s.niter++;
        if (s.niter >= d_params.maxiter) { w_finish(s, kd, inv_N, done); break; }
        s.reset_h = reset;
        w_begin_iter(s, kd, s.carry_fdoth, s.carry_hmax, inv_N, done);
        break;
    }
    case WP_RESET_STOP:
        w_finish(s, kd, inv_N, done);
        break;
    }
    *st = s;
}

// programmatic dependent launch (PDL): kernels of the CG tick chain / BD step chain are launched with programmatic stream
// serialization, so the next grid is launched while the previous one drains; each such kernel starts with
// cudaGridDependencySynchronize() (waits for the previous grid's completion and memory flush; a no-op without the attribute).
// WCM_PDL_OFF=1, or a GPU below compute capability 9.0 = ordinary launches (include/wcm_arch.cuh).
static bool wcm_pdl_on() { return wcm_pdl_enabled(); }   // include/wcm_arch.cuh: compute capability >= 9.0, WCM_PDL_OFF unset
template <typename... KArgs, typename... Args>
static void wcm_launch(void (*k)(KArgs...), dim3 g, dim3 b, Args... args) {
    cudaLaunchConfig_t cfg = {};
    cfg.gridDim = g; cfg.blockDim = b; cfg.dynamicSmemBytes = 0; cfg.stream = 0;
    cudaLaunchAttribute at[1];
    at[0].id = cudaLaunchAttributeProgrammaticStreamSerialization;
    at[0].val.programmaticStreamSerializationAllowed = 1;
    cfg.attrs = at; cfg.numAttrs = wcm_pdl_on() ? 1 : 0;
    cudaLaunchKernelEx(&cfg, k, static_cast<KArgs>(args)...);
}

static void force_lanes_launch(int k, int BF, const KernelData &kd, const WCGState *st, double *part, int P) {
    switch (k) {
        case 1:  wcm_launch(w_k_force_lanes<1, 256>, dim3(BF), dim3(256), kd, st, part, P); break;
        case 20: wcm_launch(w_k_force_lanes<2, 512>, dim3(BF), dim3(512), kd, st, part, P); break;
        case 21: wcm_launch(w_k_force_lanes<2, 256>, dim3(BF), dim3(256), kd, st, part, P); break;
        case 4:  wcm_launch(w_k_force_lanes<4, 1024, 1>, dim3(BF), dim3(1024), kd, st, part, P); break;
        case 42: wcm_launch(w_k_force_lanes<4, 512>, dim3(BF), dim3(512), kd, st, part, P); break;
        default: wcm_launch(w_k_force_lanes<2, 512, 2, 2>, dim3(BF), dim3(512), kd, st, part, P); break;
    }
}

// the segment's start state written by a one-thread kernel (was a host-to-device copy before every segment)
__global__ void w_k_init_state(WCGState *st, int start_neval) {
    wcm_grid_dependency_sync();   // after the previous segment's queued ticks
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    WCGState s0;
    memset(&s0, 0, sizeof(s0));
    s0.mode = WM_EVAL; s0.phase = WP_INIT; s0.a = 0.0; s0.neval = start_neval; s0.stop = STOP_RUNNING;
    *st = s0;
}

// host: one segment (between neighbour builds), same contract as the cooperative launch of fused_cg_kernel
// fewer host round trips around the CG segments (mapped stop/iteration/evaluation words read when the done flag is
// seen, no device-wide synchronisations after a segment and after the neighbour build, the start state copied asynchronously
// from pinned memory, neighbour-build timing events read at the end). Same kernels, same inputs, same order in the stream.
// WCM_DEVICE_LOOPS_OFF=1 or WCM_DEVICE_LOOPS_CG_OFF=1: the synchronous host loop.
static bool w_device_loops_cg_on() { static const bool on = !std::getenv("WCM_DEVICE_LOOPS_OFF") && !std::getenv("WCM_DEVICE_LOOPS_CG_OFF"); return on; }
static WCGState *w_d_st = nullptr;
static double *w_d_part = nullptr; static int w_part_cap = 0;
static int *w_h_done = nullptr, *w_d_done = nullptr;

static void w_split_cg_segment(KernelData &kd) {
    const int N = kd.N, N3 = 3 * N;
    const int T = 256, BV = (N3 + T - 1) / T, BF = (N + T - 1) / T;
    const int P = std::max(BV, BF);
    if (!w_d_st) {
        CUDA_CHECK(cudaMalloc(&w_d_st, sizeof(WCGState)));
        CUDA_CHECK(cudaHostAlloc((void **)&w_h_done, sizeof(int), cudaHostAllocMapped));
        CUDA_CHECK(cudaHostGetDevicePointer((void **)&w_d_done, w_h_done, 0));
    }
    if (P > w_part_cap) {
        if (w_d_part) cudaFree(w_d_part);
        w_part_cap = P + P / 2;
        CUDA_CHECK(cudaMalloc(&w_d_part, 4 * (size_t)w_part_cap * sizeof(double)));
    }
    static WCGState *w_device_loops_h_s0 = nullptr;   // pinned start state (the previous copy has completed: its segment is done)
    if (!w_device_loops_h_s0) CUDA_CHECK(cudaHostAlloc((void **)&w_device_loops_h_s0, sizeof(WCGState), cudaHostAllocDefault));
    WCGState &s0 = *w_device_loops_h_s0;
    memset(&s0, 0, sizeof(s0));
    s0.mode = WM_EVAL; s0.phase = WP_INIT; s0.a = 0.0; s0.neval = kd.start_neval; s0.stop = STOP_RUNNING;
    *(volatile int *)w_h_done = 0;
    if (w_device_loops_cg_on()) wcm_launch(w_k_init_state, dim3(1), dim3(32), w_d_st, kd.start_neval);
    else              CUDA_CHECK(cudaMemcpy(w_d_st, &s0, sizeof(s0), cudaMemcpyHostToDevice));
    const int PP = w_part_cap;
    static cudaEvent_t ev[2];
    static bool ev_made = false;
    if (!ev_made) {
        ev_made = true;
        CUDA_CHECK(cudaEventCreateWithFlags(&ev[0], cudaEventDisableTiming));
        CUDA_CHECK(cudaEventCreateWithFlags(&ev[1], cudaEventDisableTiming));
    }
    static const int TICKS = std::getenv("WCM_DEVICE_LOOPS_TICKS") ? std::max(1, atoi(std::getenv("WCM_DEVICE_LOOPS_TICKS"))) : 3;   // probes
    static const bool off_cg_atom_vector = std::getenv("WCM_CG_ATOM_VECTOR_OFF") != nullptr;
    // lanes per atom in the force kernel. Default: 2 lanes, 512 threads, 2 loads per lane in flight, 2 blocks/SM (64
    // registers). Probes (WCM_CG_FORCE_LANES_K): 20 = 2 lanes, 4 loads/lane, 1 block/SM; 21 = 2 lanes, 256 threads; 4 = 4 lanes, 1024
    // threads; 42 = 4 lanes, 512 threads; 1 = one lane (shuffle-free control). At t500 / t2996 / t6000 the default was the
    // fastest of these. WCM_CG_FORCE_LANES_VERIFY=1: both kernels on each segment's first state, forces and block partials compared.
    static const int k_cg_force_lanes_env = std::getenv("WCM_CG_FORCE_LANES_OFF") ? 0 : (std::getenv("WCM_CG_FORCE_LANES_K") ? atoi(std::getenv("WCM_CG_FORCE_LANES_K")) : 2);
    const int k_cg_force_lanes = force_lanes_ok ? k_cg_force_lanes_env : 0;
    if (k_cg_force_lanes != 0 && std::getenv("WCM_CG_FORCE_LANES_VERIFY")) {   // same state, both force kernels: per-atom forces and block partials
        std::vector<double> fa(N3), fb(N3), pa(4 * (size_t)PP), pb(4 * (size_t)PP);
        wcm_launch(w_k_force, dim3(BF), dim3(T), kd, (const WCGState *)w_d_st, w_d_part, PP);
        CUDA_CHECK(cudaDeviceSynchronize());
        CUDA_CHECK(cudaMemcpy(fa.data(), kd.f, N3 * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(pa.data(), w_d_part, pa.size() * sizeof(double), cudaMemcpyDeviceToHost));
        force_lanes_launch(k_cg_force_lanes, BF, kd, (const WCGState *)w_d_st, w_d_part, PP);
        CUDA_CHECK(cudaDeviceSynchronize());
        CUDA_CHECK(cudaMemcpy(fb.data(), kd.f, N3 * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(pb.data(), w_d_part, pb.size() * sizeof(double), cudaMemcpyDeviceToHost));
        long df = 0, dp = 0; int first = -1;
        for (int i = 0; i < N3; i++) if (memcmp(&fa[i], &fb[i], 8)) { if (first < 0) first = i; df++; }
        for (int q = 0; q < 4; q++) for (int b = 0; b < BF; b++) if (memcmp(&pa[q * PP + b], &pb[q * PP + b], 8)) dp++;
        fprintf(stderr, "WCM_CG_FORCE_LANES_VERIFY: %ld of %d force components and %ld of %d block partials differ", df, N3, dp, 4 * BF);
        if (first >= 0) fprintf(stderr, " (first: atom %d comp %d %.17g vs %.17g)", first / 3, first % 3, fa[first], fb[first]);
        fprintf(stderr, "\n");
    }
    for (long batch = 0; ; batch++) {
        for (int t = 0; t < TICKS; t++) {
            if (off_cg_atom_vector) wcm_launch(w_k_vec, dim3(BV), dim3(T), kd, (const WCGState *)w_d_st, w_d_part, PP);
            else        wcm_launch(w_k_vec3, dim3(BF), dim3(T), kd, (const WCGState *)w_d_st, w_d_part, PP);
            if (k_cg_force_lanes == 0) wcm_launch(w_k_force, dim3(BF), dim3(T), kd, (const WCGState *)w_d_st, w_d_part, PP);
            else force_lanes_launch(k_cg_force_lanes, BF, kd, (const WCGState *)w_d_st, w_d_part, PP);
            wcm_launch(w_k_ctrl, dim3(1), dim3(1024), kd, w_d_st, (const double *)w_d_part, PP, BF, off_cg_atom_vector ? BV : BF, (volatile int *)w_d_done);
        }
        CUDA_CHECK(cudaEventRecord(ev[batch & 1]));
        if (batch > 0) CUDA_CHECK(cudaEventSynchronize(ev[(batch - 1) & 1]));
        if (*(volatile int *)w_h_done) break;
    }
    if (!w_device_loops_cg_on()) CUDA_CHECK(cudaDeviceSynchronize());   // the next build is stream-ordered after the queued ticks
    CUDA_CHECK(cudaGetLastError());
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

// csr-exclusions verification (WCM_CSR_EXCLUSIONS_VERIFY=1): the former builder, compared with the new one on every call
static ExclusionCSR build_exclusion_csr_ref(const FusedMinSystem *sys) {
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

static ExclusionCSR build_exclusion_csr(const FusedMinSystem *sys) {
    int N = sys->N;

    // Match LAMMPS "special_bonds angle yes" with default LJ weights (0 0 0):
    //   1-2: always excluded
    //   1-3: excluded ONLY if the pair is the first/last atoms of an angle
    //   1-4: always excluded
    // Pre-replication all 1-3 pairs are in angles (linear chain), so behavior
    // matches the old code. Post-replication, fork atoms have 3 bonded neighbors
    // creating 1-3 pairs NOT in any angle — those must get pair interactions.
    //
    // CSR adjacency, a sorted key vector for the angle end pairs and one sort+unique per atom instead of per-atom
    // std::set / a std::set of pairs (~1.5 M node allocations per call); every atom's excluded set (sorted, unique) is the same.
    std::vector<int> aoff(N + 1, 0);
    for (int b = 0; b < sys->nbonds; b++) { aoff[sys->bond_i[b] + 1]++; aoff[sys->bond_j[b] + 1]++; }
    for (int i = 0; i < N; i++) aoff[i + 1] += aoff[i];
    std::vector<int> adj(aoff[N]), fill(aoff.begin(), aoff.end() - 1);
    for (int b = 0; b < sys->nbonds; b++) {
        adj[fill[sys->bond_i[b]]++] = sys->bond_j[b];
        adj[fill[sys->bond_j[b]]++] = sys->bond_i[b];
    }

    std::vector<unsigned long long> angle_13;
    angle_13.reserve(2 * (size_t)sys->nangles);
    for (int a = 0; a < sys->nangles; a++) {
        const unsigned long long i = (unsigned)sys->angle_i[a], k = (unsigned)sys->angle_k[a];
        angle_13.push_back((i << 32) | k);
        angle_13.push_back((k << 32) | i);
    }
    std::sort(angle_13.begin(), angle_13.end());
    auto in_angle = [&angle_13](int i, int k) {
        return std::binary_search(angle_13.begin(), angle_13.end(), ((unsigned long long)(unsigned)i << 32) | (unsigned)k);
    };

    ExclusionCSR csr;
    csr.offsets.resize(N + 1);
    csr.offsets[0] = 0;
    csr.indices.reserve(8 * (size_t)aoff[N]);
    std::vector<int> ex;
    for (int i = 0; i < N; i++) {
        ex.clear();
        for (int p = aoff[i]; p < aoff[i + 1]; p++) {
            const int j = adj[p];
            ex.push_back(j);                                   // 1-2: always excluded
            for (int q = aoff[j]; q < aoff[j + 1]; q++) {
                const int k = adj[q];
                if (k == i) continue;
                if (in_angle(i, k))
                    ex.push_back(k);                           // 1-3: only if in an angle
                for (int r = aoff[k]; r < aoff[k + 1]; r++) {
                    const int l = adj[r];
                    if (l == j) continue;
                    ex.push_back(l);                           // 1-4: always excluded
                }
            }
        }
        std::sort(ex.begin(), ex.end());
        ex.erase(std::unique(ex.begin(), ex.end()), ex.end());
        csr.indices.insert(csr.indices.end(), ex.begin(), ex.end());
        csr.offsets[i + 1] = (int)csr.indices.size();
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
    p->skin    = 6.0;   // 3.0 -> 6.0 A: neighbour list rebuilt every ~30 CG iterations instead of ~15 (2704 relaunches per call at t = 2996 s)
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
    // boundary atoms (type 1) get an empty list. compute_atom_forces returns zero force and energy for type 1
    // before it reads the list, and nothing else reads it; each of them used to scan 27 cells of the 223 A grid.
    if (ti == 1) { numneigh[i] = 0; return; }
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
            // distance test first, exclusion binary search only for pairs inside the cutoff (same accepted set, same order)
            double ddx = xi - x[3*j], ddy = yi - x[3*j+1], ddz = zi - x[3*j+2];
            double d2 = ddx*ddx + ddy*ddy + ddz*ddz;
            if (d2 < rcut2 && d2 < d_nl_cut2[ti][type[j]] && !is_excluded(excl_list, e_start, e_end, j)) {
                if (nn < max_neighs)
                    neighlist[(size_t)nn * nlist_stride + i] = j;
                nn++;
            }
        }
    }
    numneigh[i] = (nn <= max_neighs) ? nn : max_neighs;
}

// balanced atom-to-thread assignment for the force loop. Each thread handles grid-stride slots k = tid, tid+n, ...
// (n = grid threads, N > n here: 149,638 atoms on 75,776 threads). The force loop is latency-bound on its busiest thread, and
// with the identity order a third of the threads get two DNA atoms, so two membrane-adjacent atoms (~100 boundary neighbours
// each) can land on one thread. After each neighbour build the atoms are sorted by work (numneigh + 1; boundary atoms 0: they
// return before reading their list) and dealt so that slot k < n gets the k-th heaviest and slot n + t gets the t-th lightest:
// the heaviest atoms share their thread with the lightest. Every atom's force is computed exactly as before; only the order
// of the per-thread energy partial sums changes.
static int *w_key_in = nullptr, *w_key_out = nullptr, *w_idx_in = nullptr, *w_idx_out = nullptr, *w_order = nullptr;
static int w_cap = 0;
static void *w_tmp = nullptr;
static size_t w_tmp_cap = 0;

__global__ void wcm_work_kernel(const int *type, const int *numneigh, int N, int *key, int *idx) {
    int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= N) return;
    key[a] = (type[a] == 1) ? 0 : numneigh[a] + 1;
    idx[a] = a;
}

// chunk c = atoms 32c .. 32c+31; key = heaviest atom of the chunk (numneigh + 1; boundary atoms 0)
__global__ void wcm_chunk_work_kernel(const int *type, const int *numneigh, int nchunks, int *key, int *idx) {
    int c = blockIdx.x * blockDim.x + threadIdx.x;
    if (c >= nchunks) return;
    int w = 0;
    for (int l = 0; l < 32; l++) {
        int a = 32 * c + l;
        int wa = (type[a] == 1) ? 0 : numneigh[a] + 1;
        w = max(w, wa);
    }
    key[c] = w;
    idx[c] = c;
}

// slot chunk s (slots 32s .. 32s+31) gets atom chunk sorted[s] if s < nc (the heaviest, first grid-stride pass), else
// sorted[nchunks - 1 + nc - s] (the lightest pair with the heaviest); slots past the full chunks keep their own atoms.
__global__ void wcm_chunk_deal_kernel(const int *sorted, int N, int nchunks, int nc, int *order) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= N) return;
    int s = k >> 5, lane = k & 31;
    if (s >= nchunks) { order[k] = k; return; }
    int sc = (s < nc) ? s : (nchunks - 1 + nc - s);
    order[k] = 32 * sorted[sc] + lane;
}

__global__ void wcm_deal_kernel(const int *sorted, int N, int n, int *order) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= N) return;
    order[k] = (k < n) ? sorted[k] : sorted[N - 1 + n - k];
}

static const int *wcm_balance_order(const int *d_type, const int *d_numneigh, int N, int nthreads) {
    if (N > w_cap) {
        if (w_key_in) { cudaFree(w_key_in); cudaFree(w_key_out); cudaFree(w_idx_in); cudaFree(w_idx_out); cudaFree(w_order); }
        int cap = N + N / 2;
        CUDA_CHECK(cudaMalloc(&w_key_in,  cap * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&w_key_out, cap * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&w_idx_in,  cap * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&w_idx_out, cap * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&w_order,   cap * sizeof(int)));
        w_cap = cap;
    }
    int threads = 256, blocks = (N + threads - 1) / threads;
    // deal whole warps' worth of consecutive atoms (chunks of 32) instead of single atoms, so the 32 lanes of a warp
    // read consecutive entries of every per-atom array and neighbour-list row (coalesced); chunks are ranked by their heaviest
    // atom (a warp waits for its slowest lane). The last partial chunk stays in the last slots.
    int nchunks = N / 32;
    int cblocks = (nchunks + threads - 1) / threads;
    if (nchunks > 0) {
        wcm_chunk_work_kernel<<<cblocks, threads>>>(d_type, d_numneigh, nchunks, w_key_in, w_idx_in);
        size_t bytes = 0;
        CUDA_CHECK(cub::DeviceRadixSort::SortPairsDescending(nullptr, bytes, w_key_in, w_key_out, w_idx_in, w_idx_out, nchunks, 0, 16));
        if (bytes > w_tmp_cap) {
            if (w_tmp) CUDA_CHECK(cudaFree(w_tmp));
            CUDA_CHECK(cudaMalloc(&w_tmp, bytes));
            w_tmp_cap = bytes;
        }
        CUDA_CHECK(cub::DeviceRadixSort::SortPairsDescending(w_tmp, bytes, w_key_in, w_key_out, w_idx_in, w_idx_out, nchunks, 0, 16));
    }
    wcm_chunk_deal_kernel<<<blocks, threads>>>(w_idx_out, N, nchunks, nthreads / 32, w_order);
    CUDA_CHECK(cudaGetLastError());
    return w_order;
}

// ============================================================================
//  two-grid neighbour build
// ============================================================================
// The single cell grid is sized by the largest list cutoff (boundary-DNA 217 + skin = 223 A), so each DNA atom gathered the
// ~730 atoms of 27 such cells although DNA-DNA pairs are listed only within 34 + 6 = 40 A. Now DNA atoms (type >= 3) are binned
// in a fine grid whose cell is the largest DNA-DNA list cutoff, and the other atoms (boundary, ribosomes) in the coarse grid of
// the global cutoff. A DNA atom takes DNA partners from its 27 fine cells and non-DNA partners from its 27 coarse cells; a
// ribosome takes DNA partners from a fine stencil wide enough for the ribosome-DNA cutoff. Same acceptance test (global and
// per-type-pair list cutoff, exclusions), so every atom gets the same neighbour set; the order within a list changes (DNA
// partners first), which only reorders force summation. Boundary atoms keep empty lists (024).
static double w_fine_cut = -1.0;     // largest DNA-DNA list cutoff (A), -1: no DNA-DNA pairs (old path)
static double w_ribo_dna_cut = -1.0; // largest ribosome-DNA list cutoff (A)
static int *w_fcell = nullptr, *w_fsorted = nullptr; static int w_fN = 0;
static int *w_fstart = nullptr, *w_fcount = nullptr; static int w_fcap = 0;
static void *w_scan_tmp2 = nullptr; static size_t w_scan_cap2 = 0;

__global__ void w_cell_idx_kernel(const double *x, const int *type, int N, int want_dna,
    int ncx, int ncy, int ncz, double cx, double cy, double cz, double blo_x, double blo_y, double blo_z,
    int *cell_of_atom, int *cell_count)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    int dna = type[i] >= 3;
    if (dna != want_dna) { cell_of_atom[i] = -1; return; }
    int ix = min(max((int)((x[3*i]   - blo_x) / cx), 0), ncx - 1);
    int iy = min(max((int)((x[3*i+1] - blo_y) / cy), 0), ncy - 1);
    int iz = min(max((int)((x[3*i+2] - blo_z) / cz), 0), ncz - 1);
    int c = ix * ncy * ncz + iy * ncz + iz;
    cell_of_atom[i] = c;
    atomicAdd(&cell_count[c], 1);
}

__global__ void w_scatter_kernel(const int *cell_of_atom, int N, const int *cell_start, int *cell_fill, int *sorted)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    int c = cell_of_atom[i];
    if (c < 0) return;
    sorted[cell_start[c] + atomicAdd(&cell_fill[c], 1)] = i;
}

__device__ __forceinline__
void w_scan_cells(const double *x, const int *type, int i, int ti, double xi, double yi, double zi,
                  const int *sorted, const int *cstart, int ix, int iy, int iz, int S, int ncx, int ncy, int ncz,
                  double rcut2, const int *excl_list, int e_start, int e_end,
                  int *neighlist, int max_neighs, int nlist_stride, int &nn)
{
    for (int dx = -S; dx <= S; dx++)
    for (int dy = -S; dy <= S; dy++)
    for (int dz = -S; dz <= S; dz++) {
        int jx = ix + dx, jy = iy + dy, jz = iz + dz;
        if (jx < 0 || jx >= ncx || jy < 0 || jy >= ncy || jz < 0 || jz >= ncz) continue;
        int cidx = jx * ncy * ncz + jy * ncz + jz;
        int ce = cstart[cidx + 1];
        for (int si = cstart[cidx]; si < ce; si++) {
            int j = sorted[si];
            if (j == i) continue;
            double ddx = xi - x[3*j], ddy = yi - x[3*j+1], ddz = zi - x[3*j+2];
            double d2 = ddx*ddx + ddy*ddy + ddz*ddz;
            if (d2 < rcut2 && d2 < d_nl_cut2[ti][type[j]] && !is_excluded(excl_list, e_start, e_end, j)) {
                if (nn < max_neighs) neighlist[(size_t)nn * nlist_stride + i] = j;
                nn++;
            }
        }
    }
}

__global__ void w_build_nlist2_kernel(
    const double *x, const int *type, int N, double blo_x, double blo_y, double blo_z, double rcut2,
    const int *f_sorted, const int *f_start, int nfx, int nfy, int nfz, double fcx, double fcy, double fcz, int S_ribo,
    const int *c_sorted, const int *c_start, int ncx, int ncy, int ncz, double cx, double cy, double cz,
    const int *excl_off, const int *excl_list, int *numneigh, int *neighlist, int max_neighs, int nlist_stride)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    int ti = type[i];
    if (ti == 1) { numneigh[i] = 0; return; }
    double xi = x[3*i], yi = x[3*i+1], zi = x[3*i+2];
    int e_start = excl_off[i], e_end = excl_off[i + 1];
    int nn = 0;
    // DNA partners: fine grid (stencil 1 for DNA atoms, S_ribo for ribosomes)
    {
        int ix = min(max((int)((xi - blo_x) / fcx), 0), nfx - 1);
        int iy = min(max((int)((yi - blo_y) / fcy), 0), nfy - 1);
        int iz = min(max((int)((zi - blo_z) / fcz), 0), nfz - 1);
        w_scan_cells(x, type, i, ti, xi, yi, zi, f_sorted, f_start, ix, iy, iz, (ti >= 3) ? 1 : S_ribo, nfx, nfy, nfz,
                     rcut2, excl_list, e_start, e_end, neighlist, max_neighs, nlist_stride, nn);
    }
    // boundary / ribosome partners: coarse grid
    {
        int ix = min(max((int)((xi - blo_x) / cx), 0), ncx - 1);
        int iy = min(max((int)((yi - blo_y) / cy), 0), ncy - 1);
        int iz = min(max((int)((zi - blo_z) / cz), 0), ncz - 1);
        w_scan_cells(x, type, i, ti, xi, yi, zi, c_sorted, c_start, ix, iy, iz, 1, ncx, ncy, ncz,
                     rcut2, excl_list, e_start, e_end, neighlist, max_neighs, nlist_stride, nn);
    }
    numneigh[i] = (nn <= max_neighs) ? nn : max_neighs;
}

#include "wcm_nlist_warp.cuh"

static void w_exclusive_scan(int *d_in, int *d_out, int n)
{
    size_t bytes = 0;
    CUDA_CHECK(cub::DeviceScan::ExclusiveSum(nullptr, bytes, d_in, d_out, n));
    if (bytes > w_scan_cap2) {
        if (w_scan_tmp2) CUDA_CHECK(cudaFree(w_scan_tmp2));
        CUDA_CHECK(cudaMalloc(&w_scan_tmp2, bytes));
        w_scan_cap2 = bytes;
    }
    CUDA_CHECK(cub::DeviceScan::ExclusiveSum(w_scan_tmp2, bytes, d_in, d_out, n));
}

static void w_rebuild_two_grid(
    double *d_x, int *d_type, int N, int *d_excl_off, int *d_excl_list, double rcut2,
    const double box_lo[3], const double box_hi[3],
    int ncx, int ncy, int ncz, double cx, double cy, double cz, int ncells,
    int *d_cell_of_atom, int *d_sorted_atoms, int *d_cell_start, int *d_cell_count,
    int *d_numneigh, int *d_neighlist, int max_neighs, int nlist_stride)
{
    int threads = 256, blocks_N = (N + threads - 1) / threads;
    double Lx = box_hi[0] - box_lo[0], Ly = box_hi[1] - box_lo[1], Lz = box_hi[2] - box_lo[2];
    int nfx = std::max(1, (int)(Lx / w_fine_cut)), nfy = std::max(1, (int)(Ly / w_fine_cut)), nfz = std::max(1, (int)(Lz / w_fine_cut));
    int nf = nfx * nfy * nfz;
    double fcx = Lx / nfx, fcy = Ly / nfy, fcz = Lz / nfz;
    int S_ribo = (w_ribo_dna_cut > 0.0) ? (int)ceil(w_ribo_dna_cut / std::min(fcx, std::min(fcy, fcz))) : 1;
    if (N > w_fN) {
        if (w_fcell) { cudaFree(w_fcell); cudaFree(w_fsorted); }
        w_fN = N + N / 2;
        CUDA_CHECK(cudaMalloc(&w_fcell, w_fN * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&w_fsorted, w_fN * sizeof(int)));
    }
    if (nf + 1 > w_fcap) {
        if (w_fstart) { cudaFree(w_fstart); cudaFree(w_fcount); }
        w_fcap = nf + 1;
        CUDA_CHECK(cudaMalloc(&w_fstart, w_fcap * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&w_fcount, w_fcap * sizeof(int)));
    }
    // coarse grid: boundary + ribosomes
    CUDA_CHECK(cudaMemset(d_cell_count, 0, (ncells + 1) * sizeof(int)));
    w_cell_idx_kernel<<<blocks_N, threads>>>(d_x, d_type, N, 0, ncx, ncy, ncz, cx, cy, cz, box_lo[0], box_lo[1], box_lo[2],
                                             d_cell_of_atom, d_cell_count);
    w_exclusive_scan(d_cell_count, d_cell_start, ncells + 1);
    CUDA_CHECK(cudaMemset(d_cell_count, 0, (ncells + 1) * sizeof(int)));
    w_scatter_kernel<<<blocks_N, threads>>>(d_cell_of_atom, N, d_cell_start, d_cell_count, d_sorted_atoms);
    // fine grid: DNA
    CUDA_CHECK(cudaMemset(w_fcount, 0, (nf + 1) * sizeof(int)));
    w_cell_idx_kernel<<<blocks_N, threads>>>(d_x, d_type, N, 1, nfx, nfy, nfz, fcx, fcy, fcz, box_lo[0], box_lo[1], box_lo[2],
                                             w_fcell, w_fcount);
    w_exclusive_scan(w_fcount, w_fstart, nf + 1);
    CUDA_CHECK(cudaMemset(w_fcount, 0, (nf + 1) * sizeof(int)));
    w_scatter_kernel<<<blocks_N, threads>>>(w_fcell, N, w_fstart, w_fcount, w_fsorted);
    // warp-per-atom build (same lists, same order); WCM_CG_WARP_NLIST_OFF=1: the thread-per-atom kernel
    static const bool off_cg_warp_nlist = std::getenv("WCM_CG_WARP_NLIST_OFF") != nullptr, verify_cg_warp_nlist = std::getenv("WCM_CG_WARP_NLIST_VERIFY") != nullptr;
    if (off_cg_warp_nlist || verify_cg_warp_nlist)
        w_build_nlist2_kernel<<<blocks_N, threads>>>(d_x, d_type, N, box_lo[0], box_lo[1], box_lo[2], rcut2,
            w_fsorted, w_fstart, nfx, nfy, nfz, fcx, fcy, fcz, S_ribo,
            d_sorted_atoms, d_cell_start, ncx, ncy, ncz, cx, cy, cz,
            d_excl_off, d_excl_list, d_numneigh, d_neighlist, max_neighs, nlist_stride);
    if (!off_cg_warp_nlist) {
        int *nn_out = d_numneigh, *nl_out = d_neighlist;
        static int *v_nn = nullptr, *v_nl = nullptr, *v_bad = nullptr; static size_t v_cap = 0;
        if (verify_cg_warp_nlist) {
            const size_t need = (size_t)max_neighs * nlist_stride;
            if (need > v_cap) {
                if (v_nl) { cudaFree(v_nl); cudaFree(v_nn); }
                CUDA_CHECK(cudaMalloc(&v_nl, need * sizeof(int))); CUDA_CHECK(cudaMalloc(&v_nn, (size_t)nlist_stride * sizeof(int)));
                if (!v_bad) CUDA_CHECK(cudaMalloc(&v_bad, sizeof(int)));
                v_cap = need;
            }
            nn_out = v_nn; nl_out = v_nl;
        }
        const long long wthreads = 32LL * N;
        w_build_nlist3_kernel<<<(unsigned)((wthreads + threads - 1) / threads), threads>>>(d_x, d_type, N, box_lo[0], box_lo[1], box_lo[2], rcut2,
            w_fsorted, w_fstart, nfx, nfy, nfz, fcx, fcy, fcz, S_ribo,
            d_sorted_atoms, d_cell_start, ncx, ncy, ncz, cx, cy, cz,
            d_excl_off, d_excl_list, nn_out, nl_out, max_neighs, nlist_stride);
        if (verify_cg_warp_nlist) {
            CUDA_CHECK(cudaMemset(v_bad, 0, sizeof(int)));
            w_nlist_compare_kernel<<<blocks_N, threads>>>(N, d_numneigh, d_neighlist, v_nn, v_nl, nlist_stride, v_bad);
            int bad = 0; CUDA_CHECK(cudaMemcpy(&bad, v_bad, sizeof(int), cudaMemcpyDeviceToHost));
            static long builds = 0, bad_builds = 0; builds++; if (bad) bad_builds++;
            if (bad || builds % 200 == 1) fprintf(stderr, "WCM_CG_WARP_NLIST_VERIFY: build %ld: %d of %d rows differ (%ld of %ld builds with differences)\n", builds, bad, N, bad_builds, builds);
        }
    }
    if (!w_device_loops_cg_on()) CUDA_CHECK(cudaDeviceSynchronize());   // nothing on the host reads the lists
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
    if (w_fine_cut > 0.0) {
        w_rebuild_two_grid(d_x, d_type, N, d_excl_off, d_excl_list, rcut2, box_lo, box_hi, ncx, ncy, ncz, cx, cy, cz, ncells,
                           d_cell_of_atom, d_sorted_atoms, d_cell_start, d_cell_count, d_numneigh, d_neighlist, max_neighs,
                           nlist_stride);
        return;
    }
    int threads = 256;
    int blocks_N = (N + threads - 1) / threads;

    gpu_compute_cell_idx_kernel<<<blocks_N, threads>>>(
        d_x, N, ncx, ncy, ncz, cx, cy, cz,
        box_lo[0], box_lo[1], box_lo[2], d_cell_of_atom);

    CUDA_CHECK(cudaMemset(d_cell_count, 0, (ncells + 1) * sizeof(int)));
    gpu_count_cell_atoms_kernel<<<blocks_N, threads>>>(d_cell_of_atom, N, d_cell_count);

    // cub scan with a kept temp buffer; thrust::exclusive_scan cudaMalloc'd and cudaFree'd its temp storage each rebuild
    static void *d_scan_tmp = nullptr;
    static size_t scan_tmp_cap = 0;
    size_t scan_tmp_bytes = 0;
    CUDA_CHECK(cub::DeviceScan::ExclusiveSum(nullptr, scan_tmp_bytes, d_cell_count, d_cell_start, ncells + 1));
    if (scan_tmp_bytes > scan_tmp_cap) {
        if (d_scan_tmp) CUDA_CHECK(cudaFree(d_scan_tmp));
        CUDA_CHECK(cudaMalloc(&d_scan_tmp, scan_tmp_bytes));
        scan_tmp_cap = scan_tmp_bytes;
    }
    CUDA_CHECK(cub::DeviceScan::ExclusiveSum(d_scan_tmp, scan_tmp_bytes, d_cell_count, d_cell_start, ncells + 1));

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

static FusedMinResult fused_cg_minimize_impl(
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

            // 8 partial slots per block in each of the two alternating reduction buffers (grid_reduce_k)
            CUDA_CHECK(cudaMalloc(&cache->d_block_sums_a,
                                  8 * cache->max_blocks * sizeof(double)));
            CUDA_CHECK(cudaMalloc(&cache->d_block_sums_b,
                                  8 * cache->max_blocks * sizeof(double)));
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
        if (getenv("WCM_CSR_EXCLUSIONS_VERIFY") != nullptr) {
            clock::time_point t_ref0 = clock::now();
            ExclusionCSR ref = build_exclusion_csr_ref(sys);
            const double ref_ms = fused_ms_since(t_ref0);
            clock::time_point t_new0 = clock::now();
            ExclusionCSR again = build_exclusion_csr(sys);
            const double new_ms = fused_ms_since(t_new0);
            fprintf(stderr, "WCM_CSR_EXCLUSIONS_VERIFY: exclusion CSR %s (%d entries; former builder %.1f ms, new %.1f ms)\n",
                    (ref.offsets == excl_csr.offsets && ref.indices == excl_csr.indices && again.indices == excl_csr.indices) ? "IDENTICAL" : "DIFFERENT",
                    excl_csr.offsets[N], ref_ms, new_ms);
        }

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
    {   // FP32 pair tables
        static const int fp64 = std::getenv("WCM_FP64") != nullptr;
        float csq[9][9], eps[9][9], sig6[9][9], A[9][9], invrc[9][9], fpre[9][9];
        for (int i = 0; i < 9; i++)
            for (int j = 0; j < 9; j++) {
                csq[i][j] = (float)params->cutoff_sq[i][j]; eps[i][j] = (float)params->epsilon[i][j];
                sig6[i][j] = (float)params->sigma_6[i][j]; A[i][j] = (float)params->soft_A[i][j];
                const double rc = params->soft_rc[i][j];
                invrc[i][j] = rc > 0.0 ? (float)(1.0 / rc) : 0.0f;
                fpre[i][j] = rc > 0.0 ? (float)(params->soft_A[i][j] * FUSED_PI / rc) : 0.0f;
            }
        CUDA_CHECK(cudaMemcpyToSymbol(d_pf_csq, csq, sizeof(csq))); CUDA_CHECK(cudaMemcpyToSymbol(d_pf_eps, eps, sizeof(eps)));
        CUDA_CHECK(cudaMemcpyToSymbol(d_pf_sig6, sig6, sizeof(sig6))); CUDA_CHECK(cudaMemcpyToSymbol(d_pf_A, A, sizeof(A)));
        CUDA_CHECK(cudaMemcpyToSymbol(d_pf_invrc, invrc, sizeof(invrc))); CUDA_CHECK(cudaMemcpyToSymbol(d_pf_fpre, fpre, sizeof(fpre)));
        CUDA_CHECK(cudaMemcpyToSymbol(d_pair_fp64, &fp64, sizeof(int)));
        force_lanes_ok = !fp64;
        for (int t = 0; t < 5; t++) if (params->bond_is_fene[t]) force_lanes_ok = false;
    }
    {
        double nl2[9][9];
        for (int i = 0; i < 9; i++)
            for (int j = 0; j < 9; j++) {
                double rc = (params->pair_mode[i][j] != 0 && params->cutoff_sq[i][j] > 0.0) ? sqrt(params->cutoff_sq[i][j]) + params->skin : -1.0;
                nl2[i][j] = (rc > 0.0) ? rc * rc : -1.0;
            }
        CUDA_CHECK(cudaMemcpyToSymbol(d_nl_cut2, nl2, sizeof(nl2)));
        // fine-grid cell = largest DNA-DNA list cutoff; ribosome stencil from the largest ribosome-DNA list cutoff
        double fc2 = -1.0, rd2 = -1.0;
        for (int i = 3; i < 9; i++) {
            for (int j = 3; j < 9; j++) fc2 = std::max(fc2, nl2[i][j]);
            rd2 = std::max(rd2, nl2[2][i]);
        }
        w_fine_cut = (fc2 > 0.0) ? sqrt(fc2) : -1.0;
        w_ribo_dna_cut = (rd2 > 0.0) ? sqrt(rd2) : -1.0;
    }
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
    static int *w_device_loops_h_words = nullptr, *w_device_loops_d_words = nullptr;   // stop, iterations, evaluations in mapped memory
    if (w_device_loops_cg_on()) {
        if (!w_device_loops_h_words) {
            CUDA_CHECK(cudaHostAlloc((void **)&w_device_loops_h_words, 3 * sizeof(int), cudaHostAllocMapped));
            CUDA_CHECK(cudaHostGetDevicePointer((void **)&w_device_loops_d_words, w_device_loops_h_words, 0));
        }
        kd.d_stop = w_device_loops_d_words; kd.d_niter = w_device_loops_d_words + 1; kd.d_neval = w_device_loops_d_words + 2;
    }
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

    std::vector<cudaEvent_t> nl_events;
    while (h_stop == STOP_NEEDS_REBUILD) {
        cudaEvent_t nl_start, nl_stop;
        CUDA_CHECK(cudaEventCreate(&nl_start));
        CUDA_CHECK(cudaEventCreate(&nl_stop));
        CUDA_CHECK(cudaEventRecord(nl_start));

        // cg-two-grid check (env WCM_NLIST_CHECK=1, first build of the process only): build with the one-grid path into scratch
        // lists, then the two-grid path, and compare every atom's neighbour set on the host.
        static bool w_checked = false;
        int *chk_nn = nullptr, *chk_nl = nullptr;
        bool do_check = !w_checked && w_fine_cut > 0.0 && std::getenv("WCM_NLIST_CHECK");
        if (do_check) {
            w_checked = true;
            CUDA_CHECK(cudaMalloc(&chk_nn, (size_t)N * sizeof(int)));
            CUDA_CHECK(cudaMalloc(&chk_nl, (size_t)kd.max_neighs * cache->N_cap * sizeof(int)));
            double keep = w_fine_cut; w_fine_cut = -1.0;
            gpu_rebuild_neighbor_list(kd.x, kd.type, N, cache->d_excl_off, cache->d_excl_list, rcut2, box_lo, box_hi,
                ncx, ncy, ncz, cx, cy, cz, ncells, cache->d_cell_of_atom, cache->d_sorted_atoms,
                cache->d_cell_start, cache->d_cell_count, chk_nn, chk_nl, kd.max_neighs, cache->N_cap);
            w_fine_cut = keep;
        }
        gpu_rebuild_neighbor_list(
            kd.x, kd.type, N,
            cache->d_excl_off, cache->d_excl_list,
            rcut2, box_lo, box_hi,
            ncx, ncy, ncz, cx, cy, cz, ncells,
            cache->d_cell_of_atom, cache->d_sorted_atoms,
            cache->d_cell_start, cache->d_cell_count,
            kd.numneigh, kd.neighlist, kd.max_neighs,
            cache->N_cap);
        if (do_check) {
            std::vector<int> an(N), bn(N);
            CUDA_CHECK(cudaMemcpy(an.data(), chk_nn, N * sizeof(int), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(bn.data(), kd.numneigh, N * sizeof(int), cudaMemcpyDeviceToHost));
            int rows = 0; for (int i = 0; i < N; i++) rows = std::max(rows, std::max(an[i], bn[i]));
            std::vector<int> al((size_t)rows * cache->N_cap), bl((size_t)rows * cache->N_cap);
            CUDA_CHECK(cudaMemcpy(al.data(), chk_nl, al.size() * sizeof(int), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(bl.data(), kd.neighlist, bl.size() * sizeof(int), cudaMemcpyDeviceToHost));
            long long pairs = 0; int bad = 0;
            for (int i = 0; i < N; i++) {
                std::vector<int> x, y;
                for (int k = 0; k < an[i]; k++) x.push_back(al[(size_t)k * cache->N_cap + i]);
                for (int k = 0; k < bn[i]; k++) y.push_back(bl[(size_t)k * cache->N_cap + i]);
                std::sort(x.begin(), x.end()); std::sort(y.begin(), y.end());
                pairs += x.size(); if (x != y) bad++;
            }
            fprintf(stderr, "[WCM_NLIST_CHECK] atoms %d, one-grid pairs %lld, atoms whose neighbour set differs: %d\n", N, pairs, bad);
            cudaFree(chk_nn); cudaFree(chk_nl);
        }
        kd.order = wcm_balance_order(kd.type, kd.numneigh, N, max_blocks * block_size);

        CUDA_CHECK(cudaEventRecord(nl_stop));
        if (w_device_loops_cg_on()) {
            nl_events.push_back(nl_start); nl_events.push_back(nl_stop);   // read after the last segment
        } else {
        CUDA_CHECK(cudaEventSynchronize(nl_stop));
        float nl_ms = 0;
        CUDA_CHECK(cudaEventElapsedTime(&nl_ms, nl_start, nl_stop));
        nlist_time_ms += nl_ms;
        CUDA_CHECK(cudaEventDestroy(nl_start));
        CUDA_CHECK(cudaEventDestroy(nl_stop));
        }

        kd.start_iter  = total_iter;
        kd.start_neval = total_neval;
        kd.resume      = (launch_count > 0) ? 1 : 0;

        static const bool w_persistent = std::getenv("WCM_CG_PERSISTENT") != nullptr;
        if (w_persistent) {
            dim3 grid_dim(max_blocks), block_dim(block_size);
            void *args[] = { &kd };
            CUDA_CHECK(cudaLaunchCooperativeKernel(
                (void*)fused_cg_kernel, grid_dim, block_dim,
                args, block_size * sizeof(double), 0));
            CUDA_CHECK(cudaDeviceSynchronize());
        } else {
            w_split_cg_segment(kd);
        }

        if (launch_count == 0)
            CUDA_CHECK(cudaMemcpy(&result.energy_initial, kd.d_energy_init,
                                  sizeof(double), cudaMemcpyDeviceToHost));

        if (w_device_loops_cg_on() && !w_persistent) {   // written before the done flag the segment loop has seen
            h_stop = ((volatile int *)w_device_loops_h_words)[0]; total_iter = ((volatile int *)w_device_loops_h_words)[1]; total_neval = ((volatile int *)w_device_loops_h_words)[2];
        } else {
        CUDA_CHECK(cudaMemcpy(&h_stop,       kd.d_stop,         sizeof(int),    cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(&total_iter,   kd.d_niter,        sizeof(int),    cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(&total_neval,  kd.d_neval,        sizeof(int),    cudaMemcpyDeviceToHost));
        }
        launch_count++;

        if (h_stop == STOP_NEEDS_REBUILD)
            fprintf(stderr, "  Neighbor rebuild at iter %d (neval=%d)\n", total_iter, total_neval);
    }

    CUDA_CHECK(cudaEventRecord(t_stop));
    CUDA_CHECK(cudaEventSynchronize(t_stop));
    float ms = 0;
    CUDA_CHECK(cudaEventElapsedTime(&ms, t_start, t_stop));
    for (size_t e = 0; e + 1 < nl_events.size(); e += 2) {
        float nl_ms = 0;
        CUDA_CHECK(cudaEventElapsedTime(&nl_ms, nl_events[e], nl_events[e + 1]));
        nlist_time_ms += nl_ms;
        cudaEventDestroy(nl_events[e]); cudaEventDestroy(nl_events[e + 1]);
    }

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

// ============================================================================
//  CUDA context pre-warm
// ============================================================================
static std::thread w_prewarm;
void fused_prewarm_cuda_context_async()
{
    if (w_prewarm.joinable()) return;
    w_prewarm = std::thread([] { if (cudaSetDevice(0) == cudaSuccess) cudaFree(0); });
}
void fused_prewarm_join()
{
    if (w_prewarm.joinable()) w_prewarm.join();
}

// ============================================================================
//  boundary atoms in space-filling-curve order
// ============================================================================
// Boundary atoms (type 1, no bonds or angles) come in icosphere generation order, which is spatially scattered, so a DNA atom near
// the membrane gathered its ~100 boundary partners from as many cache lines. Their slots are reordered among themselves along a
// Morton curve before the minimisation and restored after it. They carry no topology or exclusions, keep type 1, and have zero
// force, so they do not move (x0 + alpha * 0): restoring the saved positions is exact. Other atoms keep their indices.
static inline unsigned long long w_morton_spread(unsigned int v) {
    unsigned long long x = v & 0x1fffff;
    x = (x | x << 32) & 0x1f00000000ffffULL; x = (x | x << 16) & 0x1f0000ff0000ffULL;
    x = (x | x << 8) & 0x100f00f00f00f00fULL; x = (x | x << 4) & 0x10c30c30c30c30c3ULL; x = (x | x << 2) & 0x1249249249249249ULL;
    return x;
}

FusedMinResult fused_cg_minimize(
    FusedMinSystem *sys,
    const FusedMinParams *params,
    const double box_lo[3], const double box_hi[3],
    int gpu_id,
    FusedDeviceCache *cache)
{
    const int N = sys->N;
    std::vector<char> topo(N, 0);
    for (int b = 0; b < sys->nbonds; b++) { topo[sys->bond_i[b]] = 1; topo[sys->bond_j[b]] = 1; }
    for (int a = 0; a < sys->nangles; a++) { topo[sys->angle_i[a]] = 1; topo[sys->angle_j[a]] = 1; topo[sys->angle_k[a]] = 1; }
    std::vector<int> slots;
    for (int i = 0; i < N; i++) if (sys->type[i] == 1 && !topo[i]) slots.push_back(i);
    if (slots.size() < 2) return fused_cg_minimize_impl(sys, params, box_lo, box_hi, gpu_id, cache);
    std::vector<unsigned long long> code(slots.size());
    for (size_t k = 0; k < slots.size(); k++) {
        unsigned int q[3];
        for (int d = 0; d < 3; d++) {
            double L = box_hi[d] - box_lo[d], t = (L > 0.0) ? (sys->x[3*slots[k] + d] - box_lo[d]) / L : 0.0;
            t = std::min(std::max(t, 0.0), 1.0);
            q[d] = (unsigned int)(t * 2097151.0);
        }
        code[k] = w_morton_spread(q[0]) | (w_morton_spread(q[1]) << 1) | (w_morton_spread(q[2]) << 2);
    }
    std::vector<int> order(slots.size());
    for (size_t k = 0; k < order.size(); k++) order[k] = (int)k;
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return code[a] < code[b]; });
    // slot slots[k] receives the boundary atom originally in slots[order[k]]
    std::vector<double> x_orig(3 * slots.size()), xref_orig(3 * slots.size());
    for (size_t k = 0; k < slots.size(); k++)
        for (int d = 0; d < 3; d++) { x_orig[3*k+d] = sys->x[3*slots[k]+d]; xref_orig[3*k+d] = sys->xref[3*slots[k]+d]; }
    for (size_t k = 0; k < slots.size(); k++)
        for (int d = 0; d < 3; d++) { sys->x[3*slots[k]+d] = x_orig[3*order[k]+d]; sys->xref[3*slots[k]+d] = xref_orig[3*order[k]+d]; }
    FusedMinResult r = fused_cg_minimize_impl(sys, params, box_lo, box_hi, gpu_id, cache);
    for (size_t k = 0; k < slots.size(); k++)
        for (int d = 0; d < 3; d++) { sys->x[3*slots[k]+d] = x_orig[3*k+d]; sys->xref[3*slots[k]+d] = xref_orig[3*k+d]; }
    return r;
}
