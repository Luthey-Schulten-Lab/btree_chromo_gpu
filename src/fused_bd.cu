// Author(s): Ron Acda (using an iterative LLM-guided workflow, https://github.com/quarkron/iterative-hillclimber/tree/main)
// fused Brownian dynamics for the soft-harmonic DNA run (see include/fused_bd.h).
//
// Why: the LAMMPS Kokkos BD run was ~190 us per step at t=2996 (149,638 atoms, 23,333 steps per hook), of which ~85 us were GPU
// kernels (soft pair, bond, angle, brownian, neighbour check, ...) and the rest per-kernel launch/fence/host bookkeeping. Here a
// step is ONE kernel (forces + integration for one atom per thread, positions double-buffered), launched back to back from the
// host; the host synchronises only on check steps to read the largest displacement since the last list build.
#include "wcm_arch.cuh"
#include "fused_bd.h"

#include <cuda_runtime.h>
#include <cub/device/device_scan.cuh>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cstring>

#define BD_CHECK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
    fprintf(stderr, "fused_bd CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e_)); return -1; } } while (0)

static constexpr double BD_PI = 3.14159265358979323846;

__constant__ double c_bd_A[FUSED_BD_MAXT][FUSED_BD_MAXT];
__constant__ double c_bd_rc[FUSED_BD_MAXT][FUSED_BD_MAXT];
__constant__ double c_bd_rc2[FUSED_BD_MAXT][FUSED_BD_MAXT];     // interaction cutoff^2, <= 0: none
__constant__ double c_bd_nl2[FUSED_BD_MAXT][FUSED_BD_MAXT];     // (rc + skin)^2, <= 0: never listed
__constant__ double c_bd_bK[FUSED_BD_MAXB], c_bd_br0[FUSED_BD_MAXB], c_bd_aK[FUSED_BD_MAXB];
// the mixed-fp32 change (mixed precision): FP32 pair term (cutoff^2, 1/rc, A pi / rc) on the FP64 coordinate difference, FP64 accumulation;
// WCM_FP64=1 keeps the FP64 term.
__constant__ float c_bdf_rc2[FUSED_BD_MAXT][FUSED_BD_MAXT], c_bdf_invrc[FUSED_BD_MAXT][FUSED_BD_MAXT], c_bdf_fpre[FUSED_BD_MAXT][FUSED_BD_MAXT];
__constant__ int c_bd_fp64;

// ---------------------------------------------------------------------------------------------------------------------------
// Philox4x32-10 (Random123), stateless: counter (step, atom, 0, 0), key (seed)
// ---------------------------------------------------------------------------------------------------------------------------
__device__ __forceinline__ uint4 bd_philox(uint4 c, uint2 k) {
    const uint32_t M0 = 0xD2511F53u, M1 = 0xCD9E8D57u, W0 = 0x9E3779B9u, W1 = 0xBB67AE85u;
    #pragma unroll
    for (int r = 0; r < 10; r++) {
        uint32_t hi0 = __umulhi(M0, c.x), lo0 = M0 * c.x;
        uint32_t hi1 = __umulhi(M1, c.z), lo1 = M1 * c.z;
        c = make_uint4(hi1 ^ c.y ^ k.x, lo1, hi0 ^ c.w ^ k.y, lo0);
        k.x += W0; k.y += W1;
    }
    return c;
}
__device__ __forceinline__ double bd_u01(uint32_t r) { return ((double)r + 0.5) * 2.3283064365386963e-10; }   // (0,1)

// ---------------------------------------------------------------------------------------------------------------------------
// neighbour list: one uniform cell grid (cell >= largest list cutoff among the types present), lists only for mobile atoms
// ---------------------------------------------------------------------------------------------------------------------------
// positions are double4 (x, y, z, atom type) so a neighbour is one aligned 32-byte load; same arithmetic.
__global__ void bd_cell_kernel(const double4 *x, int N, int ncx, int ncy, int ncz, double cx, double cy, double cz,
                               double lx, double ly, double lz, int *cell_of, int *count) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    const double4 q = x[i];
    int ix = min(max((int)((q.x - lx) / cx), 0), ncx - 1);
    int iy = min(max((int)((q.y - ly) / cy), 0), ncy - 1);
    int iz = min(max((int)((q.z - lz) / cz), 0), ncz - 1);
    int c = (ix * ncy + iy) * ncz + iz;
    cell_of[i] = c;
    atomicAdd(&count[c], 1);
}
__global__ void bd_scatter_kernel(const int *cell_of, int N, const int *start, int *fill, int *sorted) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    int c = cell_of[i];
    sorted[start[c] + atomicAdd(&fill[c], 1)] = i;
}
__device__ __forceinline__ bool bd_excluded(const int *list, int s, int e, int j) {
    while (s < e) { int m = (s + e) >> 1; int v = list[m]; if (v == j) return true; if (v < j) s = m + 1; else e = m; }
    return false;
}
__global__ void bd_build_kernel(const double4 *x, const int *type, const unsigned char *mobile, int N,
                                const int *sorted, const int *start, int ncx, int ncy, int ncz, double cx, double cy, double cz,
                                double lx, double ly, double lz, const int *eoff, const int *elist,
                                int *numneigh, int *nl, int max_neighs, int stride, int *overflow) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    if (!mobile[i]) { numneigh[i] = 0; return; }
    int ti = type[i];
    const double4 qi = x[i];
    double xi = qi.x, yi = qi.y, zi = qi.z;
    int ix = min(max((int)((xi - lx) / cx), 0), ncx - 1);
    int iy = min(max((int)((yi - ly) / cy), 0), ncy - 1);
    int iz = min(max((int)((zi - lz) / cz), 0), ncz - 1);
    int es = eoff[i], ee = eoff[i + 1], nn = 0;
    for (int dx = -1; dx <= 1; dx++) for (int dy = -1; dy <= 1; dy++) for (int dz = -1; dz <= 1; dz++) {
        int jx = ix + dx, jy = iy + dy, jz = iz + dz;
        if (jx < 0 || jx >= ncx || jy < 0 || jy >= ncy || jz < 0 || jz >= ncz) continue;
        int c = (jx * ncy + jy) * ncz + jz;
        for (int s = start[c], e = start[c + 1]; s < e; s++) {
            int j = sorted[s];
            if (j == i) continue;
            const double4 qj = x[j];
            double nl2 = c_bd_nl2[ti][(int)qj.w];
            if (nl2 <= 0.0) continue;
            double ddx = xi - qj.x, ddy = yi - qj.y, ddz = zi - qj.z;
            if (ddx*ddx + ddy*ddy + ddz*ddz < nl2 && !bd_excluded(elist, es, ee, j)) {
                if (nn < max_neighs) nl[(size_t)nn * stride + i] = j;
                nn++;
            }
        }
    }
    numneigh[i] = min(nn, max_neighs);
    if (nn > max_neighs) atomicMax(overflow, nn);
}

#include "wcm_bd_nlist_warp.cuh"

// ---------------------------------------------------------------------------------------------------------------------------
// one BD step: forces on atom a from its list, bonds and angles, then x_out = x_in + dt (g1 f + g2 (u - 0.5))
// ---------------------------------------------------------------------------------------------------------------------------
struct BDStepArgs {
    const double4 *xin; double4 *xout; const double4 *xbuild;
    const int *type; const unsigned char *mobile; int N;
    const int *numneigh, *nl; int stride;
    const int *boff, *bdat;       // per-atom bonds: (type, other)
    const int *aoff, *adat;       // per-atom angles: (type, i, j, k)
    double dt, g1, g2;
    uint2 key; uint32_t step_lo, step_hi;
    int do_check; unsigned long long *maxd2_bits;
    const float4 *xinf; float4 *xoutf;   // FP32 copy of the positions (x, y, z, type) for every gathered load
    const int *halt;                     // set on the device when a list rebuild is due; steps then do nothing
};

__device__ __forceinline__ void bd_pair(int ta, int tj, double dx, double dy, double dz, double &fx, double &fy, double &fz) {
    if (!c_bd_fp64) {
        const float x = (float)dx, y = (float)dy, z = (float)dz, r2 = x*x + y*y + z*z;
        if (r2 >= c_bdf_rc2[ta][tj] || r2 <= 0.0f) return;
        const float r = sqrtf(r2), fpair = c_bdf_fpre[ta][tj] * sinpif(r * c_bdf_invrc[ta][tj]) / r;
        fx += (double)(fpair * x); fy += (double)(fpair * y); fz += (double)(fpair * z);
        return;
    }
    double r2 = dx*dx + dy*dy + dz*dz;
    if (r2 >= c_bd_rc2[ta][tj] || r2 <= 0.0) return;
    double r = sqrt(r2), rc = c_bd_rc[ta][tj];
    double fpair = c_bd_A[ta][tj] * sin(BD_PI * r / rc) * BD_PI / rc / r;       // LAMMPS pair soft
    fx += fpair * dx; fy += fpair * dy; fz += fpair * dz;
}

// 3 blocks/SM (80 registers, 64 B spill) instead of 2 (106 registers)
__global__ void __launch_bounds__(256, 3) bd_step_kernel(BDStepArgs A) {
    int a = blockIdx.x * blockDim.x + threadIdx.x;
    double d2 = 0.0;
    if (a < A.N) {
        const double4 qa = A.xin[a];
        double xa = qa.x, ya = qa.y, za = qa.z;
        if (!A.mobile[a]) {
            A.xout[a] = qa;
        } else {
            int ta = (int)qa.w;
            double fx = 0.0, fy = 0.0, fz = 0.0;
            int nn = A.numneigh[a];
            const int *nl = A.nl + a;
            for (int jj = 0; jj < nn; jj += 4) {
                int m = nn - jj;
                int j0 = nl[(size_t)jj * A.stride];
                int j1 = (m > 1) ? nl[(size_t)(jj + 1) * A.stride] : j0;
                int j2 = (m > 2) ? nl[(size_t)(jj + 2) * A.stride] : j0;
                int j3 = (m > 3) ? nl[(size_t)(jj + 3) * A.stride] : j0;
                const double4 q0 = A.xin[j0], q1 = A.xin[j1], q2 = A.xin[j2], q3 = A.xin[j3];
                int t0 = (int)q0.w, t1 = (int)q1.w, t2 = (int)q2.w, t3 = (int)q3.w;
                double x0 = q0.x, y0 = q0.y, z0 = q0.z;
                double x1 = q1.x, y1 = q1.y, z1 = q1.z;
                double x2 = q2.x, y2 = q2.y, z2 = q2.z;
                double x3 = q3.x, y3 = q3.y, z3 = q3.z;
                bd_pair(ta, t0, xa - x0, ya - y0, za - z0, fx, fy, fz);
                if (m > 1) bd_pair(ta, t1, xa - x1, ya - y1, za - z1, fx, fy, fz);
                if (m > 2) bd_pair(ta, t2, xa - x2, ya - y2, za - z2, fx, fy, fz);
                if (m > 3) bd_pair(ta, t3, xa - x3, ya - y3, za - z3, fx, fy, fz);
            }
            // bonds: LAMMPS bond harmonic, E = K (r - r0)^2, fbond = -2 K (r - r0) / r
            if (!c_bd_fp64) {   // bonds and angles in FP32 on FP64 differences, FP64 accumulation
                for (int b = A.boff[a], be = A.boff[a + 1]; b < be; b++) {
                    int bt = A.bdat[2*b], o = A.bdat[2*b + 1];
                    const double4 qo = A.xin[o];
                    const float x = (float)(xa - qo.x), y = (float)(ya - qo.y), z = (float)(za - qo.z);
                    const float r = sqrtf(x*x + y*y + z*z);
                    if (r > 0.0f) {
                        const float fb = -2.0f * (float)c_bd_bK[bt] * (r - (float)c_bd_br0[bt]) / r;
                        fx += (double)(fb * x); fy += (double)(fb * y); fz += (double)(fb * z);
                    }
                }
                for (int q = A.aoff[a], qe = A.aoff[a + 1]; q < qe; q++) {
                    int at = A.adat[4*q], i1 = A.adat[4*q+1], i2 = A.adat[4*q+2], i3 = A.adat[4*q+3];
                    const double4 p1 = A.xin[i1], p2 = A.xin[i2], p3 = A.xin[i3];
                    const float d1x = (float)(p1.x - p2.x), d1y = (float)(p1.y - p2.y), d1z = (float)(p1.z - p2.z);
                    const float d2x = (float)(p3.x - p2.x), d2y = (float)(p3.y - p2.y), d2z = (float)(p3.z - p2.z);
                    const float rsq1 = d1x*d1x + d1y*d1y + d1z*d1z, rsq2 = d2x*d2x + d2y*d2y + d2z*d2z;
                    const float r1 = sqrtf(rsq1), r2 = sqrtf(rsq2);
                    if (r1 <= 0.0f || r2 <= 0.0f) continue;
                    float c = (d1x*d2x + d1y*d2y + d1z*d2z) / (r1 * r2);
                    c = fmaxf(-1.0f, fminf(1.0f, c));
                    const float K = (float)c_bd_aK[at];
                    const float a11 = K * c / rsq1, a12 = -K / (r1 * r2), a22 = K * c / rsq2;
                    const float f1x = a11 * d1x + a12 * d2x, f1y = a11 * d1y + a12 * d2y, f1z = a11 * d1z + a12 * d2z;
                    const float f3x = a22 * d2x + a12 * d1x, f3y = a22 * d2y + a12 * d1y, f3z = a22 * d2z + a12 * d1z;
                    if (a == i1)      { fx += (double)f1x; fy += (double)f1y; fz += (double)f1z; }
                    else if (a == i3) { fx += (double)f3x; fy += (double)f3y; fz += (double)f3z; }
                    else              { fx -= (double)(f1x + f3x); fy -= (double)(f1y + f3y); fz -= (double)(f1z + f3z); }
                }
            } else {
            for (int b = A.boff[a], be = A.boff[a + 1]; b < be; b++) {
                int bt = A.bdat[2*b], o = A.bdat[2*b + 1];
                const double4 qo = A.xin[o];
                double dx = xa - qo.x, dy = ya - qo.y, dz = za - qo.z;
                double r = sqrt(dx*dx + dy*dy + dz*dz);
                if (r > 0.0) {
                    double fb = -2.0 * c_bd_bK[bt] * (r - c_bd_br0[bt]) / r;
                    fx += fb * dx; fy += fb * dy; fz += fb * dz;
                }
            }
            // angles: LAMMPS angle cosine, E = K (1 + c); f1 = a11 d1 + a12 d2, f3 = a22 d2 + a12 d1, f2 = -(f1 + f3)
            for (int q = A.aoff[a], qe = A.aoff[a + 1]; q < qe; q++) {
                int at = A.adat[4*q], i1 = A.adat[4*q+1], i2 = A.adat[4*q+2], i3 = A.adat[4*q+3];
                const double4 p1 = A.xin[i1], p2 = A.xin[i2], p3 = A.xin[i3];
                double d1x = p1.x - p2.x, d1y = p1.y - p2.y, d1z = p1.z - p2.z;
                double d2x = p3.x - p2.x, d2y = p3.y - p2.y, d2z = p3.z - p2.z;
                double rsq1 = d1x*d1x + d1y*d1y + d1z*d1z, rsq2 = d2x*d2x + d2y*d2y + d2z*d2z;
                double r1 = sqrt(rsq1), r2 = sqrt(rsq2);
                if (r1 <= 0.0 || r2 <= 0.0) continue;
                double c = (d1x*d2x + d1y*d2y + d1z*d2z) / (r1 * r2);
                c = fmax(-1.0, fmin(1.0, c));
                double K = c_bd_aK[at];
                double a11 = K * c / rsq1, a12 = -K / (r1 * r2), a22 = K * c / rsq2;
                double f1x = a11 * d1x + a12 * d2x, f1y = a11 * d1y + a12 * d2y, f1z = a11 * d1z + a12 * d2z;
                double f3x = a22 * d2x + a12 * d1x, f3y = a22 * d2y + a12 * d1y, f3z = a22 * d2z + a12 * d1z;
                if (a == i1)      { fx += f1x; fy += f1y; fz += f1z; }
                else if (a == i3) { fx += f3x; fy += f3y; fz += f3z; }
                else              { fx -= f1x + f3x; fy -= f1y + f3y; fz -= f1z + f3z; }
            }
            }   // the fp32-bonded change FP64 path
            // fix brownian (uniform noise)
            uint4 r = bd_philox(make_uint4(A.step_lo, A.step_hi, (uint32_t)a, 0u), A.key);
            double nx = qa.x + A.dt * (A.g1 * fx + A.g2 * (bd_u01(r.x) - 0.5));
            double ny = qa.y + A.dt * (A.g1 * fy + A.g2 * (bd_u01(r.y) - 0.5));
            double nz = qa.z + A.dt * (A.g1 * fz + A.g2 * (bd_u01(r.z) - 0.5));
            A.xout[a] = make_double4(nx, ny, nz, qa.w);
            if (A.do_check) {
                const double4 qb = A.xbuild[a];
                double ex = nx - qb.x, ey = ny - qb.y, ez = nz - qb.z;
                d2 = ex*ex + ey*ey + ez*ez;
            }
        }
    }
    if (A.do_check) {   // warp max, one atomic per warp (non-negative doubles order like their bit patterns)
        for (int o = 16; o > 0; o >>= 1) d2 = fmax(d2, __shfl_down_sync(0xffffffffu, d2, o));
        if ((threadIdx.x & 31) == 0 && d2 > 0.0) atomicMax(A.maxd2_bits, (unsigned long long)__double_as_longlong(d2));
    }
}

// the bd-fp32-gather change (mixed precision): every gathered position (neighbours, bond partners, angle atoms, the atom itself for its
// differences) comes from an FP32 copy (float4 = 16 B instead of 32 B), U neighbours in flight per round; forces as 117/118
// (FP32 terms, FP64 accumulation); the FP64 positions are still integrated and written, and the FP32 copy is rewritten from them
// every step, so no rounding accumulates in the trajectory. Coordinates |x| < 2200 A -> copy resolution <= 2.4e-4 A.
template <int U, bool HALT = false>
__global__ void __launch_bounds__(256, 3) bd_step_kernel_f(BDStepArgs A) {
    wcm_grid_dependency_sync();
    if (HALT && *A.halt) return;       // a rebuild is pending; the host re-enqueues this step after it
    int a = blockIdx.x * blockDim.x + threadIdx.x;
    double d2 = 0.0;
    if (a < A.N) {
        const double4 qa = A.xin[a];
        if (!A.mobile[a]) {
            A.xout[a] = qa; A.xoutf[a] = A.xinf[a];
        } else {
            const float4 fa = A.xinf[a];
            int ta = (int)qa.w;
            double fx = 0.0, fy = 0.0, fz = 0.0;
            int nn = A.numneigh[a];
            const int *nl = A.nl + a;
            for (int jj = 0; jj < nn; jj += U) {
                int m = nn - jj;
                int j[U]; float4 q[U];
                #pragma unroll
                for (int u = 0; u < U; u++) j[u] = (u == 0 || m > u) ? nl[(size_t)(jj + u) * A.stride] : -1;
                #pragma unroll
                for (int u = 0; u < U; u++) q[u] = A.xinf[j[u] >= 0 ? j[u] : j[0]];
                #pragma unroll
                for (int u = 0; u < U; u++) {
                    if (u > 0 && m <= u) break;
                    const int tj = (int)q[u].w;
                    const float x = fa.x - q[u].x, y = fa.y - q[u].y, z = fa.z - q[u].z, r2 = x*x + y*y + z*z;
                    if (r2 >= c_bdf_rc2[ta][tj] || r2 <= 0.0f) continue;
                    const float r = sqrtf(r2), fpair = c_bdf_fpre[ta][tj] * sinpif(r * c_bdf_invrc[ta][tj]) / r;
                    fx += (double)(fpair * x); fy += (double)(fpair * y); fz += (double)(fpair * z);
                }
            }
            for (int b = A.boff[a], be = A.boff[a + 1]; b < be; b++) {
                int bt = A.bdat[2*b], o = A.bdat[2*b + 1];
                const float4 qo = A.xinf[o];
                const float x = fa.x - qo.x, y = fa.y - qo.y, z = fa.z - qo.z;
                const float r = sqrtf(x*x + y*y + z*z);
                if (r > 0.0f) {
                    const float fb = -2.0f * (float)c_bd_bK[bt] * (r - (float)c_bd_br0[bt]) / r;
                    fx += (double)(fb * x); fy += (double)(fb * y); fz += (double)(fb * z);
                }
            }
            for (int q = A.aoff[a], qe = A.aoff[a + 1]; q < qe; q++) {
                int at = A.adat[4*q], i1 = A.adat[4*q+1], i2 = A.adat[4*q+2], i3 = A.adat[4*q+3];
                const float4 p1 = A.xinf[i1], p2 = A.xinf[i2], p3 = A.xinf[i3];
                const float d1x = p1.x - p2.x, d1y = p1.y - p2.y, d1z = p1.z - p2.z;
                const float d2x = p3.x - p2.x, d2y = p3.y - p2.y, d2z = p3.z - p2.z;
                const float rsq1 = d1x*d1x + d1y*d1y + d1z*d1z, rsq2 = d2x*d2x + d2y*d2y + d2z*d2z;
                const float r1 = sqrtf(rsq1), r2 = sqrtf(rsq2);
                if (r1 <= 0.0f || r2 <= 0.0f) continue;
                float c = (d1x*d2x + d1y*d2y + d1z*d2z) / (r1 * r2);
                c = fmaxf(-1.0f, fminf(1.0f, c));
                const float K = (float)c_bd_aK[at];
                const float a11 = K * c / rsq1, a12 = -K / (r1 * r2), a22 = K * c / rsq2;
                const float f1x = a11 * d1x + a12 * d2x, f1y = a11 * d1y + a12 * d2y, f1z = a11 * d1z + a12 * d2z;
                const float f3x = a22 * d2x + a12 * d1x, f3y = a22 * d2y + a12 * d1y, f3z = a22 * d2z + a12 * d1z;
                if (a == i1)      { fx += (double)f1x; fy += (double)f1y; fz += (double)f1z; }
                else if (a == i3) { fx += (double)f3x; fy += (double)f3y; fz += (double)f3z; }
                else              { fx -= (double)(f1x + f3x); fy -= (double)(f1y + f3y); fz -= (double)(f1z + f3z); }
            }
            uint4 r = bd_philox(make_uint4(A.step_lo, A.step_hi, (uint32_t)a, 0u), A.key);
            double nx = qa.x + A.dt * (A.g1 * fx + A.g2 * (bd_u01(r.x) - 0.5));
            double ny = qa.y + A.dt * (A.g1 * fy + A.g2 * (bd_u01(r.y) - 0.5));
            double nz = qa.z + A.dt * (A.g1 * fz + A.g2 * (bd_u01(r.z) - 0.5));
            A.xout[a] = make_double4(nx, ny, nz, qa.w);
            A.xoutf[a] = make_float4((float)nx, (float)ny, (float)nz, fa.w);
            if (A.do_check) {
                const double4 qb = A.xbuild[a];
                double ex = nx - qb.x, ey = ny - qb.y, ez = nz - qb.z;
                d2 = ex*ex + ey*ey + ez*ez;
            }
        }
    }
    if (A.do_check) {
        for (int o = 16; o > 0; o >>= 1) d2 = fmax(d2, __shfl_down_sync(0xffffffffu, d2, o));
        if ((threadIdx.x & 31) == 0 && d2 > 0.0) atomicMax(A.maxd2_bits, (unsigned long long)__double_as_longlong(d2));
    }
}

// the displacement test of a check step on the device (the host used to copy maxd2 back after every check, a full
// round trip every 5 steps). If no rebuild is pending: when the step's maximum squared displacement exceeds the trigger, set the
// halt flag (the following steps become no-ops) and publish the step to resume at and maxd2 to mapped host memory; always reset
// maxd2 for the next check window, as the host's memset did.
__global__ void bd_check_kernel(unsigned long long *maxd2_bits, double thr2, int *halt, long long resume_step, volatile long long *h_info) {
    wcm_grid_dependency_sync();
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    if (*halt) return;
    const unsigned long long bits = *maxd2_bits;
    double md2; memcpy(&md2, &bits, sizeof(md2));
    if (md2 > thr2) {
        *halt = 1;
        h_info[1] = (long long)bits;
        __threadfence_system();
        h_info[0] = resume_step;
    }
    *maxd2_bits = 0ull;
}

// list overflow of a build decided on the device: halt the following steps and publish the overflow count (the host
// used to copy it back after every build); the host then grows the list and builds again before resuming, as the build loop did.
__global__ void bd_ovf_kernel(const int *ovf, int *halt, long long resume_step, volatile long long *h_info) {
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    const int o = *ovf;
    if (o > 0) {
        *halt = 1;
        h_info[2] = o;
        __threadfence_system();
        h_info[0] = resume_step;
    }
}

// ---------------------------------------------------------------------------------------------------------------------------
// host
// ---------------------------------------------------------------------------------------------------------------------------
// programmatic dependent launch (PDL): kernels of the CG tick chain / BD step chain are launched with programmatic stream
// serialization, so the next grid is launched while the previous one drains; each such kernel starts with
// cudaGridDependencySynchronize() (waits for the previous grid's completion and memory flush; a no-op without the attribute).
// WCM_PDL_OFF=1, or a GPU below compute capability 9.0 = ordinary launches (include/wcm_arch.cuh).
static bool bd_pdl_on() { return wcm_pdl_enabled(); }   // include/wcm_arch.cuh: compute capability >= 9.0, WCM_PDL_OFF unset
template <typename... KArgs, typename... Args>
static void bd_launch(void (*k)(KArgs...), dim3 g, dim3 b, Args... args) {
    cudaLaunchConfig_t cfg = {};
    cfg.gridDim = g; cfg.blockDim = b; cfg.dynamicSmemBytes = 0; cfg.stream = 0;
    cudaLaunchAttribute at[1];
    at[0].id = cudaLaunchAttributeProgrammaticStreamSerialization;
    at[0].val.programmaticStreamSerializationAllowed = 1;
    cfg.attrs = at; cfg.numAttrs = bd_pdl_on() ? 1 : 0;
    cudaLaunchKernelEx(&cfg, k, static_cast<KArgs>(args)...);
}

int fused_bd_run(FusedBDSystem *s, long nsteps, FusedBDStats *st) {
    using clk = std::chrono::steady_clock;
    auto t_all = clk::now();
    const int N = s->N;
    *st = FusedBDStats{0, 0.0, 0.0, 0.0, 0.0};
    if (N <= 0 || nsteps <= 0) return 0;
    BD_CHECK(cudaSetDevice(0));

    // parameters -> constant memory; list cutoff over the types present
    std::vector<char> present(FUSED_BD_MAXT, 0);
    for (int i = 0; i < N; i++) if (s->type[i] > 0 && s->type[i] < FUSED_BD_MAXT) present[s->type[i]] = 1;
    double rc2[FUSED_BD_MAXT][FUSED_BD_MAXT], nl2[FUSED_BD_MAXT][FUSED_BD_MAXT];
    double cutmax = 0.0;
    for (int i = 0; i < FUSED_BD_MAXT; i++) for (int j = 0; j < FUSED_BD_MAXT; j++) {
        bool on = s->pair_A[i][j] != 0.0 && s->pair_rc[i][j] > 0.0;
        rc2[i][j] = on ? s->pair_rc[i][j] * s->pair_rc[i][j] : -1.0;
        nl2[i][j] = on ? (s->pair_rc[i][j] + s->skin) * (s->pair_rc[i][j] + s->skin) : -1.0;
        if (on && present[i] && present[j]) cutmax = std::max(cutmax, s->pair_rc[i][j] + s->skin);
    }
    if (cutmax <= 0.0) cutmax = s->skin > 0.0 ? s->skin : 1.0;
    BD_CHECK(cudaMemcpyToSymbol(c_bd_A, s->pair_A, sizeof(rc2)));
    {
        static const int fp64 = std::getenv("WCM_FP64") != nullptr;
        float frc2[FUSED_BD_MAXT][FUSED_BD_MAXT], finv[FUSED_BD_MAXT][FUSED_BD_MAXT], fpre[FUSED_BD_MAXT][FUSED_BD_MAXT];
        for (int i = 0; i < FUSED_BD_MAXT; i++) for (int j = 0; j < FUSED_BD_MAXT; j++) {
            const bool on = rc2[i][j] > 0.0;
            frc2[i][j] = on ? (float)rc2[i][j] : -1.0f;
            finv[i][j] = on ? (float)(1.0 / s->pair_rc[i][j]) : 0.0f;
            fpre[i][j] = on ? (float)(s->pair_A[i][j] * BD_PI / s->pair_rc[i][j]) : 0.0f;
        }
        BD_CHECK(cudaMemcpyToSymbol(c_bdf_rc2, frc2, sizeof(frc2))); BD_CHECK(cudaMemcpyToSymbol(c_bdf_invrc, finv, sizeof(finv)));
        BD_CHECK(cudaMemcpyToSymbol(c_bdf_fpre, fpre, sizeof(fpre))); BD_CHECK(cudaMemcpyToSymbol(c_bd_fp64, &fp64, sizeof(int)));
    }
    BD_CHECK(cudaMemcpyToSymbol(c_bd_rc, s->pair_rc, sizeof(rc2)));
    BD_CHECK(cudaMemcpyToSymbol(c_bd_rc2, rc2, sizeof(rc2)));
    BD_CHECK(cudaMemcpyToSymbol(c_bd_nl2, nl2, sizeof(nl2)));
    BD_CHECK(cudaMemcpyToSymbol(c_bd_bK, s->bond_K, sizeof(s->bond_K)));
    BD_CHECK(cudaMemcpyToSymbol(c_bd_br0, s->bond_r0, sizeof(s->bond_r0)));
    BD_CHECK(cudaMemcpyToSymbol(c_bd_aK, s->angle_K, sizeof(s->angle_K)));

    // per-atom bond and angle tables (only mobile atoms need them; built for all, in list order)
    std::vector<int> boff(N + 1, 0), aoff(N + 1, 0);
    for (int b = 0; b < s->nbonds; b++) { boff[s->bond_i[b] + 1]++; boff[s->bond_j[b] + 1]++; }
    for (int q = 0; q < s->nangles; q++) { aoff[s->angle_i[q] + 1]++; aoff[s->angle_j[q] + 1]++; aoff[s->angle_k[q] + 1]++; }
    for (int i = 0; i < N; i++) { boff[i + 1] += boff[i]; aoff[i + 1] += aoff[i]; }
    std::vector<int> bdat(2 * (size_t)boff[N]), adat(4 * (size_t)aoff[N]), bpos(boff.begin(), boff.end() - 1), apos(aoff.begin(), aoff.end() - 1);
    for (int b = 0; b < s->nbonds; b++) {
        int i = s->bond_i[b], j = s->bond_j[b], t = s->bond_type[b];
        int p = bpos[i]++; bdat[2*p] = t; bdat[2*p+1] = j;
        p = bpos[j]++;     bdat[2*p] = t; bdat[2*p+1] = i;
    }
    for (int q = 0; q < s->nangles; q++) {
        int ids[3] = {s->angle_i[q], s->angle_j[q], s->angle_k[q]};
        for (int w = 0; w < 3; w++) {
            int p = apos[ids[w]]++;
            adat[4*p] = s->angle_type[q]; adat[4*p+1] = ids[0]; adat[4*p+2] = ids[1]; adat[4*p+3] = ids[2];
        }
    }

    // device buffers
    double4 *d_x0, *d_x1, *d_xb; int *d_type, *d_nn, *d_nl = nullptr, *d_boff, *d_bdat, *d_aoff, *d_adat, *d_eoff, *d_elist;
    unsigned char *d_mob; int *d_cellof, *d_sorted, *d_start, *d_count, *d_ovf; unsigned long long *d_maxd2;
    const size_t n4 = (size_t)N * sizeof(double4);
    BD_CHECK(cudaMalloc(&d_x0, n4)); BD_CHECK(cudaMalloc(&d_x1, n4)); BD_CHECK(cudaMalloc(&d_xb, n4));
    std::vector<double4> hx((size_t)N);
    for (int i = 0; i < N; i++) hx[i] = make_double4(s->x[3*i], s->x[3*i+1], s->x[3*i+2], (double)s->type[i]);
    BD_CHECK(cudaMalloc(&d_type, N * sizeof(int))); BD_CHECK(cudaMalloc(&d_mob, N)); BD_CHECK(cudaMalloc(&d_nn, N * sizeof(int)));
    BD_CHECK(cudaMalloc(&d_boff, (N + 1) * sizeof(int))); BD_CHECK(cudaMalloc(&d_bdat, std::max<size_t>(1, bdat.size()) * sizeof(int)));
    BD_CHECK(cudaMalloc(&d_aoff, (N + 1) * sizeof(int))); BD_CHECK(cudaMalloc(&d_adat, std::max<size_t>(1, adat.size()) * sizeof(int)));
    const int ne = s->excl_off[N];
    BD_CHECK(cudaMalloc(&d_eoff, (N + 1) * sizeof(int))); BD_CHECK(cudaMalloc(&d_elist, std::max(1, ne) * sizeof(int)));
    BD_CHECK(cudaMalloc(&d_cellof, N * sizeof(int))); BD_CHECK(cudaMalloc(&d_sorted, N * sizeof(int)));
    BD_CHECK(cudaMalloc(&d_ovf, sizeof(int))); BD_CHECK(cudaMalloc(&d_maxd2, sizeof(unsigned long long)));
    BD_CHECK(cudaMemcpy(d_x0, hx.data(), n4, cudaMemcpyHostToDevice));
    // FP32 position copies (ping-pong with d_x0/d_x1); WCM_BD_FP32_GATHER_OFF=1 or WCM_FP64=1 = FP64 gathers (bd_step_kernel)
    static const int u_bd_fp32_gather = std::getenv("WCM_BD_FP32_GATHER_U") ? atoi(std::getenv("WCM_BD_FP32_GATHER_U")) : 4;
    const bool use_bd_fp32_gather = std::getenv("WCM_BD_FP32_GATHER_OFF") == nullptr && std::getenv("WCM_FP64") == nullptr;
    float4 *d_xf0 = nullptr, *d_xf1 = nullptr;
    if (use_bd_fp32_gather) {
        std::vector<float4> hf((size_t)N);
        for (int i = 0; i < N; i++) hf[i] = make_float4((float)hx[i].x, (float)hx[i].y, (float)hx[i].z, (float)hx[i].w);
        BD_CHECK(cudaMalloc(&d_xf0, (size_t)N * sizeof(float4))); BD_CHECK(cudaMalloc(&d_xf1, (size_t)N * sizeof(float4)));
        BD_CHECK(cudaMemcpy(d_xf0, hf.data(), (size_t)N * sizeof(float4), cudaMemcpyHostToDevice));
    }
    BD_CHECK(cudaMemcpy(d_type, s->type, N * sizeof(int), cudaMemcpyHostToDevice));
    BD_CHECK(cudaMemcpy(d_mob, s->mobile, N, cudaMemcpyHostToDevice));
    BD_CHECK(cudaMemcpy(d_boff, boff.data(), (N + 1) * sizeof(int), cudaMemcpyHostToDevice));
    if (!bdat.empty()) BD_CHECK(cudaMemcpy(d_bdat, bdat.data(), bdat.size() * sizeof(int), cudaMemcpyHostToDevice));
    BD_CHECK(cudaMemcpy(d_aoff, aoff.data(), (N + 1) * sizeof(int), cudaMemcpyHostToDevice));
    if (!adat.empty()) BD_CHECK(cudaMemcpy(d_adat, adat.data(), adat.size() * sizeof(int), cudaMemcpyHostToDevice));
    BD_CHECK(cudaMemcpy(d_eoff, s->excl_off, (N + 1) * sizeof(int), cudaMemcpyHostToDevice));
    if (ne > 0) BD_CHECK(cudaMemcpy(d_elist, s->excl_list, ne * sizeof(int), cudaMemcpyHostToDevice));

    // cell grid over the box
    double L[3], cs[3]; int nc[3];
    for (int d = 0; d < 3; d++) { L[d] = s->box_hi[d] - s->box_lo[d]; nc[d] = std::max(1, (int)(L[d] / cutmax)); cs[d] = L[d] / nc[d]; }
    const int ncells = nc[0] * nc[1] * nc[2];
    BD_CHECK(cudaMalloc(&d_start, (ncells + 1) * sizeof(int))); BD_CHECK(cudaMalloc(&d_count, (ncells + 1) * sizeof(int)));
    void *d_tmp = nullptr; size_t tmp_bytes = 0;
    BD_CHECK(cub::DeviceScan::ExclusiveSum(nullptr, tmp_bytes, d_count, d_start, ncells + 1));
    BD_CHECK(cudaMalloc(&d_tmp, std::max<size_t>(1, tmp_bytes)));

    const int T = 256, B = (N + T - 1) / T;
    int max_neighs = 256;
    // async_ovf != nullptr: no host copy of the overflow count; bd_ovf_kernel halts the steps instead (see the loop)
    int *ovf_halt = nullptr; volatile long long *ovf_info = nullptr; long long ovf_resume = 0;
    auto build = [&](const double4 *dx) -> int {
        for (;;) {
            if (d_nl == nullptr) BD_CHECK(cudaMalloc(&d_nl, (size_t)max_neighs * N * sizeof(int)));
            BD_CHECK(cudaMemsetAsync(d_count, 0, (ncells + 1) * sizeof(int)));
            bd_cell_kernel<<<B, T>>>(dx, N, nc[0], nc[1], nc[2], cs[0], cs[1], cs[2], s->box_lo[0], s->box_lo[1], s->box_lo[2], d_cellof, d_count);
            BD_CHECK(cub::DeviceScan::ExclusiveSum(d_tmp, tmp_bytes, d_count, d_start, ncells + 1));
            BD_CHECK(cudaMemsetAsync(d_count, 0, (ncells + 1) * sizeof(int)));
            bd_scatter_kernel<<<B, T>>>(d_cellof, N, d_start, d_count, d_sorted);
            BD_CHECK(cudaMemsetAsync(d_ovf, 0, sizeof(int)));
            // warp-per-atom build (same lists, same order, same overflow); WCM_BD_WARP_NLIST_OFF=1: the thread-per-atom kernel
            static const bool off_bd_warp_nlist = std::getenv("WCM_BD_WARP_NLIST_OFF") != nullptr, verify_bd_warp_nlist = std::getenv("WCM_BD_WARP_NLIST_VERIFY") != nullptr;
            if (off_bd_warp_nlist || verify_bd_warp_nlist)
                bd_build_kernel<<<B, T>>>(dx, d_type, d_mob, N, d_sorted, d_start, nc[0], nc[1], nc[2], cs[0], cs[1], cs[2],
                                          s->box_lo[0], s->box_lo[1], s->box_lo[2], d_eoff, d_elist, d_nn, d_nl, max_neighs, N, d_ovf);
            if (!off_bd_warp_nlist) {
                static int *v_nn = nullptr, *v_nl = nullptr, *v_ovf = nullptr, *v_bad = nullptr; static size_t v_cap = 0;
                int *nn_out = d_nn, *nl_out = d_nl, *ovf_out = d_ovf;
                if (verify_bd_warp_nlist) {
                    const size_t need = (size_t)max_neighs * N;
                    if (need > v_cap) {
                        if (v_nl) { cudaFree(v_nl); cudaFree(v_nn); }
                        BD_CHECK(cudaMalloc(&v_nl, need * sizeof(int))); BD_CHECK(cudaMalloc(&v_nn, (size_t)N * sizeof(int)));
                        if (!v_ovf) { BD_CHECK(cudaMalloc(&v_ovf, sizeof(int))); BD_CHECK(cudaMalloc(&v_bad, sizeof(int))); }
                        v_cap = need;
                    }
                    BD_CHECK(cudaMemsetAsync(v_ovf, 0, sizeof(int)));
                    nn_out = v_nn; nl_out = v_nl; ovf_out = v_ovf;
                }
                const long long wthreads = 32LL * N;
                bd_build_warp_kernel<<<(unsigned)((wthreads + T - 1) / T), T>>>(dx, d_type, d_mob, N, d_sorted, d_start, nc[0], nc[1], nc[2],
                    cs[0], cs[1], cs[2], s->box_lo[0], s->box_lo[1], s->box_lo[2], d_eoff, d_elist, nn_out, nl_out, max_neighs, N, ovf_out);
                if (verify_bd_warp_nlist) {
                    BD_CHECK(cudaMemsetAsync(v_bad, 0, sizeof(int)));
                    bd_nlist_compare_kernel<<<B, T>>>(N, d_nn, d_nl, v_nn, v_nl, N, v_bad);
                    int bad = 0, o1 = 0, o2 = 0;
                    BD_CHECK(cudaMemcpy(&bad, v_bad, sizeof(int), cudaMemcpyDeviceToHost));
                    BD_CHECK(cudaMemcpy(&o1, d_ovf, sizeof(int), cudaMemcpyDeviceToHost)); BD_CHECK(cudaMemcpy(&o2, v_ovf, sizeof(int), cudaMemcpyDeviceToHost));
                    static long builds = 0, bad_builds = 0; builds++; if (bad || o1 != o2) bad_builds++;
                    if (bad || o1 != o2 || builds % 200 == 1)
                        fprintf(stderr, "WCM_BD_WARP_NLIST_VERIFY: build %ld: %d of %d rows differ, overflow %d / %d (%ld of %ld builds with differences)\n",
                                builds, bad, N, o1, o2, bad_builds, builds);
                }
            }
            if (ovf_halt) { bd_ovf_kernel<<<1, 32>>>(d_ovf, ovf_halt, ovf_resume, ovf_info); break; }
            int ovf = 0;
            BD_CHECK(cudaMemcpy(&ovf, d_ovf, sizeof(int), cudaMemcpyDeviceToHost));
            if (ovf == 0) break;
            cudaFree(d_nl); d_nl = nullptr; max_neighs = ovf + ovf / 4 + 16;   // grow and rebuild
        }
        BD_CHECK(cudaMemcpyAsync(d_xb, dx, n4, cudaMemcpyDeviceToDevice));
        return 0;
    };

    double build_ms = 0.0;
    auto tb = clk::now();
    if (build(d_x0) != 0) return -1;
    BD_CHECK(cudaDeviceSynchronize());
    build_ms += std::chrono::duration<double, std::milli>(clk::now() - tb).count();
    st->builds = 1;

    BDStepArgs A;
    A.type = d_type; A.mobile = d_mob; A.N = N; A.numneigh = d_nn; A.nl = d_nl; A.stride = N;
    A.boff = d_boff; A.bdat = d_bdat; A.aoff = d_aoff; A.adat = d_adat; A.xbuild = d_xb;
    A.dt = s->dt; A.g1 = s->g1; A.g2 = s->g2;
    A.key = make_uint2((uint32_t)(s->seed & 0xffffffffu), (uint32_t)(s->seed >> 32));
    A.maxd2_bits = d_maxd2;
    // Rebuild when an atom has moved skin/2 - 6 A since the last build: the displacement is only looked at every check_every
    // steps, and a pair can be missed only if two atoms together moved more than the skin; 6 A covers the motion between checks
    // (probes at t=2996, checks every 5 steps: at most 3.9-5.0 A past the trigger for skins 30-65 A). the fused-bd-skin change (was 0.4 x skin).
    const double trig = std::max(0.25 * s->skin, 0.5 * s->skin - 6.0);
    const double half_skin2 = trig * trig;
    double4 *xin = d_x0, *xout = d_x1;
    float4 *xinf = d_xf0, *xoutf = d_xf1;
    BD_CHECK(cudaMemset(d_maxd2, 0, sizeof(unsigned long long)));
    auto tk = clk::now();
    // pipelined steps with the rebuild decision on the device (use_bd_fp32_gather path; WCM_DEVICE_LOOPS_OFF=1 or WCM_DEVICE_LOOPS_BD_OFF=1: the loop below). The host
    // enqueues steps in chunks without waiting; each check step is followed by bd_check_kernel. When a rebuild is due the steps
    // after the check do nothing; the host sees the mapped flag one chunk later, waits for the (empty) queue, builds the lists from
    // the positions of the resume step (ping-pong parity of the step index) and continues from that step. The executed steps,
    // their inputs and the rebuild points are those of the loop below.
    static const bool off_device_loops = std::getenv("WCM_DEVICE_LOOPS_OFF") != nullptr || std::getenv("WCM_DEVICE_LOOPS_BD_OFF") != nullptr;
    static const int ch_device_loops = std::getenv("WCM_DEVICE_LOOPS_CHUNK") ? std::max(1, atoi(std::getenv("WCM_DEVICE_LOOPS_CHUNK"))) : 4;
    if (use_bd_fp32_gather && !off_device_loops && (u_bd_fp32_gather == 2 || u_bd_fp32_gather == 4 || u_bd_fp32_gather == 8)) {
        static int *d_halt = nullptr; static long long *h_info = nullptr, *d_info = nullptr;
        if (!d_halt) {
            BD_CHECK(cudaMalloc(&d_halt, sizeof(int)));
            BD_CHECK(cudaHostAlloc((void **)&h_info, 3 * sizeof(long long), cudaHostAllocMapped));
            BD_CHECK(cudaHostGetDevicePointer((void **)&d_info, h_info, 0));
        }
        BD_CHECK(cudaMemset(d_halt, 0, sizeof(int)));
        ((volatile long long *)h_info)[0] = -1; ((volatile long long *)h_info)[2] = 0;
        A.halt = d_halt;
        ovf_halt = d_halt; ovf_info = (volatile long long *)d_info;
        // probe: WCM_DEVICE_LOOPS_OVF_TEST=1 quarters the list capacity after the first build, so the next build overflows and takes the
        // device-side overflow path (the list in use stays valid: it is only freed when that path grows it)
        if (std::getenv("WCM_DEVICE_LOOPS_OVF_TEST")) max_neighs = std::max(8, max_neighs / 4);
        cudaEvent_t ev[2];
        BD_CHECK(cudaEventCreateWithFlags(&ev[0], cudaEventDisableTiming));
        BD_CHECK(cudaEventCreateWithFlags(&ev[1], cudaEventDisableTiming));
        long step = 0; long chunk = 0;
        for (;;) {
            if (step < nsteps) {
                for (int c = 0; c < ch_device_loops && step < nsteps; c++, step++) {
                    const bool check = s->check_every > 0 && ((step + 1) % s->check_every == 0) && step + 1 < nsteps;
                    A.xin = (step & 1) ? d_x1 : d_x0; A.xout = (step & 1) ? d_x0 : d_x1;
                    A.xinf = (step & 1) ? d_xf1 : d_xf0; A.xoutf = (step & 1) ? d_xf0 : d_xf1;
                    A.nl = d_nl;
                    A.step_lo = (uint32_t)(step & 0xffffffff); A.step_hi = (uint32_t)((unsigned long long)step >> 32);
                    A.do_check = check ? 1 : 0;
                    if (u_bd_fp32_gather == 8) bd_launch(bd_step_kernel_f<8, true>, dim3(B), dim3(T), A);
                    else if (u_bd_fp32_gather == 2) bd_launch(bd_step_kernel_f<2, true>, dim3(B), dim3(T), A);
                    else bd_launch(bd_step_kernel_f<4, true>, dim3(B), dim3(T), A);
                    if (check) bd_launch(bd_check_kernel, dim3(1), dim3(32), d_maxd2, half_skin2, d_halt, (long long)(step + 1), (volatile long long *)d_info);
                }
                BD_CHECK(cudaEventRecord(ev[chunk & 1]));
                if (chunk > 0) BD_CHECK(cudaEventSynchronize(ev[(chunk - 1) & 1]));
                chunk++;
            } else {
                BD_CHECK(cudaDeviceSynchronize());
            }
            long long rs = ((volatile long long *)h_info)[0];
            if (rs < 0) { if (step >= nsteps) break; continue; }
            BD_CHECK(cudaDeviceSynchronize());   // the steps enqueued after the halt have run (as no-ops)
            rs = ((volatile long long *)h_info)[0];
            const long long ovf = ((volatile long long *)h_info)[2];
            if (ovf > 0) {   // the last build overflowed: grow and build again from the same positions (the build loop's retry)
                if (std::getenv("WCM_DEVICE_LOOPS_OVF_TEST")) fprintf(stderr, "list overflow (%lld > %d) handled on the device path at step %lld\n", ovf, max_neighs, rs);
                cudaFree(d_nl); d_nl = nullptr; max_neighs = (int)ovf + (int)ovf / 4 + 16;
                ((volatile long long *)h_info)[2] = 0;
            } else {
                unsigned long long bits = (unsigned long long)((volatile long long *)h_info)[1];
                double md2; memcpy(&md2, &bits, sizeof(md2));
                st->max_disp_at_build = std::max(st->max_disp_at_build, std::sqrt(md2));
                st->builds++;
            }
            auto t0 = clk::now();
            BD_CHECK(cudaMemsetAsync(d_halt, 0, sizeof(int)));
            ((volatile long long *)h_info)[0] = -1;
            ovf_resume = rs;
            if (build((rs & 1) ? d_x1 : d_x0) != 0) return -1;
            build_ms += std::chrono::duration<double, std::milli>(clk::now() - t0).count();
            step = rs; chunk = 0;
        }
        cudaEventDestroy(ev[0]); cudaEventDestroy(ev[1]);
        ovf_halt = nullptr; ovf_info = nullptr;
        xin = (nsteps & 1) ? d_x1 : d_x0;   // the positions after the last step
    } else
    for (long step = 0; step < nsteps; step++) {
        const bool check = s->check_every > 0 && ((step + 1) % s->check_every == 0) && step + 1 < nsteps;
        A.xin = xin; A.xout = xout; A.nl = d_nl;
        A.step_lo = (uint32_t)(step & 0xffffffff); A.step_hi = (uint32_t)((unsigned long long)step >> 32);
        A.do_check = check ? 1 : 0;
        if (use_bd_fp32_gather) {
            A.xinf = xinf; A.xoutf = xoutf;
            if (u_bd_fp32_gather == 8) bd_launch(bd_step_kernel_f<8>, dim3(B), dim3(T), A); else if (u_bd_fp32_gather == 2) bd_launch(bd_step_kernel_f<2>, dim3(B), dim3(T), A); else bd_launch(bd_step_kernel_f<4>, dim3(B), dim3(T), A);
            std::swap(xinf, xoutf);
        } else {
            bd_step_kernel<<<B, T>>>(A);
        }
        std::swap(xin, xout);
        if (check) {
            unsigned long long bits = 0;
            BD_CHECK(cudaMemcpy(&bits, d_maxd2, sizeof(bits), cudaMemcpyDeviceToHost));   // syncs with the step
            double md2; memcpy(&md2, &bits, sizeof(md2));
            if (md2 > half_skin2) {
                st->max_disp_at_build = std::max(st->max_disp_at_build, std::sqrt(md2));
                auto t0 = clk::now();
                if (build(xin) != 0) return -1;
                BD_CHECK(cudaDeviceSynchronize());
                build_ms += std::chrono::duration<double, std::milli>(clk::now() - t0).count();
                st->builds++;
            }
            BD_CHECK(cudaMemsetAsync(d_maxd2, 0, sizeof(unsigned long long)));
        }
    }
    BD_CHECK(cudaGetLastError());
    BD_CHECK(cudaDeviceSynchronize());
    st->kernel_ms = std::chrono::duration<double, std::milli>(clk::now() - tk).count() - build_ms;
    st->build_ms = build_ms;
    BD_CHECK(cudaMemcpy(hx.data(), xin, n4, cudaMemcpyDeviceToHost));
    for (int i = 0; i < N; i++) { s->x[3*i] = hx[i].x; s->x[3*i+1] = hx[i].y; s->x[3*i+2] = hx[i].z; }

    if (d_xf0) { cudaFree(d_xf0); cudaFree(d_xf1); }
    cudaFree(d_x0); cudaFree(d_x1); cudaFree(d_xb); cudaFree(d_type); cudaFree(d_mob); cudaFree(d_nn); cudaFree(d_nl);
    cudaFree(d_boff); cudaFree(d_bdat); cudaFree(d_aoff); cudaFree(d_adat); cudaFree(d_eoff); cudaFree(d_elist);
    cudaFree(d_cellof); cudaFree(d_sorted); cudaFree(d_start); cudaFree(d_count); cudaFree(d_tmp); cudaFree(d_ovf); cudaFree(d_maxd2);
    st->total_ms = std::chrono::duration<double, std::milli>(clk::now() - t_all).count();
    return 0;
}
