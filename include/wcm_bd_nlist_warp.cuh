// Author(s): Ron Acda (using an iterative LLM-guided workflow, https://github.com/quarkron/iterative-hillclimber/tree/main)
#pragma once
// warp-per-atom neighbour build for the fused BD, producing the same lists as bd_build_kernel (same neighbours, same
// order in every row, same numneigh, same max_neighs truncation and the same overflow maximum).
//
// bd_build_kernel walks each mobile atom's 27 cells candidate by candidate in one thread (392 us per build at t = 2996 s, 737
// builds per hook). As in the cg-warp-nlist change (wcm_nlist_warp.cuh) one warp takes one atom: the 27 cells' occupancies are prefix-summed
// across the warp, their candidates are tested 32 at a time in the serial order (cell order dx, dy, dz; within a cell the sorted
// order) with the identical test, and a ballot + popcount places each accepted candidate in its serial slot.
// Needs c_bd_nl2 and bd_excluded() from fused_bd.cu (included after them).

__global__ void bd_build_warp_kernel(const double4 *x, const int *type, const unsigned char *mobile, int N,
                                     const int *sorted, const int *start, int ncx, int ncy, int ncz, double cx, double cy, double cz,
                                     double lx, double ly, double lz, const int *eoff, const int *elist,
                                     int *numneigh, int *nl, int max_neighs, int stride, int *overflow)
{
    const int i = (int)((blockIdx.x * (size_t)blockDim.x + threadIdx.x) >> 5), lane = threadIdx.x & 31;
    if (i >= N) return;                                  // uniform in a warp
    if (!mobile[i]) { if (lane == 0) numneigh[i] = 0; return; }
    const int ti = type[i];
    const double4 qi = x[i];
    const double xi = qi.x, yi = qi.y, zi = qi.z;
    const int ix = min(max((int)((xi - lx) / cx), 0), ncx - 1);
    const int iy = min(max((int)((yi - ly) / cy), 0), ncy - 1);
    const int iz = min(max((int)((zi - lz) / cz), 0), ncz - 1);
    const int es = eoff[i], ee = eoff[i + 1];
    // lane c < 27: cell c of the stencil in the serial order (dx outer, dy, dz inner)
    int cnt = 0, first = 0;
    if (lane < 27) {
        const int jx = ix + lane / 9 - 1, jy = iy + (lane / 3) % 3 - 1, jz = iz + lane % 3 - 1;
        if (jx >= 0 && jx < ncx && jy >= 0 && jy < ncy && jz >= 0 && jz < ncz) {
            const int c = (jx * ncy + jy) * ncz + jz;
            first = start[c]; cnt = start[c + 1] - first;
        }
    }
    int incl = cnt;
    #pragma unroll
    for (int o = 1; o < 32; o <<= 1) { const int v = __shfl_up_sync(0xffffffffu, incl, o); if (lane >= o) incl += v; }
    const int total = __shfl_sync(0xffffffffu, incl, 31);
    int nn = 0;
    for (int base = 0; base < total; base += 32) {
        const int k = base + lane, kk = min(k, total - 1);
        int lo = 0, hi = 31;
        #pragma unroll
        for (int step = 0; step < 5; step++) {
            const int mid = (lo + hi) >> 1;
            const int pm = __shfl_sync(0xffffffffu, incl, mid);
            if (pm > kk) hi = mid; else lo = mid + 1;
        }
        const int pc = __shfl_sync(0xffffffffu, incl, lo), cc = __shfl_sync(0xffffffffu, cnt, lo), fc = __shfl_sync(0xffffffffu, first, lo);
        bool ok = false; int j = -1;
        if (k < total) {
            j = sorted[fc + (k - (pc - cc))];
            if (j != i) {
                const double4 qj = x[j];
                const double nl2 = c_bd_nl2[ti][(int)qj.w];
                if (nl2 > 0.0) {
                    const double ddx = xi - qj.x, ddy = yi - qj.y, ddz = zi - qj.z;
                    ok = ddx*ddx + ddy*ddy + ddz*ddz < nl2 && !bd_excluded(elist, es, ee, j);
                }
            }
        }
        const unsigned m = __ballot_sync(0xffffffffu, ok);
        const int pos = nn + __popc(m & ((1u << lane) - 1u));
        if (ok && pos < max_neighs) nl[(size_t)pos * stride + i] = j;
        nn += __popc(m);
    }
    if (lane == 0) {
        numneigh[i] = min(nn, max_neighs);
        if (nn > max_neighs) atomicMax(overflow, nn);
    }
}

// WCM_BD_WARP_NLIST_VERIFY: rows that differ between two builds (numneigh, or any of the first numneigh entries)
__global__ void bd_nlist_compare_kernel(int N, const int *nn_a, const int *nl_a, const int *nn_b, const int *nl_b, int stride, int *bad)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    bool diff = nn_a[i] != nn_b[i];
    for (int k = 0; !diff && k < nn_a[i]; k++) diff = nl_a[(size_t)k * stride + i] != nl_b[(size_t)k * stride + i];
    if (diff) atomicAdd(bad, 1);
}
