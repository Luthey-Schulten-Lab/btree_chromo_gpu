// Author(s): Ron Acda (using an iterative LLM-guided workflow, https://github.com/quarkron/iterative-hillclimber/tree/main)
#pragma once
// warp-per-atom two-grid neighbour build for the fused CG minimiser, producing the same lists as w_build_nlist2_kernel
// (same neighbours, same order in every row, same numneigh and the same max_neighs truncation).
//
// w_build_nlist2_kernel gives each atom one thread that walks its 27 fine cells (DNA partners) and then its 27 coarse cells
// (boundary / ribosome partners) candidate by candidate. A DNA bead next to the membrane has ~350 boundary candidates in its
// coarse cells, so a few thousand such threads set the whole kernel (275 us per build at t = 2996 s, 1 387 builds per hook).
// Here one warp takes one atom: the cells of a stencil are taken 32 at a time, a warp prefix sum of their occupancies flattens
// their candidates into one index range in the original order (cell order dx, dy, dz; within a cell the sorted order), the
// lanes test 32 consecutive candidates at once with the identical acceptance test, and a ballot + popcount gives each accepted
// candidate its slot, so the row is written in exactly the order the serial loop writes it.
// Needs is_excluded() and d_nl_cut2 from fused_cg_minimize.cu (included after them).

__device__ __forceinline__ void w3_scan_stencil(
    const double *x, const int *type, int i, int ti, double xi, double yi, double zi,
    const int *sorted, const int *cstart, int ix, int iy, int iz, int S, int ncx, int ncy, int ncz,
    double rcut2, const int *excl_list, int e_start, int e_end,
    int *neighlist, int max_neighs, int nlist_stride, int lane, int &nn)
{
    const int W = 2 * S + 1, ncell = W * W * W;
    for (int g = 0; g < ncell; g += 32) {
        // this lane's cell of the group (the serial loop's order: dx outer, dy, dz inner), its occupancy and first slot
        const int c = g + lane;
        int cnt = 0, first = 0;
        if (c < ncell) {
            const int dx = c / (W * W) - S, dy = (c / W) % W - S, dz = c % W - S;
            const int jx = ix + dx, jy = iy + dy, jz = iz + dz;
            if (jx >= 0 && jx < ncx && jy >= 0 && jy < ncy && jz >= 0 && jz < ncz) {
                const int cidx = jx * ncy * ncz + jy * ncz + jz;
                first = cstart[cidx]; cnt = cstart[cidx + 1] - first;
            }
        }
        // inclusive prefix of the occupancies over the group
        int incl = cnt;
        #pragma unroll
        for (int o = 1; o < 32; o <<= 1) { const int v = __shfl_up_sync(0xffffffffu, incl, o); if (lane >= o) incl += v; }
        const int total = __shfl_sync(0xffffffffu, incl, 31);
        for (int base = 0; base < total; base += 32) {
            const int k = base + lane;
            // the group cell holding flat candidate k: the first lane whose inclusive prefix exceeds k (every lane searches,
            // past-the-end lanes with k clamped, so that the shuffles run on the full warp)
            const int kk = min(k, total - 1);
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
                    const double ddx = xi - x[3*j], ddy = yi - x[3*j+1], ddz = zi - x[3*j+2];
                    const double d2 = ddx*ddx + ddy*ddy + ddz*ddz;
                    ok = d2 < rcut2 && d2 < d_nl_cut2[ti][type[j]] && !is_excluded(excl_list, e_start, e_end, j);
                }
            }
            const unsigned m = __ballot_sync(0xffffffffu, ok);
            const int pos = nn + __popc(m & ((1u << lane) - 1u));
            if (ok && pos < max_neighs) neighlist[(size_t)pos * nlist_stride + i] = j;
            nn += __popc(m);
        }
    }
}

__global__ void w_build_nlist3_kernel(
    const double *x, const int *type, int N, double blo_x, double blo_y, double blo_z, double rcut2,
    const int *f_sorted, const int *f_start, int nfx, int nfy, int nfz, double fcx, double fcy, double fcz, int S_ribo,
    const int *c_sorted, const int *c_start, int ncx, int ncy, int ncz, double cx, double cy, double cz,
    const int *excl_off, const int *excl_list, int *numneigh, int *neighlist, int max_neighs, int nlist_stride)
{
    const int i = (int)((blockIdx.x * (size_t)blockDim.x + threadIdx.x) >> 5), lane = threadIdx.x & 31;
    if (i >= N) return;                                  // whole warps leave together (i is uniform in a warp)
    const int ti = type[i];
    if (ti == 1) { if (lane == 0) numneigh[i] = 0; return; }
    const double xi = x[3*i], yi = x[3*i+1], zi = x[3*i+2];
    const int e_start = excl_off[i], e_end = excl_off[i + 1];
    int nn = 0;
    {
        const int ix = min(max((int)((xi - blo_x) / fcx), 0), nfx - 1);
        const int iy = min(max((int)((yi - blo_y) / fcy), 0), nfy - 1);
        const int iz = min(max((int)((zi - blo_z) / fcz), 0), nfz - 1);
        w3_scan_stencil(x, type, i, ti, xi, yi, zi, f_sorted, f_start, ix, iy, iz, (ti >= 3) ? 1 : S_ribo, nfx, nfy, nfz,
                        rcut2, excl_list, e_start, e_end, neighlist, max_neighs, nlist_stride, lane, nn);
    }
    {
        const int ix = min(max((int)((xi - blo_x) / cx), 0), ncx - 1);
        const int iy = min(max((int)((yi - blo_y) / cy), 0), ncy - 1);
        const int iz = min(max((int)((zi - blo_z) / cz), 0), ncz - 1);
        w3_scan_stencil(x, type, i, ti, xi, yi, zi, c_sorted, c_start, ix, iy, iz, 1, ncx, ncy, ncz,
                        rcut2, excl_list, e_start, e_end, neighlist, max_neighs, nlist_stride, lane, nn);
    }
    if (lane == 0) numneigh[i] = (nn <= max_neighs) ? nn : max_neighs;
}

// WCM_CG_WARP_NLIST_VERIFY: rows that differ between two builds (numneigh, or any of the first numneigh entries)
__global__ void w_nlist_compare_kernel(int N, const int *nn_a, const int *nl_a, const int *nn_b, const int *nl_b, int stride, int *bad)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    bool diff = nn_a[i] != nn_b[i];
    for (int k = 0; !diff && k < nn_a[i]; k++) diff = nl_a[(size_t)k * stride + i] != nl_b[(size_t)k * stride + i];
    if (diff) atomicAdd(bad, 1);
}
