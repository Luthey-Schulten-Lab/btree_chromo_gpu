# btree_chromo v2.0.0

GPU chromosome engine (C++ / Kokkos–LAMMPS) for replicating bacterial DNA.

**This is version 2.** It inherits the persistent SMC loop-extrusion model from
the `protein_science` branch (Maytin *et al.*, *Protein Science*, 2026) and adds
GPU performance work plus three correctness fixes. The physics of loop extrusion
is unchanged; what is new is how fast it runs and how reliably it writes output.

## Which branch should I use?

| Branch | What it is | Use it when |
|---|---|---|
| **`btree_chromo_v2.0.0`** | This branch. v1 physics + GPU optimizations + fixes. | You want the current engine. |
| `protein_science` | Maytin *et al.* as published, frozen at `4635d2a`. | You are reproducing that paper. |
| `main` | Separate public lineage — **no shared history** with either branch above. | Not a drop-in replacement; do not mix. |

Inherited from `protein_science`: persistent SMC loop state across DNA hooks via
`load_loops` / `write_loops` / `translocate`, which lets SMC dwell times run much
longer than a single hook interval (~4 s).

## What v2.0.0 adds

### Performance

**In-place LAMMPS coordinate update** (`ce1d839`) — `sys_update_lammps_inplace`
and `scatter_coords_from_sys` scatter bead coordinates directly, removing the
second per-hook write / clear / `read_data` round-trip.

**Fused GPU conjugate-gradient minimizer** (`ce1d839`) —
`src/fused_cg_minimize.cu` replaces the stock LAMMPS CPU minimize for the
loop-bond relaxation run each hook. The Makefile compiles `src/*.cu` with the
rest of the tree.

### Correctness

**Hybrid loop bond fix** (`b4c0f9a`) — the `bond_style hybrid` reset left the
replication-fork bond (type 3) unset, aborting `simulator_run_soft_FENE` during
chromosome partitioning. Adds the missing `bond_coeff` to `lmp.loop_bond_hybrid`.

**Boundary radius for soft bdry–DNA pairs** (`7e21ce0`) — the topoisomerase pair
model sized the soft boundary–DNA cutoff with the ribosome radius
(`sigma_mono_ribo`, 117 Å) instead of the boundary radius (`sigma_mono_bdry`,
217 Å). Because `soft/kk` is finite at contact, the undersized cutoff left a
cheap path through the membrane and DNA leaked out during replication. The
weakest point on the shell rises from 13.8 to 47.7 kBT, and a directional
containment test on seed 17 goes from up to 797 escaped beads per frame to zero.
Replication is unaffected — both runs reach exactly 108676 beads.

**Atomic output publication** (`d5dcada`) — `write_bin` and `dump_topology`
opened the final path directly, publishing the filename before the contents. A
reader polling for those files could open one mid-write and get a truncated
chromosome that still looked valid. They now write to a temporary name and
`rename(2)` it into place, so the final name only ever resolves to a finished
file. Note this is for visibility to a concurrent reader, not crash durability —
there is no `fsync`.

## Reproducing the JCP runs

The ~26.8 h full-cycle production run in the JCP manuscript was built from
`1725d1b`, which has the two performance changes but **not** the boundary or
output fixes. To reproduce it exactly, build that commit. For new work, use the
branch tip.

## Build

Follow the lab Kokkos / LAMMPS DNA model install for this tree, then:

```bash
git clone -b btree_chromo_v2.0.0 \
  https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu.git
cd btree_chromo_gpu
make
```

Pin a SHA rather than a branch name in automated builds — both branches move.

## Whole-cell integration

For 4DWCM, `numSmc` comes from `loop_params.txt` rather than replication-length
scaling; that handling is already upstream here, so no separate patch step is
needed. The wrapper, pinned SHAs, and container build are in
[Minimal_Cell_4DWCM](https://github.com/Luthey-Schulten-Lab/Minimal_Cell_4DWCM/tree/1.1.0).
