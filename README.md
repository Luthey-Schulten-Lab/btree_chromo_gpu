# btree_chromo (`protein_science` branch)

C++ / Kokkos–LAMMPS chromosome engine for replicating bacterial DNA with
**persistent SMC loop state** across DNA hooks (`load_loops` / `write_loops` /
`translocate`). This branch is the chromosome backend for the 4DWCM
**protein_science** performance stack (~30 h full cycle after host-side
isolation). Science context: Maytin et al. 2026 / Minimal_Cell_ChromosomeSegregation.

> **Note:** GitHub `main` is a separate public lineage (no shared merge-base with
> this branch). Do not treat `main` as a drop-in replacement for
> `protein_science`.

## Recent changes (performance / WCM coupling)

Relative to the previous `protein_science` tip on origin:

1. **In-place LAMMPS coordinate update** — `sys_update_lammps_inplace` /
   `scatter_coords_from_sys` replaces the second per-hook
   write / clear / `read_data` round-trip when refreshing bead coordinates.
2. **Fused GPU conjugate-gradient minimizer** — `src/fused_cg_minimize.cu`
   (headers under `include/fused_cg_*.h`) replaces the stock LAMMPS CPU
   minimize for loop-bond relaxation each hook.
3. **Hybrid loop bond fix** — set `bond_coeff` for replication-fork bond type 3
   in `lmp.loop_bond_hybrid` so `simulator_run_soft_FENE` does not abort during
   chromosome partitioning.

Makefile builds the new CUDA sources with the rest of `src/`.

These BD-side changes are what the JCP protein_science ~29.4 h 4DWCM runs use,
together with host CPU isolation in Optimize_4DWCM_Minimal_Cell.

## Whole-cell integration

For 4DWCM, also apply the small WCM patches so `numSmc` comes from
`loop_params.txt` (not replication-length scaling):

- Optimize repo: [`btree_chromo_wcm/`](https://github.com/luthey-schulten-chemistry/Optimize_4DWCM_Minimal_Cell/tree/protein_science/btree_chromo_wcm)

Pinned SHAs and tags: see Optimize [`VERSIONS.md`](https://github.com/luthey-schulten-chemistry/Optimize_4DWCM_Minimal_Cell/blob/protein_science/VERSIONS.md).

## Build (sketch)

Follow the lab Kokkos / LAMMPS DNA model install for this tree, then:

```bash
git checkout protein_science
# optional: apply btree_chromo_wcm patches from Optimize
make   # or your usual build that compiles src/*.cu
```

## Related

- Public default engine (different history): [`main`](https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu/tree/main)
- 4DWCM performance wrapper: [Optimize_4DWCM_Minimal_Cell](https://github.com/luthey-schulten-chemistry/Optimize_4DWCM_Minimal_Cell) (`protein_science`)
