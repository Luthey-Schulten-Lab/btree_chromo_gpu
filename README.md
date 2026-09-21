# btree_chromo v2.0.0

GPU chromosome engine (C++ / Kokkos–LAMMPS) for replicating bacterial DNA.

**This is version 2.** It inherits the persistent SMC loop-extrusion model from
the `protein_science` branch (Maytin *et al.*, *Protein Science*, 2026) and adds
optimizations such as fused minimizer and related fixes. 

## Which branch should I use?

| Branch | What it is | Use it when |
|---|---|---|
| **`btree_chromo_v2.0.0`** | This branch. v1 physics + GPU optimizations + fixes. | You want to run 4DWCM `1.1.0`. |
| `protein_science` | Maytin *et al.* as published, frozen at `4635d2a`. | You are reproducing that paper. |
| `main` | btree_chromo version 1 | You want to run 4DWCM `main`. |

Inherited from `protein_science`: persistent SMC loop state across DNA hooks via
`load_loops` / `write_loops` / `translocate`, which lets SMC dwell times run much
longer than a single hook interval (~4 s).

## What v2.0.0 adds

### Performance

**Fused GPU conjugate-gradient minimizer** (`ce1d839`) —
`src/fused_cg_minimize.cu` replaces the stock LAMMPS CPU minimize for the
loop-bond relaxation run each hook. The Makefile compiles `src/*.cu` with the
rest of the tree.

## Build

Follow the lab Kokkos / LAMMPS DNA model install for this tree, then:

```bash
git clone -b btree_chromo_v2.0.0 \
  https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu.git
cd btree_chromo_gpu
make
```

## Whole-cell integration

For 4DWCM, the wrapper, pinned SHAs, and container build are in
[Minimal_Cell_4DWCM](https://github.com/Luthey-Schulten-Lab/Minimal_Cell_4DWCM/tree/1.1.0).
