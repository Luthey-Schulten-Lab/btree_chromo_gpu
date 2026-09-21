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

## Parameters relating to SMC behavior

Following Maytin *et al.* 2026, three parameters control loop extrusion, and
tuning them lets chromosome segregation happen without any external force. Their
product divided by the chromosome length is roughly the fraction of the
chromosome extruded into loops, and needs to be greater than ~1 for segregation
to be reliable.

| Parameter | Default | Set by |
|---|---|---|
| Translocation speed | 500 bp/s (2000 bp per 4 s hook) | `translocate:100,T` |
| Dwell time | ~200 s | `basal_death_prob=0.0002` |
| Active SMC complexes | ~50 at *t*=0 | `numSmc` (from the RDME Smc count) |

**Translocation speed** — how fast the complexes extrude loops, here 500 bp/s.
At the start of each 4 s hook both sides of the SMC translocate by 100 beads in
the `loop_simulator` object, written `translocate:100,T` in the chromosome
operations file. The LAMMPS system is then re-initialized and the
`loop_simulator` state written into it with `simulator_form_loops:F`; the `F`
means bond only the SMC heads together, without adding extra bonds inside the
loop.

**Dwell time** — how long a complex stays bound on average. This is set as the
probability that the SMC unbinds each time the `loop_simulator` advances one
step. At the default `basal_death_prob=0.0002` it takes ~5000 steps to unbind,
a dwell time of 200 s.

**Number of active SMC** — how many complexes are bound and actively extruding.
The total comes from the Smc protein count halved, since Smc is a homodimer,
then scaled by the bound fraction: `int((P_0415 / 2) * dna_smc_bound_fraction)`,
with `dna_smc_bound_fraction=0.5` by default. Proteomics gives 200 Smc proteins
at the start, so 100 dimers and 50 active complexes. That roughly doubles to 100
by the end of the cell cycle, though growth varies between replicates — in one
replicate with little *0415* expression it rose only from 50 to 66 after
105 minutes.

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
