# btree_chromo

## Description
Program to model theta structures of replicating bacterial chromosomes.

## Installation

**make all**

The executable (*program*) will be in /build/apps.

## Usage

Prepare a *directives.inp* file containing the directives to be executed by the binary tree program. Lines beginning with '#' are ignoreed.

Prepare any input files needed for the chosen directives.

Run with: **./program (some location)/directives.inp**

Use the testcase for an example: **./program /home/ben/Workspace/btree_chromo/test_case/directives.inp**

*Possible Directives (parameters are comma-separated and following ':' when needed)*
 - input_state:input_file - *creates state from input_file*
 - print - *prints binary tree state in terminal*
 - btree_prng_seed:seed - *seeds the btree's prng*
 - replicator_prng_seed:seed - *seeds the replicator's prng*
 - transform:(b)\_(r_cw)\_(r_ccw) - *applies single transform to branch (b), with replication extents (r\_cw) and (r\_ccw)*
 - transforms_file:transforms_file - *applies transforms stored in transforms_file*
 - random_transforms:N - *applies random transforms until (N) units are added or the maximum size is reached*
 - load_replication_model:rep_model_file - *loads the parameters and initial conditions for a replication model*
 - replicate:t - *replicates the current btree using loaded replication model for (t) seconds*
 - output_state:output_file - *writes state to output_file*
 - update_topology - *solves bond topology of system*
 - dump_topology:topology_file,idx - *dumps topology to topology_file with selected indexing convention (idx)*
 - update_CG_map:f_CG - *update coarse-graining with selected factor (f_CG)*
 - dump_CG_map:CG_map_file,f_CG,idx - *dumps CG_map to CG_map_file with selected factor (f_CG) and indexing convention (idx)*
 - regions_file:regions_file,idx - *reads chromosome regions from regions_file with selected indexing convention (idx)*
 - dump_regions:regions_count_file - *updates regions counts given current state and dumps counts to regions_count_file*
 - load_BD_lengths:BD_length_file - *reads lengths for Brownian dynamics simulation*
 - load_mono_coords:coords_file,order - *reads binary file with monomer coordinates (doubles) using data ordering convention (row/col)*
 - load_ribo_coords:coords_file,order - *reads binary file with ribosome coordinates (doubles) using data ordering convention (row/col)*
 - write_LAMMPS_data:LAMMPS_data_file - *write a LAMMPS file (data.*) using the current mono, ribo, and bdry coordinates, and the current replication state for the bond/angle topology*
 - write_mono_xyz:mono_file_xyz - *write the current monomer coordinates as an .xyz file to load into VMD*
 - set_initial_state - *set the initial state of the mapper to the current replication state*
 - set_final_state - *set the final state of the mapper to the current replication state*
 - map_replication - *based on the difference in final and initial replication states, determine new monomer coordinates and add new monomers to the LAMMPS system*

## Support
brg4@illinois.edu

## Authors and acknowledgment
Benjamin R. Gilbert - brg4@illinois.edu

## Project status
This project is under development.
