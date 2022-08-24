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

## Support
brg4@illinois.edu

## Authors and acknowledgment
Benjamin R. Gilbert - brg4@illinois.edu

## Project status
This project is under development.
