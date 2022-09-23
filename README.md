# btree_chromo

## Description

**C++ program to model theta structures of replicating bacterial chromosomes using binary trees.**

The program is organized in the following manner:

    Upon execution a **btree_driver** executes a series of directives stored in a file provided by the user (a full list of directives is provided in the 'Usage' section). The **btree_driver** contains the following objects.
    
       - btree - *a specialized binary tree class representing replicating circular dsDNA in nested theta structures that can only be manipulated using public member functions representing processes physically possible for circular dsDNA*
       	 - determines topology when system is represented as circular(/theta structure) polymers for arbitrary replication states
	 - counts genome features at nucleotide resolution for arbitrary replication states
	 - prepares coarse-graining into chromosomal loci for contact map calculations for arbitrary replication states
	 
       - replicator - *a class containing a stochastic chemical kinetics model for DNA replication initiation due to DnaA, which then simulates trajectories in the space of replication states using the Gillespie method*
       
       - LAMMPS_sys - *a class storing the spatial (and other) information of a system of replicating circular dsDNA (monomers of ellipsoidal), ribosomes (ellipsoidal particles), and boundary particles (point-like particles) for use in Brownian dynamics simulations using LAMMPS*
       
       - mapper - *a class that governs the creation of new DNA monomers in the spatial model given an 'initial replication state' and a 'final replication state' using the train-track model of replication*
       
       - LAMMPS_simulator - *a class that runs Brownian dynamics and accessory routines using LAMMPS to simulate the spatial model*

By combining a series of directives together a user could perform the following example protocol.

   1) Define a circular chromosome of a specific size.
   2) Read a file containing defined genome features.
   3) Load a model of DnaA replication initiation.
   4) Load the spatial model of a DNA structure in the polymer representation matching the specified size.
   5) Run the replicator to determine a final replication state at halfway through the cell cycle.
   6) Use the mapper to create new DNA monomers in the spatial model given the change in replication states.
   7) Simulate the spatial model using Brownian dynamics in LAMMPS.
   8) Output the counts of genome features.
   9) Output the coarse-graining for contact map calculations.

   The user now has A) a replication state halfway through the cell cycle, B) a matching spatial model, C) the counts of genome features in this replication state, and D) the coarse-graining necessary to calculate chromosome contact maps for replicating chromosomes.

## Repository directory structure

 - /src - *source files for program*
 - /include - *header files for program*
 - /test_case - *test cases to demonstrate program*
 - /LAMMPS_src_additions - *source files to be added to LAMMPS build for DNA polymer model*
 - /LAMMPS_DNA_model - *files to be included in LAMMPS simulations running DNA polymer model*
 - /LAMMPS_basic_input_scripts - *basic input scripts for LAMMPS simulations of DNA polymer model*

## Installation

   1) Install **OpenMPI** - *I built v4.1.4 from source using GCC-v12.1.0*
   2) Install **fmt Library** - *I built v9.1.0 from source using GCC-v12.1.0*
   3) Follow instructions in '/LAMMPS_src_additions' to make the additions to the LAMMPS source code
   4) Build and install the modified version of **LAMMPS** - *I built the 2022/02/17 release using GCC-v12.1.0*
   5) Edit the following variables in the Makefile to match your installations from the previous steps. - *I used environment variables to specify these local installations on my machine in a bash scripts, but you can type them in manually if you would prefer.*
      - **OPENMPI_LDFLAGS**
      - **FMT_LDFLAGS**
      - **LAMMPS_LDFLAGS**
      - **OPENMPI_INCLUDE**
      - **FMT_INCLUDE**
      - **LAMMPS_INCLUDE**

    6) **make all**

The executable (*program*) will be in /build/apps.

## Usage

Prepare a *directives.inp* file containing the directives to be executed by the binary tree program. Lines beginning with '#' are ignored.

Prepare any input files needed for the chosen directives.

Run with: **./program (some location)/directives.inp**

Use the testcase for an example: **./program /home/ben/Workspace/btree_chromo/test_case/directives.inp**

*Possible Directives (parameters are comma-separated and following ':' when needed)*
 - input_state:input_file - *creates state from input_file*
 - print - *prints binary tree state in terminal*
 - btree_prng_seed:seed - *seeds the btree's prng*
 - replicator_prng_seed:seed - *seeds the replicator's prng*
 - transform:(b)\_cw(r_cw)\_ccw(r_ccw) - *applies single transform to branch (b), with replication extents (r\_cw) and (r\_ccw) along clockwise and counter-clockwise directions, respectively*
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
 - write_LAMMPS_data:LAMMPS_data_file - *write a LAMMPS file (data.-) using the current mono, ribo, and bdry coordinates, and the current replication state for the bond/angle topology*
 - write_mono_xyz:mono_file_xyz - *write the current monomer coordinates as an .xyz file to load into VMD*
 - set_initial_state - *set the initial state of the mapper to the current replication state*
 - set_final_state - *set the final state of the mapper to the current replication state*
 - map_replication - *based on the difference in final and initial replication states, determine new monomer coordinates and add new monomers to the LAMMPS system*
 - prepare_simulator:log_file - *initialize MPI and a LAMMPS object that writes its output to log_file, all further commands with 'simulator' in their name will use this LAMMPS object*
 - simulator_run_file:run_file - *executes the 'include' command to run the LAMMPS commands stored in run_file*
 - sync_simulator_and_system - *copies the current simulation state to the LAMMPS_sys object used to control the topology*
 - clear_simulator - *execute the 'clear' command to clear the LAMMPS object*

## Support
brg4@illinois.edu

## Authors and acknowledgment
Benjamin R. Gilbert - brg4@illinois.edu

## Project status
This project is under development.
