# Diagnostics tool: diag

Post-processing tool for GKV binary output



### What is diag?
One difficulty of users may read GKV binary output which are decomposed by MPI. Post-processing program **diag** helps to read GKV binary output of a desired quantity at a desired time step. Since the read quantity is constructed as a global variable, e.g., *phi(-nx:nx,0:global ny,-global nz:global nz-1)* (not a local variable decomposed by MPI, *phi_local(-nx:nx,0:ny,-nz:nz-1)*), users do not need to be conscious of MPI parallelization of GKV. The main program *diag_main.f90* calls each diagnostics module *out_\*\*\*\*\*\*.f90*, which should be encapsulated so as to avoid interference and misuse. Reading GKV binary file is done by calling *diag_rb* module in each diagnostics module. Although there are some diagnostics modules implemented, users can design a new diagnostics module by themselves.


### How to use diag
Expanding source codes of **diag**:
- diag/
  - README.md
  - Makefile
  - go.diag (Batch script)
  - backup/ (Makefile and go.diag in some environments)
  - plot_gn/ (Sample plot files for gnuplot)
  - plot_ipynb/ (Sample plot files for python)
  - src
    - diag_header.f90 (Module for setting grid resolutions and MPI processes in GKV)
    - diag_main.f90 (Main program calling each diagnostics module)
    - diag_rb.f90 (Module for reading GKV binary output)
    - diag_\*\*\*\*\*\*.f90 (Module for other settings)
    - ...
    - out_\*\*\*\*\*\*.f90 (Module for each diagnostics)
    - ...

How to use **diag** is in the following steps:
1. Setting parameters in *diag/src/diag_header.f90*
    ```
    %%% DIAG parameters %%%
    integer, parameter :: snum = 1 ! begining of simulation runs
    integer, parameter :: enum = 1 ! end of simulation runs
    !%%%%%%%%%%%%%%%% ! Set run numbers covering diagnosed time range.

    !%%% GKV parameters %%%
    integer, parameter :: nxw = 2, nyw = 8
    integer, parameter :: nx = 0, global ny = 5 ! 2/3 de-aliasing rule
    integer, parameter :: global nz = 64, global nv = 24, global nm = 15
    integer, parameter :: nzb = 3, & ! the number of ghost grids in z
                          nvb = 3 ! the number of ghost grids in v and m

    integer, parameter :: nprocw = 1, nprocz = 4, nprocv = 2, nprocm = 2, nprocs = 2
    !%%%%%%%%%%%%%%%% ! These should be same as gkvp_header.f90.
    ```

2. Calling diagnostics modules in *diag/src/diag_main.f90*
3. Setting the output directory of GKV, *DIR*, in *go.diag*
4. Compile & Execution
5. Output data is dumped in *$DIR/post/*

