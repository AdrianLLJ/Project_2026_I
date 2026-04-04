---
Author: J-dot-Barrientos
Date: 2026-04-01
---

# Directory parallel

Contains the same module subdirectories as `sequential`: `initialization`, `energy` and `MC_update`.
Only the first one has changes from the sequential version. io_module contains a new subroutine that
allows the broadcasting of the variables of the input file.

`initialization` also contains a test program for initialization only, `init_replicas.f90`, and the
input file for the same purpose.

Files added in this directory: `mainglobal.f90`, which now uses mpi and allows to make several simulations
as independent replicas starting from different initial configurations (thus, with initRandom input
variable with value 1) using different seeds per replica. `parallel.mk`, for compilation, adapted from
 ../sequential/sequential.mk. `input.dat` which contains the simulation parameters.

## Compilation

```bash
make -f parallel.mk
```

## Run

```bash
mpirun -np <number.of.processors> mainglobal_mpi.exe results < input.dat
```

## Clean

```bash
make -f parallel.mk clean
```

or

```bash
rm -rf build
```

**Note**: a results subdirectory will be created if it does not exist.
It will contain subdirectories named `replica_<number>` with the results
of each parallel simulation. One of the next steps in the parallel version
is the gathering of the output data.
