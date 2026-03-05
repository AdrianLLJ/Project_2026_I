# Polymer Initialization Program

## Description
Generates an initial 3D configuration of a polymer chain based on user-defined parameters. Outputs an XYZ file for visualization.

## Files
- `constants.f90` – Physical constants.  
- `system.f90` – Polymer system variables.  
- `io_module.f90` – Input/output routines.  
- `init_config.f90` – Allocation and polymer initialization.  
- `main.f90` – Main program.  
- `input.dat` – User input (polymer size, bond length, bond angle).  

## Compilation
Requires `gfortran`. Run in the project directory:

```bash
make
```
## Running
```bash
./polymer
```
## Cleaning
```bash
make clean
```
