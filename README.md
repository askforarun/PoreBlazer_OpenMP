# PoreBlazer v4.0

PoreBlazer (v4.0) source code, examples, and reference geometric properties of porous materials.

**Code Refactoring**: The core codebase has been refactored for better maintainability and performance.

## Recent changes

The latest update adds OpenMP-based parallel performance optimizations in core routines:

- `lattice_calculations` (major runtime contributor)
- `surface_area` (Monte Carlo surface area sampling)
- `volumes` (helium pore-volume calculations)

Build scripts were also updated to compile with OpenMP by default:

- `src/Makefile`: GNU default and OpenMP flags enabled
- `src/Makefile_gfort`: `-fopenmp`
- `src/Makefile_intel`: `-qopenmp`

Detailed optimization and benchmarking notes are in:

- `src/README_OPENMP.txt`

## Repository contents

- `src/`: Fortran source code, Makefiles, prebuilt executables, and docs
- `data/`: MOF subset property datasets (`MOFsubsetPB4.dat`, `MOFsubsetZeo++.dat`, `MOFsubsetRASPA.dat`)
- `case_studies/`: Example structures and reference outputs
- `Windows/`: Windows executable and HKUST-1 example run setup
- `PB4.0_vs_Zeo++_vs_RASPA.zip`: Comparative case-study package
- `case_studies.zip`: Extended case-study package

## Build

From repository root:

```bash
cd src
make
```

This generates:

- `poreblazer.exe`

## Compilation instructions

From the `src/` directory, use one of the following:

### 1. Default build (uses `src/Makefile`)

```bash
cd src
make clean
make
```

This follows `FORTRAN_COMPILER` in `src/Makefile` (currently `gfortran`) and enables OpenMP.

### 2. GNU Fortran build

```bash
cd src
make -f Makefile_gfort clean
make -f Makefile_gfort
```

Uses `-fopenmp`.

### 3. Intel Fortran build

```bash
cd src
make -f Makefile_intel clean
make -f Makefile_intel
```

Uses `-qopenmp`.

### 4. Override compiler without editing files

```bash
cd src
make clean
make FORTRAN_COMPILER=ifort
```

Or:

```bash
cd src
make clean
make FORTRAN_COMPILER=gfortran
```

## Run

From your run directory, make sure these files are present:

- `input.dat`
- `defaults.dat`
- `UFF.atoms`
- target structure file (for example `HKUST1.xyz`)

Run:

```bash
./poreblazer.exe < input.dat
```

## OpenMP usage

Control thread count with `OMP_NUM_THREADS`:

```bash
export OMP_NUM_THREADS=16
./poreblazer.exe < input.dat
```

Optional affinity settings for multi-core systems:

```bash
export OMP_PROC_BIND=close
export OMP_PLACES=cores
```

For validation, compare 1-thread and multi-thread runs; results should match (or differ only at tiny floating-point roundoff level).

## Input modes

### Basic mode

Specify structure file and unit-cell parameters in `input.dat`.

Example:

```text
HKUST1.xyz
26.28791 26.28791 26.28791
90 90 90
```

### Advanced mode

Edit `defaults.dat` to control:

- probe parameters
- cutoff distance
- lattice grid size
- PSD binning settings
- random seed
- optional network-visualization output

## Main outputs

Key summary appears in `summary.dat`, including:

- density
- pore limiting diameter (PLD)
- largest cavity diameter (LCD)
- total and network-accessible surface area
- total and network-accessible pore volumes

Common generated files:

- `Total_psd_cumulative.txt`
- `Total_psd.txt`
- `Network-accessible_psd_cumulative.txt`
- `Network-accessible_psd.txt`
- `probe_occupiable_volume.xyz`
- `nitrogen_network.xyz` (if enabled)
- `nitrogen_network.grd` (if enabled)

## Case studies

Included examples cover MOFs, zeolites, and slit pore systems, with reference outputs for comparison.

## Windows usage

Windows users can run the precompiled binary in `Windows/`.

A ready HKUST-1 example is provided in `Windows/HKUST1`; run `run.bat` there to execute and write output to `results.txt`.

## Performance and Large Systems

PoreBlazer (v4.0) is designed to handle systems ranging from small unit cells to large structures with **100,000+ atoms**.

### Linked Cell Optimization
For large systems, PoreBlazer implements a **linked cell (cell list)** algorithm. 
- **Activation**: Automatically enabled as an overlap prefilter for orthorhombic cells with 2,000 or more atoms.
- **Benefit**: Significantly reduces the computational cost of detecting atom-lattice overlaps.

### Handling 100,000+ Atoms
- **Technical Capacity**: Uses 4-byte integers for atom counts, supporting up to ~2 billion atoms.
- **Memory Efficiency**: Memory usage for atom storage is optimized (e.g., 100,000 atoms require ~10 MB).
- **Complexity**: The core distance calculation is $O(N_{atoms} \times N_{grid})$. For very large systems, the performance is best maintained using:
  - **Orthorhombic cells** (to take advantage of linked cell pre-filtering).
  - **OpenMP Parallelization** (to distribute the workload across CPU cores).

## Citation

If you use PB v4.0, please cite the associated PoreBlazer v4.0 publication(s) documented in `src/README_PB_v4.0.txt`.

## License

GNU General Public License v3.0 or later.

## Contact

Lev Sarkisov  
lev.sarkisov@manchester.ac.uk
