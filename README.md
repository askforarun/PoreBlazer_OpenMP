# Overview

PoreBlazer v4.0 is a Fortran tool for geometric characterization of porous materials. It computes pore limiting diameter, largest cavity diameter, helium pore volume, geometric pore volume, probe-occupiable free volume, pore size distributions, and surface area metrics from crystal structure inputs.

This repository keeps the original PoreBlazer workflow while adding performance-focused improvements:

- OpenMP parallelization in core routines
- linked-cell / linked-list acceleration for large orthorhombic systems
- a PSD/free-volume-only mode that skips surface-area calculations
- validation scripts for comparing against the original Sarkisov code

Original PoreBlazer v4.0 code by Lev Sarkisov: https://github.com/SarkisovGitHub/PoreBlazer . This repository is an optimized fork and is maintained separately from the original upstream code.

This repository is distributed under the GNU General Public License v3.0 or later (`GPL-3.0-or-later`). See [LICENSE](./LICENSE).

## Key Features

- Original-code compatibility: the optimized serial code has been validated against the original `SarkisovGitHub/PoreBlazer` case studies.
- OpenMP acceleration: `lattice_calculations`, `surface_area`, and `volumes` can use multiple CPU cores.
- Linked-cell neighbor search: large orthorhombic systems use per-cell linked lists of atoms to reduce atom-lattice distance work.
- Percolation-axis reporting: the limiting-diameter analysis reports whether percolation occurs along `x`, `y`, and/or `z`, and writes this as `D_axes` in `summary.dat`.
- Volume/PSD-only mode: an optional input flag skips both surface-area phases while preserving downstream RNG state for PSD and free-volume calculations.
- Bundled examples: the repository includes reference case studies and a Windows example.

## What This Tool Does

Typical PoreBlazer workflow:

1. Read a structure file and simulation-cell parameters from `input.dat`.
2. Read force-field, probe, lattice, PSD, and RNG settings from `defaults.dat`.
3. Build a geometric lattice and classify cubelets accessible to a point probe, helium, and nitrogen.
4. Compute structural descriptors such as:
   - pore limiting diameter (`PLD`)
   - largest cavity diameter (`LCD`)
   - surface area
   - helium pore volume
   - geometric pore volume
   - probe-occupiable free volume
   - total and network-accessible pore size distributions

## Table of Contents

- Citation
- Repository Contents
- What You Need
- Installation
  - GNU Fortran
  - Intel Fortran
- Running PoreBlazer
  - Basic run
  - OpenMP run
  - PSD/free-volume-only run
- Input Files
  - `input.dat`
  - `defaults.dat`
- Main Outputs
- Performance and Large Systems
- Validation
- Case Studies
- Windows Usage
- License
- Contact

## Citation

If you use PB v4.0, please cite the associated PoreBlazer v4.0 publication(s) documented in `src/README_PB_v4.0.txt`.

## Repository Contents

- `src/`: Fortran source code, Makefiles, executables, and documentation
- `case_studies/`: example structures and reference outputs
- `data/`: comparative porous-material datasets
- `Windows/`: Windows executable and example run setup
- `scripts/`: validation and benchmarking job scripts
- `PB4.0_vs_Zeo++_vs_RASPA.zip`: comparative case-study package
- `case_studies.zip`: extended case-study package

## What You Need

- A Fortran compiler:
  - `gfortran` for GNU builds
  - `ifort` or another Intel-compatible compiler for Intel builds
- Input files in your run directory:
  - `input.dat`
  - `defaults.dat`
  - `UFF.atoms`
  - the target structure file, for example `HKUST1.xyz`
- For OpenMP runs:
  - a compiler build with OpenMP enabled
  - an appropriate `OMP_NUM_THREADS` setting

## Installation

Build from the repository root:

```bash
cd src
make
```

This generates:

- `poreblazer.exe`

### GNU Fortran

```bash
cd src
make -f Makefile_gfort clean
make -f Makefile_gfort
```

Uses `-fopenmp`.

### Intel Fortran

```bash
cd src
make -f Makefile_intel clean
make -f Makefile_intel
```

Uses `-qopenmp`.

You can also override the compiler directly:

```bash
cd src
make clean
make FORTRAN_COMPILER=gfortran
```

or

```bash
cd src
make clean
make FORTRAN_COMPILER=ifort
```

## Running PoreBlazer

### Basic run

```bash
export OMP_NUM_THREADS=1
./poreblazer.exe < input.dat
```

### OpenMP run

```bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=close
export OMP_PLACES=cores
./poreblazer.exe < input.dat
```

### PSD/free-volume-only run

Use the same `input.dat` as a normal run, but set `nsample` to `0` in `defaults.dat`:

```bash
export OMP_NUM_THREADS=1
./poreblazer.exe < input.dat
```

This skips both surface-area phases while still running the lattice, PSD, and volume calculations.

### Example case-study run

From the repository root:

```bash
cd case_studies/HKUST1
export OMP_NUM_THREADS=1
../../src/poreblazer.exe < input.dat
```

### Example `sbatch` run

For a threaded case-study run on a cluster:

```bash
sbatch --job-name=pb-hkust1 \
  --time=02:00:00 \
  --cpus-per-task=8 \
  --chdir=/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/case_studies/HKUST1 \
  --output=/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/case_studies/HKUST1/slurm-%j.out \
  --error=/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/case_studies/HKUST1/slurm-%j.err \
  --wrap="bash -lc 'export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores; /users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/src/poreblazer.exe < input.dat'"
```

### Example polymer `sbatch` runs

Upstream serial:

```bash
sbatch --job-name=pb-polymer-upstream \
  --time=24:00:00 \
  --cpus-per-task=1 \
  --chdir=/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/polymers/polymer_upstream \
  --output=/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/polymers/polymer_upstream/slurm-%j.out \
  --error=/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/polymers/polymer_upstream/slurm-%j.err \
  --wrap="bash -lc 'export OMP_NUM_THREADS=1; /usr/bin/time -p -o time.txt ./poreblazer_upstream.exe < input.dat > run.out'"
```

Optimized serial with linked cells:

```bash
sbatch --job-name=pb-polymer-serial \
  --time=24:00:00 \
  --cpus-per-task=1 \
  --chdir=/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/polymers/polymer_linkedcell \
  --output=/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/polymers/polymer_linkedcell/slurm-%j.out \
  --error=/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/polymers/polymer_linkedcell/slurm-%j.err \
  --wrap="bash -lc 'export OMP_NUM_THREADS=1; /usr/bin/time -p -o time.txt /users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/src/poreblazer.exe < input.dat > run.out'"
```

Optimized OpenMP with linked cells:

```bash
sbatch --job-name=pb-polymer-openmp \
  --time=24:00:00 \
  --cpus-per-task=8 \
  --chdir=/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/polymers/polymer_openmp \
  --output=/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/polymers/polymer_openmp/slurm-%j.out \
  --error=/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/polymers/polymer_openmp/slurm-%j.err \
  --wrap="bash -lc 'export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores; /usr/bin/time -p -o time.txt /users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/src/poreblazer.exe < input.dat > run.out'"
```

## Input Files

### `input.dat`

Example `input.dat` for `HKUST1`:

```text
HKUST1.xyz
26.28791 26.28791 26.28791
90 90 90
```

This format is the same for normal runs and PSD/free-volume-only runs.

### `defaults.dat`

Example `defaults.dat` for a normal run:

```text
UFF.atoms
2.58, 10.22, 298, 12.8
3.314
500
0.2
20.0, 0.25
21908391
0
```

Annotated `defaults.dat`:

```text
UFF.atoms              ! atom-type force-field file
2.58, 10.22, 298, 12.8 ! He sigma (A), He epsilon (K), temperature (K), LJ cutoff (A)
3.314                  ! N2 sigma (A)
500                    ! nsample: surface-area samples per atom; use 0 to skip surface area
0.2                    ! lattice cube size (A)
20.0, 0.25             ! largest anticipated pore diameter (A), PSD bin size (A)
21908391               ! random seed
0                      ! visualization option: 0 none, 1 xyz, 2 grd, 3 both
```

Example `defaults.dat` for a PSD/free-volume-only run with surface area skipped:

```text
UFF.atoms
2.58, 10.22, 298, 12.8
3.314
0
0.2
20.0, 0.25
21908391
0
```

Safe starting example for new systems:

```text
UFF.atoms
2.58, 10.22, 298, 12.8
3.314
100
0.25
20.0, 0.25
21908391
0
```

This is a good conservative starting point when you want a reasonable balance between runtime and output quality before tightening the grid or increasing `nsample`.

Polymer example for a large orthorhombic system:

```text
UFF.atoms
2.58, 10.22, 298, 12.8
3.314
20
1.0
20.0, 1.0
21908391
0
```

This kind of coarse setup is useful for large polymer-like systems where you want an initial validation run before moving to finer `cube_size`, smaller PSD bins, or larger `nsample`.

`defaults.dat` controls:

- probe parameters
- cutoff distance
- surface-area sampling count
- lattice cube size
- PSD binning settings
- random seed
- visualization options

Notes:

- If `nsample = 0`, surface-area calculations are skipped while preserving the downstream RNG state used by PSD and free-volume calculations.
- If `nsample = 0`, the `S_AC_*` entries in `summary.dat` are written as `0.00`, and the run becomes a volume/PSD-focused workflow.
- With `nsample = 0`, expect `summary.dat` to contain `S_AC_A^2 0.00`, `S_AC_m^2/cm^3 0.00`, and `S_AC_m^2/g 0.00` while the PSD and volume outputs are still produced normally.
- If `nsample < 0`, the code stops with an input error.
- The optimized fork uses the original 8-line `defaults.dat` format; there is no extra `surface_area_option` line.
- Very small `cube_size` and very large `nsample` can be expensive on large systems.
- If the reported `LCD` is close to the configured largest anticipated pore diameter, increase the line-6 upper bound and rerun the PSD analysis.

### Complete no-surface-area example

`input.dat`

```text
HKUST1.xyz
26.28791 26.28791 26.28791
90 90 90
```

`defaults.dat`

```text
UFF.atoms
2.58, 10.22, 298, 12.8
3.314
0
0.2
20.0, 0.25
21908391
0
```

## Main Outputs

Key results are written to `summary.dat`, including:

- density
- pore limiting diameter (`PLD`)
- largest cavity diameter (`LCD`)
- total and network-accessible surface area
- helium pore volume
- geometric pore volume
- probe-occupiable free volume

Common output files:

- `Total_psd_cumulative.txt`
- `Total_psd.txt`
- `Network-accessible_psd_cumulative.txt`
- `Network-accessible_psd.txt`
- `probe_occupiable_volume.xyz`
- `nitrogen_network.xyz` if enabled
- `nitrogen_network.grd` if enabled

## Performance and Large Systems

PoreBlazer v4.0 is designed to handle systems ranging from small unit cells to very large structures with `100,000+` atoms.

### Linked Cell Optimization

For large systems, `lattice_calculations` uses a linked-cell / cell-list algorithm backed by per-cell linked lists of atoms.

- Activation: automatically enabled for orthorhombic cells with `2,000+` atoms
- Benefit: reduces atom-lattice distance work by traversing nearby cell-linked atom lists for overlap checks, helium Lennard-Jones accumulation, and nearest-atom distances
- Note: the `2,000`-atom threshold is a current built-in heuristic used by the validated implementation for the large-system linked-cell path.

### OpenMP

OpenMP parallelization is enabled in:

- `lattice_calculations`
- `surface_area`
- `volumes`

Build details and optimization notes are in `src/README_OPENMP.txt`.

Recommended runtime settings:

- use physical cores rather than hyperthreads when possible
- set `OMP_NUM_THREADS` to a value appropriate for your machine
- optionally set:

```bash
export OMP_PROC_BIND=close
export OMP_PLACES=cores
```

### No-surface-area mode

If you only need PSD and free-volume results, skipping surface area can reduce total runtime by removing the Monte Carlo surface-sampling phases.

## Validation

The repository includes validation scripts under `scripts/`.

Available validation workflows:

- `scripts/compare_case_studies.sbatch`
  compares the optimized code against the original Sarkisov code in serial mode
- `scripts/openmp_validation.sbatch`
  compares original serial, optimized serial, and optimized OpenMP runs on the bundled case studies

### MOFs, Zeolites, and Slit Pores

Validation summary:

- the optimized serial code matched the original `SarkisovGitHub/PoreBlazer` outputs exactly across all 15 bundled case studies
- the optimized OpenMP code also matched the original and optimized serial outputs exactly across all 15 bundled case studies

Validation and timing table:

| Case | Atoms | Upstream serial (s) | Optimized serial (s) | Optimized OpenMP 8 threads (s) | Exact match |
|---|---:|---:|---:|---:|---|
| BHP | 84 | 6.74 | 6.73 | 6.31 | Yes |
| CD121 | 208 | 16.11 | 16.15 | 13.59 | Yes |
| CLO | 5136 | 129.99 | 58.73 | 42.05 | Yes |
| HKUST1 | 624 | 74.63 | 73.98 | 56.97 | Yes |
| IRMOF1 | 424 | 57.02 | 57.01 | 45.49 | Yes |
| LOV | 972 | 2.63 | 2.62 | 0.83 | Yes |
| MIL101 | 11768 | 1513.96 | 385.76 | 321.11 | Yes |
| MIL47V | 1152 | 96.37 | 95.41 | 56.44 | Yes |
| MOF180 | 926 | 555.51 | 556.27 | 315.85 | Yes |
| ROG | 240 | 3.34 | 3.30 | 2.56 | Yes |
| RON | 240 | 3.44 | 3.43 | 2.65 | Yes |
| SLIT | 882 | 1.74 | 1.81 | 0.57 | Yes |
| STO | 100 | 5.61 | 5.60 | 4.09 | Yes |
| WEI | 68 | 0.26 | 0.23 | 0.07 | Yes |
| ZIF8 | 2208 | 206.37 | 184.12 | 83.10 | Yes |

These results were generated with the original upstream-style case-study inputs and `defaults.dat` files, with surface area calculation enabled throughout. The validation artifacts are stored under [validation/319416_20260426_213105](/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/validation/319416_20260426_213105), with the summary table in [openmp_validation.csv](/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/validation/319416_20260426_213105/openmp_validation.csv). Re-run `scripts/openmp_validation.sbatch` after code changes to regenerate the table.

### Polymers

The repository also includes a separate non-bundled polymer validation, added in this optimized repository for a large polymer-like system with `79,650` atoms:

- [polymers/polymer_upstream](/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/polymers/polymer_upstream)
- [polymers/polymer_linkedcell](/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/polymers/polymer_linkedcell)
- [polymers/polymer_openmp](/users/ass2009/sharedscratch/PoreBlazer-performance-optimizations/polymers/polymer_openmp)

For this polymer case, the optimized serial and optimized OpenMP runs match exactly, and the upstream run matches all meaningful physical outputs with surface area calculation enabled:

| Quantity | Atoms | Upstream serial | Optimized serial | Optimized OpenMP |
|---|---:|---:|---:|---:|
| `PLD_A` | 79650 | 0.90 | 0.90 | 0.90 |
| `LCD_A` | 79650 | 10.12 | 10.12 | 10.12 |
| `S_AC_A^2` | 79650 | 2599.17 | 2599.17 | 2599.17 |
| `S_AC_m^2/cm^3` | 79650 | 36.06 | 36.06 | 36.06 |
| `S_AC_m^2/g` | 79650 | 33.53 | 33.53 | 33.53 |
| `V_He_A^3` | 79650 | 49922.619 | 49922.619 | 49922.619 |
| `V_G_A^3` | 79650 | 221082.309 | 221082.309 | 221082.309 |
| `V_PO_A^3` | 79650 | 14790.406 | 14790.406 | 14790.406 |
| `FV_PO` | 79650 | 0.02052 | 0.02052 | 0.02052 |

Polymer timing and performance:

| Configuration | Wall Time (s) | Speedup vs Upstream |
|---|---:|---|
| Upstream serial | 523.13 | 1.0× |
| Optimized serial (linked-cell, 1 thread) | 21.94 | 23.8× |
| Optimized OpenMP (8 threads) | 5.36 | 97.4× |

The corresponding `summary.dat`, `Total_psd.txt`, and `Total_psd_cumulative.txt` files are stored in each of those `polymers/` subdirectories. For this non-percolating polymer case, the original upstream `D` field is not used as a validation target; the optimized code reports `D = 0` and `D_axes none`.

## Windows Usage

Windows users can run the precompiled binary in `Windows/`.

A ready `HKUST1` example is included in `Windows/HKUST1`; run `run.bat` there to execute and write output to `results.txt`.

## License

GNU General Public License v3.0 or later.

## Contact

Arun Srikanth Sridhar  
askforarun@gmail.com
