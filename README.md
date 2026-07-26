# CompositionalFlow

**Reference implementations of compositional flow simulation in MATLAB and Fortran, including serial, parallel 2D, and parallel 3D solvers.**

This repository collects four implementations of a compositional-flow simulator:

1. Serial Fortran implementation
2. MPI-based 2D HPC implementation
3. MPI-based 3D HPC implementation
4. MATLAB implementation

The implementations are retained together to support numerical comparison, algorithm development, reproducibility, and performance studies across programming languages and execution models.

> **Repository status:** Research software. The original cluster configurations and compiler versions may require modification on modern systems.

---

## Contents

- [Overview](#overview)
- [Implementations](#implementations)
- [Scientific Scope](#scientific-scope)
- [Repository Organization](#repository-organization)
- [Getting Started](#getting-started)
- [Cross-Implementation Validation](#cross-implementation-validation)
- [HPC Notes](#hpc-notes)
- [Results and Visualization](#results-and-visualization)
- [Roadmap](#roadmap)
- [Citation](#citation)
- [License](#license)
- [Author](#author)

---

## Overview

Compositional-flow simulation models the transport and phase behavior of multicomponent fluids in porous media. This repository provides multiple implementations of the same broader simulation problem so that users can:

- inspect a readable MATLAB prototype;
- run a serial Fortran implementation;
- study distributed-memory parallelization in two dimensions;
- extend the parallel formulation to three dimensions;
- compare numerical results across implementations;
- evaluate computational performance and scalability.

The original project documentation states that the 2D parallel implementation supports two-phase compositional flow, a configurable number of components, single-phase operation, MPI execution, optional Sparse Grid acceleration for flash calculations, and MATLAB-based post-processing.

---

## Implementations

| Implementation | Language / model | Intended use | Directory |
|---|---|---|---|
| Serial solver | Fortran, serial | Baseline implementation and debugging | serial-fortran |
| 2D HPC solver | Fortran + MPI | Parallel two-dimensional simulation | parallel-2d-fortran |
| 3D HPC solver | Fortran + MPI | Parallel three-dimensional simulation | parallel-3d-fortran |
| MATLAB solver | MATLAB | Prototyping, analysis, and visualization | matlab |

### Recommended directory names

For a cleaner public repository, consider renaming the four folders to:

```text
serial-fortran/
parallel-2d-fortran/
parallel-3d-fortran/
matlab/
```

Renaming is optional. Avoid changing directory names until build scripts and internal paths have been checked.

---

## Scientific Scope

Based on the available project documentation, the repository includes or is intended to include:

- compositional flow in porous media;
- two-phase flow;
- configurable fluid components;
- serial and distributed-memory execution;
- two-dimensional and three-dimensional implementations;
- flash calculations;
- optional Sparse Grid surrogate acceleration;
- MATLAB post-processing.

The current documentation does not fully specify the governing equations, discretization, nonlinear solution method, linear solver, well model, boundary conditions, or thermodynamic model. These details should be documented from the source code or associated publications rather than inferred.

---

## Repository Organization

A recommended long-term layout is:

```text
CompositionalFlow/
├── README.md
├── LICENSE
├── CITATION.cff
├── CONTRIBUTING.md
├── .gitignore
│
├── serial-fortran/
│   ├── README.md
│   ├── src/
│   ├── examples/
│   └── Makefile
│
├── parallel-2d-fortran/
│   ├── README.md
│   ├── src/
│   ├── examples/
│   ├── scripts/
│   └── Makefile
│
├── parallel-3d-fortran/
│   ├── README.md
│   ├── src/
│   ├── examples/
│   ├── scripts/
│   └── Makefile
│
├── matlab/
│   ├── README.md
│   ├── src/
│   ├── examples/
│   └── plotting/
│
├── docs/
│   ├── numerical-method.md
│   ├── input-format.md
│   ├── validation.md
│   ├── performance.md
│   └── images/
│
└── tests/
    ├── reference-data/
    └── regression/
```

This target structure should be introduced gradually. A documentation-only improvement can be merged first without moving source files.

---

## Getting Started

Because the four implementations may use different build and execution procedures, each subdirectory should contain its own `README.md` with:

- prerequisites;
- compiler or MATLAB version;
- build command;
- executable or main script;
- input file selection;
- output location;
- one minimal example;
- expected result.

### Serial Fortran

Typical workflow:

```bash
cd REPLACE_WITH_SERIAL_FOLDER
make
./REPLACE_WITH_EXECUTABLE
```

### Parallel 2D Fortran

The MPI process count must be consistent between the input configuration, Makefile or launch command, and cluster script.

```bash
cd Hpc
make
mkdir -p case1
mpirun -np 4 ./CompositionalFlow_fortran_hpc
```

### Parallel 3D Fortran

Typical workflow:

```bash
cd REPLACE_WITH_3D_HPC_FOLDER
make
mkdir -p case1
mpirun -np 8 ./REPLACE_WITH_EXECUTABLE
```

The exact executable name and process-grid configuration should be copied from the existing build files.

### MATLAB

Typical workflow:

```matlab
cd('REPLACE_WITH_MATLAB_FOLDER')
% Run the documented main script.
```

Identify the main MATLAB script in the subdirectory README.

---

## Cross-Implementation Validation

The strongest way to present this repository is to demonstrate that the implementations solve the same reference problems.

Create at least one small validation case with identical physical inputs:

```text
tests/reference-case/
├── README.md
├── input/
├── serial-reference/
├── parallel-2d-reference/
├── parallel-3d-reference/
└── matlab-reference/
```

Report comparisons such as:

| Quantity | Serial Fortran | 2D MPI | 3D MPI | MATLAB | Tolerance |
|---|---:|---:|---:|---:|---:|
| Final pressure | TBD | TBD | TBD | TBD | TBD |
| Phase saturation | TBD | TBD | TBD | TBD | TBD |
| Component balance error | TBD | TBD | TBD | TBD | TBD |

Do not publish numerical values until they have been reproduced from the code.

---

## HPC Notes

The original 2D documentation references historical Shaheen and Neser environments.

### Shaheen

```bash
module load bluegene
```

### Neser

```bash
module unload openmpi/1.5.4/gcc
module load intel-compilers/11.1
module load openmpi/1.6.4/intel
```

These module names and versions are historical. Modern clusters will likely require different compiler, MPI, scheduler, and launch configurations.

For larger grids, memory pressure may be reduced by allocating more nodes and using fewer MPI processes per node.

### Sparse Grid flash calculations

The documented Makefile uses `FLASHTYPE` to select the flash-calculation mode:

```makefile
FLASHTYPE = SPARSE
```

When Sparse Grids are enabled:

1. generate the surrogate files;
2. place them in the documented `Fullgrid/` directory;
3. configure `TABLESIZE` in `RST_compositionalTwoPhaseFlow.F90`;
4. select the correct dummy flash-calculation source for the component count.

The documented Sparse Grid implementation supports up to three components.

---

## Results and Visualization

The documented HPC workflow writes results to directories named similar to:

```text
case1/
case2/
```

A generated `matlabplot.m` script may be used to create MATLAB figures.

For GitHub presentation, add two or three verified images under:

```text
docs/images/
```

Recommended figures:

- pressure field;
- phase saturation;
- component concentration;
- 2D domain decomposition;
- 3D simulation result;
- MPI strong-scaling curve.

Then display them here:

```markdown
![Pressure field](docs/images/pressure-field.png)
```

---

## Performance

Add reproducible performance data for the parallel implementations.

Recommended metadata:

- processor model;
- compiler and version;
- MPI implementation and version;
- grid size;
- component count;
- process topology;
- number of time steps;
- wall-clock time.

Example table template:

| Grid | MPI ranks | Runtime | Speedup | Parallel efficiency |
|---|---:|---:|---:|---:|
| TBD | 1 | TBD | 1.00 | 100% |
| TBD | 4 | TBD | TBD | TBD |
| TBD | 16 | TBD | TBD | TBD |

---

## Roadmap

High-value improvements, in priority order:

1. Document the exact four directory names and entry points.
2. Add one reproducible example to each implementation.
3. Add reference outputs and cross-implementation checks.
4. Document governing equations and numerical methods.
5. Modernize the Fortran build configuration.
6. Add CI compilation for the serial solver.
7. Add MPI smoke tests where supported.
8. Add verified simulation images and scaling results.
9. Archive historical cluster scripts under `scripts/legacy/`.
10. Tag a documented release such as `v1.0.0`.

---

## Citation

A provisional software citation is supplied in [`CITATION.cff`](CITATION.cff). Update it with the associated paper, DOI, and preferred release version when available.

BibTeX template:

```bibtex
@software{wu_compositional_flow,
  author = {Wu, Yuanqing},
  title  = {CompositionalFlow},
  year   = {2015},
  url    = {https://github.com/wuyuanq/CompositionalFlow}
}
```

---

## License

See [`LICENSE`](LICENSE).

Before applying an open-source license, confirm that you have permission to license code developed at KAUST or under a sponsored research project.

---

## Author

**Yuanqing Wu**  
KAUST  
Email: wuyuanq@gmail.com
