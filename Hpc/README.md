
# Parallel Compositional Flow Simulator

<div align="center">

**A High-Performance MPI-Based Reservoir Simulator for 2D Compositional Multiphase Flow**

Modern Fortran • MPI • Reservoir Simulation • High Performance Computing

</div>

---

## Overview

This repository contains a parallel **Fortran** simulator for **2D two-phase compositional flow in porous media**. The code is designed for distributed-memory high-performance computing (HPC) environments using MPI and supports an arbitrary number of fluid components through configurable input files.

In addition to compositional two-phase flow, the framework can also be configured for single-phase simulations. Optional Sparse Grid surrogate models are available to accelerate flash calculations.

---

# Features

- Parallel reservoir simulation using MPI
- 2D compositional two-phase flow
- Configurable number of fluid components
- Single-phase simulation support
- Sparse Grid surrogate for flash calculation
- Designed for HPC clusters
- MATLAB-based post-processing

---

# Requirements

- Modern Fortran compiler
- MPI implementation (OpenMPI / Intel MPI)
- Mac
- Make
- MATLAB (optional)

---

# Build

```bash
make
```

---

# Quick Start

Compile the project

```bash
make -f Makefile.mac
```

Run with MPI

```bash
mpirun -np 4 ./CompositionalFlow_fortran_hpc
```

The number of MPI processes must remain consistent in:

1. The input file (`proceAlloc(...)`)
2. The Makefile
3. The cluster job script

---

# Sparse Grid Flash Calculation

The Makefile variable

```makefile
FLASHTYPE
```

controls the flash solver.

- `FLASHTYPE = SPARSE` → Sparse Grid surrogate
- Otherwise → Direct flash calculation

Before using Sparse Grids:

1. Generate surrogate tables.
2. Copy them into the `Fullgrid` directory.
3. Configure `TABLESIZE` correctly.
4. Rename the corresponding dummy flash source according to the number of components.

---

# Simulation Workflow

```text
Input Files
      │
      ▼
Grid Initialization
      │
      ▼
Pressure Solver
      │
      ▼
Flash Calculation
      │
      ▼
Transport Solver
      │
      ▼
Output Files
      │
      ▼
MATLAB Visualization
```

---

# Numerical Methods

According to the original documentation, the simulator supports:

- Parallel MPI execution
- Compositional multiphase flow
- Flash calculations
- Sparse Grid surrogate acceleration

Further numerical implementation details should be documented from the corresponding publications or source code.

---

# Example Output

Simulation results are written into directories named

```text
case*
```

Each directory contains

```text
matlabplot.m
```

which generates visualization figures.

---

# Citation

If you use this project in academic research, please cite the associated publication (if available).

---

# Contact

**Author**

Yuanqing Wu

KAUST

Email: wuyuanq@gmail.com

---

# License

This project is licensed under the MIT License.
