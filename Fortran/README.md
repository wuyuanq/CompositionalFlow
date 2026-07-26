# Compositional Two-Phase Flow Simulator

## Overview

This repository contains a **serial Fortran simulator for 2D two-phase
compositional flow**. The simulator supports an arbitrary number of
fluid components (subject to correct input-file configuration) and can
also be used for **single-phase flow simulations**.

------------------------------------------------------------------------

## Features

-   2D two-phase compositional flow simulation
-   Support for an arbitrary number of components
-   Optional **Sparse Grid surrogate** for accelerating PT flash
    calculations
-   Conventional flash calculation mode
-   Single-phase flow simulation support

------------------------------------------------------------------------

## Repository Notes

### Selecting the Flash Method

The flash calculation method is controlled by the `FLASHTYPE` option in
the `Makefile`.

-   `FLASHTYPE = SPARSE` → use Sparse Grid surrogate
-   Otherwise → use conventional flash calculations

### Sparse Grid Setup

Before running with Sparse Grids:

1.  Generate the Sparse Grid surrogate.
2.  Copy the surrogate files into the `Fullgrid` directory.
3.  Supported systems:
    -   Up to **3 components**
    -   Sample surrogates are provided in:
        -   `Fullgrid_2c`
        -   `Fullgrid_3c`

Before running:

-   **2-component case:** copy files from `Fullgrid_2c` to `Fullgrid`
-   **3-component case:** copy files from `Fullgrid_3c` to `Fullgrid`

### TABLESIZE

Set `TABLESIZE` in:

``` text
RST_compositionalTwoPhaseFlow.F90
```

to match the size of the Sparse Grid surrogate.
Set TABLESIZE to match the total number of nodes in the Sparse Grid surrogate. For example:

2-component system: 129 pressure (p) grid points and 129 composition (z₁) grid points
TABLESIZE = 129 × 129
3-component system: 17 pressure (p) grid points, 257 z₁ grid points, and 257 z₂ grid points
TABLESIZE = 17 × (1 + 257) / 2 × 257

### Dummy Flash File

Rename the appropriate dummy flash file before compilation:

-   `RST_dummyFlashcalculation-2c.F90` (2 components)
-   `RST_dummyFlashcalculation-3c.F90` (3 components)

to

``` text
RST_dummyFlashcalculation.F90
```

### Input Files

Before running a case, remove the numbers at the end of the input-file
names.

------------------------------------------------------------------------

## Results

Simulation results are generated in the `case*` directories.

Each case contains:

-   `matlabplot.m`

which generates figures inside the `matlabplots` directory.

------------------------------------------------------------------------

## Citation

If you use this code in academic research, please cite the corresponding
publication.

------------------------------------------------------------------------

## Author

**Yuanqing Wu**\
KAUST, Saudi Arabia

Email: **wuyuanq@gmail.com**
