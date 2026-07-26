# Contributing

Contributions that improve reproducibility, documentation, portability, validation, and numerical correctness are welcome.

## Before opening an issue

Include:

- implementation used: serial Fortran, 2D MPI, 3D MPI, or MATLAB;
- operating system or HPC platform;
- compiler/MATLAB version;
- MPI implementation and version, when applicable;
- input case;
- process count and process topology;
- complete error output;
- minimal reproduction steps.

## Pull requests

1. Create a focused branch from `main`.
2. Avoid mixing source reformatting with functional changes.
3. Document numerical assumptions, physical units, and array dimensions.
4. Add or update a reproducible example when behavior changes.
5. Compare output with an existing reference case.
6. Do not commit generated binaries, module files, or large result directories.

## Fortran guidance

- Use explicit `implicit none` where possible.
- Keep module responsibilities clear.
- Document procedure inputs, outputs, and units.
- Check MPI return codes where practical.
- Avoid machine-specific absolute paths.

## MATLAB guidance

- Identify the main script or function clearly.
- Avoid relying on the current working directory when possible.
- Keep plotting separate from numerical computation.
- Document required toolboxes.
