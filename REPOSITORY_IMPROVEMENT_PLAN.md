# Repository Improvement Plan

This document separates low-risk documentation improvements from source-code restructuring.

## Phase 1 — Safe documentation upgrade

Add without moving existing source files:

- root `README.md`;
- `.gitignore`;
- `CITATION.cff`;
- `CONTRIBUTING.md`;
- one README inside each implementation folder;
- `docs/images/` for verified figures.

## Phase 2 — Reproducible examples

For each implementation, identify:

- the build command;
- executable or main MATLAB script;
- smallest working input;
- expected runtime;
- output files;
- expected numerical result.

## Phase 3 — Validation

Create common reference cases and compare:

- pressure;
- saturation;
- component quantities;
- conservation error;
- serial-versus-parallel differences.

## Phase 4 — Directory cleanup

Only after builds are reproducible:

- move source files into `src/`;
- move input cases into `examples/`;
- archive cluster-specific scripts under `scripts/legacy/`;
- place generated figures outside source directories;
- update paths in Makefiles and scripts.

## Phase 5 — Automation

- compile serial Fortran in GitHub Actions;
- run a minimal regression test;
- optionally compile MPI implementations;
- publish documentation and release artifacts.

## Phase 6 — Public release

- confirm software ownership;
- choose a license;
- add publication metadata;
- tag `v1.0.0`;
- upload a Zenodo archive if appropriate.
