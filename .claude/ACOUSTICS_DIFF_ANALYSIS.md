# Acoustics Solver: Divergence Analysis — libParanumal vs. libparanumal-dtu

**Date of analysis:** 2026-04-03  
**Original repo:** `libParanumal` (this repo)  
**Fork repo:** `libparanumal-dtu` (`~/github/libparanumal-dtu`)  
**Fork point commit:** `2f1d0a43` — "removed" (2016-12-15)

---

## 1. Overview

Both repositories share identical commit history in `solvers/acoustics/` through the OCCA 1.0 migration in mid-2018. After that point, the fork (`libparanumal-dtu`) began adding research-specific features for room acoustics simulation, while the original (`libParanumal`) continued evolving the framework infrastructure (C++ class refactoring, mesh/library reorganization, v0.5.0 release).

| | libParanumal (original) | libparanumal-dtu (fork) |
|--|--|--|
| Commits to `solvers/acoustics/` | 47 | 73 |
| Language | C++ | C |
| Last commit | 2022-06-07 (v0.5.0) | 2024-01-07 |
| Design pattern | OOP class hierarchy (`solver_t`) | Procedural structs |
| Primary purpose | General-purpose DG solver library | Room acoustics research / ML dataset generation |

---

## 2. Timeline: Changes in the Original (libParanumal)

Changes made in the original repo after the shared history ends:

| Date | Commit | Description |
|------|--------|-------------|
| 2018-07-16 – 2018-07-29 | Multiple | OCCA 1.0 migration: tagged memory, OKL macros removed, mapped pointers, MPI rank/comm changes *(shared with fork)* |
| 2018-08-13 | `b0cfc6c8` | Switch to new `gatherScatter` |
| 2018-11-15 | `edc87899` | Fix acoustics Hex3D for roofline/autoTester |
| 2018-12-07 | `b633802c` | Remove mesh `q` from acoustics |
| 2019-03-22 | `ac48202a` | New settings classes; tested in acoustics |
| 2019-03-22 | `cbcec091` | Move mesh code into its own lib, class structure |
| 2019-04-14 | `ed35abba` | **Prototype refactoring of acoustics solver** |
| 2019-05-13 | `67656695` | Sync with OCCA master; small bugfixes |
| 2019-05-19 | `a09426d4` | Add uniform BOX mesh option with boundary types |
| 2019-06-01 | `3851434f` | Add basic linear algebra launcher |
| 2019-06-01 | `c4fbc53b` | Improved norm reporting (mass matrix norm) |
| 2019-07-06 | `5586dcc8` | Move halo exchange into `ogs` lib |
| 2019-07-20 | `dfc71fb0` | Transpose options for `ogs` lib |
| 2019-07-21 | `29b6ded2` | **Add trace halo exchange to mesh + acoustics solver** |
| 2019-10-11 – 2019-10-28 | Multiple | Update other solvers; fix OCCA kernel cleanup segfaults |
| 2020-01-12 | `089c5bb2` | **Propagate settings from input files to solvers; EXTBDF/BDF with subcycling timesteppers** |
| 2020-07-19 | `6d3ec5a3` | **Clean up cubature and operator data in mesh; uniform operator naming in all kernel code** (`D` → `DT` rename) |
| 2020-07-21 | `e7ee78ef` | CI build and kernel fixes for serial runs |
| 2020-08-05 – 2020-08-11 | Multiple | Larger mesh coverage, IPDG coverage, scripted testing framework |
| 2020-08-15 | `1b03d0a6` | Remove internal BLAS/LAPACK, require user-supplied libs |
| 2020-08-24 | `d5849b63` | New platform class; reorganize library structure |
| 2020-09-11 | `60c6b891` | ParAlmond and libs updates |
| 2020-12-24 | `bfe2b785` | Tensor product interpolation to plotting nodes for Quads/Hexes |
| 2020-12-31 | `3ed0e3fe` | **Add CFL condition calculation to solvers** |
| 2022-06-07 | `cee9c5ee` | **v0.5.0 release** (final commit) |

**Net effect on acoustics:** Framework-level modernization. The solver API was cleaned up (settings classes, OGS halo exchange, platform abstraction), kernel code was made consistent across all solvers (`DT` naming, `@barrier` placement, `@global` pointer annotations), and the CFL calculation was added.

---

## 3. Timeline: Changes in the Fork (libparanumal-dtu)

The fork shares the same OCCA 1.0 migration commits but diverged with research features starting in mid-2019:

| Date | Commit | Description |
|------|--------|-------------|
| 2018-07-16 – 2018-12-07 | Multiple | OCCA 1.0 migration *(shared with original)* |
| 2019-06-06 | `ecbc2b28` | WIP: Remove pointer annotation for OCCA:CUDA compatibility |
| 2019-09-26 | `cc64b13b` | **Change to acoustics equations; add receiver interpolation support (1 receiver)** |
| 2019-11-06 | `f77d392d` | Add frequency-independent (FreqIndep) and locally reacting (LR) BCs |
| 2019-11-07 | `46f45dbc` | **Support for multiple receivers** |
| 2019-11-08 | `77ce4e44` | Fix multiple receivers for multi-CPU/GPU setups |
| 2019-11-13 | `e2dc812a` | **Implicit/Explicit time-stepping for LR BCs** |
| 2019-11-25 | `a494bfa1` | **Extended Reaction (ER) BCs — single CPU/GPU** |
| 2019-11-27 | `5fdfe97a` | Move receiver interpolation onto device (GPU kernel) |
| 2019-12-04 | `1e81c5fe` | Partial fix for multi-CPU/GPU ER BCs |
| 2019-12-05 | `42b8e87a` | WIP: Multiple GPUs for ER BCs |
| 2019-12-12 | `a5dda6db` | Fix ER multi-GPU; add snapshot support during LSERK runs |
| 2020-01-13 | `8f7299a7` | Fix CudaErrorInvalidValue |
| 2020-01-13 | `e61785f3` | Fix locked angle |
| 2020-01-17 | `4269c7a6` | **Add upwind flux — fixing stability issues** |
| 2020-01-17 | `7d3c6c2e` | Add option for ER+LR BCs |
| 2020-02-10 | `33b948f7` | Add setup info; fix crash |
| 2020-04-10 | `073b2e5f` | **Curvilinear mesh support; refactoring; acoustic example** |
| 2020-05-15 | `2672a02d` | Remove data |
| 2020-06-25 | `f6505782` | Create build scripts; move mesh files |
| 2020-07-01 | `b4793499` | **Make code compile/run on macOS** |
| 2020-07-09 | `c1b8cc1d` | Add missing `.dat` files; add simple room simulation example |
| 2020-10-02 | `3dfc3b09` | Add cube mesh example |
| 2020-10-23 | `8c42fe29` | Refactor and add first integration test |
| 2020-11-16 | `8c0c951e` | Integration tests: LR, FreqIndep, perf. refl., combined BCs |
| 2020-11-24 | `a14d98b6` | Add all integration tests (studio room, cylindrical meshes) |
| 2022-02-18 | `982b347a` | Move and update README |
| 2022-09-06 | `3820f54b` | **Implement Gaussian Random Fields (GRFs)** for initial conditions; rename settings fields; bump GCC/CUDA |
| 2022-09-16 | `30c6068b` | Merge `features/grf` |
| 2022-09-20 | `3ec43733` | **Write wave fields in XDMF format** (84 MB vs 2 GB with VTU) |
| 2022-09-21 | `51a61776` | Merge `features/xdmf-writer` |
| 2022-09-27 | `b8d542d2` | Write XDMF in actual Gauss-Lobatto nodes (fix oversampling) |
| 2022-10-07 | `68f2410e` | Generate GRFs on sparse uniform grids; interpolate at quadrature points |
| 2022-10-14 | `0e9ba776` | Option to randomly sample Gaussian source locations |
| 2022-11-15 | `4d4cf3cf` | Make `ppw` for IC mesh a setup parameter |
| 2022-11-30 | `84c9fc8e` | Add GMESH for random source positions; move helpers to utility file |
| 2022-12-04 | `00f39404` | Load additional mesh for source positions from command-line argument |
| 2022-12-18 | `22828a89` | **Implement compact HDF5 writer** (all timesteps in one dataset); remove TXT and VTU formats |
| 2023-04-08 | `3dbe7d50` | Write HDF5 data as `float` (half disk usage); fix H5Compact |
| 2023-04-08 | `e941ce57` | **Write impulse responses as WAV files** (via `tinywav`); add `[DT]` parameter |
| 2024-01-07 | `f08c22f8` | Update unit tests; document CUDA initialization issue workaround |

**Net effect on acoustics:** A full research simulation platform for room acoustics, built on top of the original solver. Key additions are physically realistic boundary conditions (ER, LR), receiver networks, stochastic initial conditions (GRF), and efficient I/O formats for large-scale simulations and ML dataset generation.

---

## 4. File-Level Comparison

### Files only in the original (libParanumal)

| File | Description |
|------|-------------|
| `okl/acousticsInitialCondition2D.okl` | Gaussian IC kernel for 2D elements |
| `okl/acousticsInitialCondition3D.okl` | Gaussian IC kernel for 3D elements |
| `src/acousticsPlotFields.cpp` | VTU/plot output |
| `src/acousticsReport.cpp` | Convergence/norm reporting |
| `src/acousticsSettings.cpp` | Settings class implementation |
| `data/acousticsGaussian2D.h` | 2D IC header |
| `data/acousticsGaussian3D.h` | 3D IC header |
| `setups/setupHex3D.rc` etc. | Standard element-type setup files |

### Files only in the fork (libparanumal-dtu)

| File | Description |
|------|-------------|
| `okl/acousticsERKernel.okl` (~23 KB) | Extended Reaction BC kernel with angle detection and interpolation |
| `okl/acousticsReceiverKernel.okl` | GPU-accelerated receiver interpolation |
| `okl/acousticsUpdate.okl` (~30 KB) | Extended update kernel with multi-stage RK variants |
| `src/acousticsMain.c` | C entry point |
| `src/acousticsSetup.c` (~66 KB) | Extensive simulation setup logic |
| `src/acousticsStep.c` | Time-stepping driver |
| `src/acousticsRun.c` | Simulation driver |
| `src/acousticsReceiver.c` (~28 KB) | Receiver element finding and interpolation |
| `src/acousticsSources.c` | Source management (GRF, random positions) |
| `src/acousticsWriters.c` | I/O: HDF5, XDMF, WAV writers |
| `src/acousticsIO.c` | JSON-based parameter output |
| `src/acousticsError.c` | Error estimation |
| `src/acousticsEstimate.c` | Estimation utilities |
| `src/acousticsUtils.c` | Miscellaneous helpers |
| `src/tinywav.c` | WAV file I/O (third-party) |
| `acousticsTests.h` | Test utilities (file comparison, pressure validation) |
| `acousticsUniform2D.h` | 2D BC macros |
| `acousticsUniform3D.h` | 3D BC macros |
| `tinywav.h` | WAV I/O header (third-party) |
| `tests/` | Catch2 unit and integration test suite |
| `simulationSetups/deeponet/` | ML training dataset generation (130+ config files) |

### Files modified in both repos (shared OKL kernels)

All 8 surface/volume kernels were modified independently in both repos after the fork. Key differences:

| File | Original changes | Fork changes |
|------|-----------------|--------------|
| `acousticsVolumeHex3D.okl` | `D` renamed to `DT` (transpose); added missing `@barrier("local")` | Retains old `D` naming; minor whitespace |
| `acousticsVolumeTet3D.okl` | Same `D`→`DT` rename; `@barrier` added | Retains `DT` name but adds `@barrier`; comment says "isothermal CNS" (legacy label) |
| `acousticsVolumeTri2D.okl` | Updated for uniform operator naming | Similar changes |
| `acousticsVolumeQuad2D.okl` | Updated for uniform operator naming | Similar changes |
| `acousticsSurfaceHex3D.okl` | **Added `@global` qualifier** to `surfaceTerms` helper function arguments (correctness fix for GPU pointer access) | Lacks `@global`; old-style annotations |
| `acousticsSurfaceTet3D.okl` | Standard upwind flux only | **Adds `upwindBC` (with velocity wall term `vn`) and `central` flux functions**; uses physical constants `p_AcConstant`, `p_rho`, `p_c` |
| `acousticsSurfaceTri2D.okl` | Standard upwind flux | Contains extended flux variants |
| `acousticsSurfaceQuad2D.okl` | Standard upwind flux | Contains extended flux variants |

---

## 5. Architectural Differences

| Aspect | Fork (libparanumal-dtu) | Original (libParanumal) |
|--------|------------------------|------------------------|
| **Language** | C (procedural) | C++ (OOP) |
| **Main struct** | `acoustics_t` struct with 100+ fields | `acoustics_t` class inheriting from `solver_t` |
| **Boundary conditions** | Extended Reaction (ER), Locally Reacting (LR), FreqIndep, PerfRefl, combined | Standard: Dirichlet / zero-flux only |
| **Receiver support** | Yes — GPU kernel; multi-receiver, multi-GPU | No |
| **Output formats** | HDF5 (float32), XDMF, WAV | VTU, plot (deprecated in fork) |
| **Initial conditions** | GRF on sparse grids, random source positions | Fixed Gaussian pulse |
| **Testing** | Catch2 unit tests + integration tests | Scripted testing framework (`testAcoustics`) |
| **Curvilinear meshes** | Yes (volume/surface kernel variants) | Not explicitly |
| **Build** | `build_acoustics.sh` + `makefile` | `makefile` only |
| **macOS support** | Added 2020-07-01 | Not specifically tested |
| **Active development** | Yes (last: 2024-01-07) | Stopped at v0.5.0 (2022-06-07) |

---

## 6. OKL Kernel Divergence Summary

The shared OKL kernel files are the most directly mergeable artifacts, but they have diverged in important ways:

### Original has (fork is missing):
- **`@global` pointer annotations** on `surfaceTerms()` helper arguments — this is a correctness fix for GPU memory access in OCCA; the fork lacks this annotation in `acousticsSurfaceHex3D.okl` and likely others.
- **`DT` naming** (transpose of differentiation matrix) — the original renamed `D` → `DT` in volume kernels as part of the v0.19 uniform operator naming (`6d3ec5a3`, 2020-07-19). The fork still uses `D` in some kernels.
- **`@barrier("local")`** — the original added a missing barrier in `acousticsVolumeHex3D.okl`; the fork's corresponding location appears to be missing it.

### Fork has (original is missing):
- **`upwindBC` flux function** — adds a wall velocity term `vn` to the upwind flux, enabling velocity-based boundary conditions (crucial for ER and frequency-dependent BCs).
- **`central` flux function** — a central (non-upwind) flux option for interior faces, not present in the original.
- **Physical constants in flux** — `p_AcConstant`, `p_rho`, `p_c` (density, speed of sound) embedded in surface kernels, reflecting a physically parametric acoustic model.

---

## 7. Features from Fork Worth Cherry-Picking

Listed in rough order of portability to the original C++ framework:

### High priority (physics correctness / bug fixes in shared code):
1. **`@global` annotation fix in surface kernels** — the original already has this; verify the fork's kernels are up to date.
2. **`upwindBC` flux variant** — enables wall velocity BCs; relatively self-contained addition to surface OKL files.
3. **`central` flux function** — adds a flux option for use in mixed BC setups.

### Medium priority (new features with moderate integration effort):
4. **Receiver system** (`acousticsReceiver.c`, `acousticsReceiverKernel.okl`) — GPU-accelerated interpolation at arbitrary receiver positions. Requires adding receiver fields to the setup and adding a kernel call in `acousticsStep`.
5. **Extended Reaction / Locally Reacting BCs** (`acousticsERKernel.okl`) — frequency-dependent boundary conditions. Requires ER/LR coefficient data files and additional fields in the solver struct.
6. **GRF initial conditions** — Gaussian Random Field ICs on sparse grids. Requires a separate IC mesh and interpolation step during setup.

### Lower priority (output/tooling with high fork-specific coupling):
7. **HDF5/XDMF output** — reduces disk usage dramatically for large simulations. Requires linking HDF5.
8. **WAV output** — useful for acoustic impulse responses. Requires `tinywav` (single-file library; easy to add).
9. **macOS build fixes** — if the original is to be built on macOS, the 2020-07-01 fixes from the fork are needed.

---

## 8. Files to Diff Carefully Before Cherry-Picking

These files exist in both repos and contain non-trivial diverged logic:

- `okl/acousticsSurfaceHex3D.okl` — `@global` vs. no annotation
- `okl/acousticsSurfaceTet3D.okl` — standard vs. `upwindBC`/`central`
- `okl/acousticsSurfaceTri2D.okl` — same as Tet3D
- `okl/acousticsSurfaceQuad2D.okl` — same as Tet3D
- `okl/acousticsVolumeHex3D.okl` — `DT` vs `D`, `@barrier` presence
- `okl/acousticsVolumeTet3D.okl` — same
- `src/acousticsSetup.cpp` vs `src/acousticsSetup.c` — completely rewritten; not directly portable but useful as reference for setup logic
- `acoustics.hpp` vs `acoustics.h` — struct fields differ significantly; cherry-picking fields requires integration into the C++ class
