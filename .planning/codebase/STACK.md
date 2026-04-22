# Technology Stack

**Analysis Date:** 2026-04-22

## Languages

**Primary:**
- C++17 - All solver logic, framework libraries, and host-side code. Compiler flag `-std=c++17` set in `make.top`.
- OKL (OCCA Kernel Language) - GPU/accelerator compute kernels. `.okl` files in `solvers/acoustics/okl/`. OKL is a thin annotation layer on top of C: `@kernel`, `@outer`, `@inner`, `@shared`, `@exclusive`, `@restrict` directives.

**Secondary:**
- C (headers) - Legacy `.h` include files for boundary/initial condition macros, e.g., `solvers/acoustics/data/acousticsGaussian2D.h`.
- Python 3 - Integration test harness only. `test/test.py`, `test/testAcoustics.py`.

## Runtime

**Environment:**
- MPI (OpenMPI or compatible) - distributed-memory parallelism across nodes. Compiler wrappers `mpicc` / `mpic++` set in `make.top`.
- OpenMP - optional shared-memory threading on CPU (enabled in release builds via `-fopenmp`).

**Package Manager:**
- None. All dependencies are either system libraries or a vendored git submodule (OCCA).
- Lockfile: Not applicable (no package manager).

## Frameworks

**HPC Compute Abstraction:**
- OCCA - GPU/accelerator portability layer. Bundled as a git submodule at `occa/` (commit `11552d0dc02fb9880f61f46e46115b6a50dada32`). Source: `https://github.com/libocca/occa`. Provides `device_t`, `kernel_t`, `deviceMemory<T>`, `pinnedMemory<T>`, and the OKL kernel language. OCCA targets CUDA, HIP, OpenCL, OpenMP, and Serial backends; the backend is selected at runtime via the `[THREAD MODEL]` setting in a `.rc` setup file.

**Build:**
- GNU Make - hierarchical makefiles. Top-level entry: `makefile`. Solver entry: `solvers/acoustics/makefile`. Shared flags/paths defined in `make.top` (included by all sub-makefiles).

**Testing:**
- Python 3 test runner (custom, no third-party framework). Tests launch solver binaries via `subprocess`, capture stdout, and compare the printed `Solution norm` against a hard-coded reference value within a tolerance of `1e-5`. Test file: `test/testAcoustics.py`. Invoked via `make test` from the solver or top-level directory.

## Key Dependencies

**Critical:**
- OCCA (submodule, version pinned by commit) - the entire GPU compute pipeline depends on it. Without OCCA the project cannot build. Located at `occa/`, linked as `-locca` from `$(OCCA_DIR)/lib`.
- MPI (system) - required for all distributed runs and even single-process builds because `mpic++` is the compiler wrapper. Linked implicitly.
- OpenBLAS (system) - linked as `-lopenblas` from `$(LIBP_BLAS_DIR)`. Default path `/usr/lib/x86_64-linux-gnu/openblas-serial`. Used by `parAlmond` and `linearSolver` layers; the acoustics solver does not call BLAS directly but links against the full `libp_libs` set.

**Infrastructure:**
- `libp::core` (`libs/core/`) - utility types, error macros (`LIBP_ABORT`, `LIBP_FORCE_ABORT`), `Factor2/3`, `RankDecomp2/3`.
- `libp::linAlg` (`libs/linAlg/`) - OCCA-backed vector operations (`innerProd`, `axpy`, `set`, etc.). Acoustics solver uses `innerProd` for norm computation.
- `libp::ogs` (`libs/ogs/`) - OCCA Gather/Scatter library for parallel halo exchange. Acoustics uses `ogs::halo_t` for trace halo communication and calls `ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add)` in setup.
- `libp::mesh` (`libs/mesh/`) - DG mesh abstraction (element geometry, quadrature, connectivity arrays, mass matrix, plot interpolation).
- `libp::timeStepper` (`libs/timeStepper/`) - explicit time integrators (`ab3`, `lserk4`, `dopri5` and PML/multi-rate variants).
- `libp::parAdogs` (`libs/parAdogs/`) - parallel mesh partitioning/reordering library.

## Configuration

**Environment:**
- No runtime environment variables are required by the solver itself.
- OCCA backend selection (`CUDA`, `HIP`, `OpenCL`, `OpenMP`, `Serial`) is set in the `.rc` setup file under `[THREAD MODEL]`.
- `OCCA_DIR` is set from `make.top` and must point to the OCCA submodule at build time.
- `LIBP_BLAS_DIR` can be overridden to point to a non-default BLAS installation.

**Build:**
- `make.top` - single authoritative file for all compiler flags, paths, and architecture selection (`LIBP_ARCH=GNU` or `INTEL`).
- Release build: `-fopenmp -O3 -mavx2 -ftree-vectorize -march=native` (GNU).
- Debug build: `make debug=1` → `-O0 -g -Wall`.
- Coverage build: `LIBP_COVERAGE=1` → adds `--coverage -fprofile-abs-path`.

## Platform Requirements

**Development:**
- C++17-capable compiler accessible as `mpic++` (e.g., GCC 8+ or Intel 19+).
- OCCA submodule must be initialized (`git submodule update --init`).
- OpenBLAS development headers and shared library on the host.
- MPI runtime (OpenMPI or MPICH).

**Production:**
- Any system supported by OCCA: x86 Linux/macOS with CUDA GPUs, AMD GPUs (HIP), or CPU-only (OpenMP/Serial).
- MPI runtime for multi-rank execution.
- `acousticsMain` binary is invoked as: `mpirun -np N ./acousticsMain setupfile.rc`

---

*Stack analysis: 2026-04-22*
