# Testing Patterns — acoustics solver

**Analysis Date:** 2026-04-22
**Scope:** `solvers/acoustics/` within the libParanumal HPC finite element library

---

## Test Framework

**Type:** Integration / regression tests only. There is no unit test framework (no GoogleTest, Catch2, etc.).

**Test runner:** Python 3 scripts in `/Users/nikolasborrel/github/libParanumal/test/`

**Invocation:**
```bash
# From the acoustics solver directory:
make test

# From the repo root test/ directory:
cd test/
./testAcoustics.py          # acoustics only
./test.py                   # all solvers (includes acoustics)

# From the root (runs all solvers):
make test
```

The acoustics `make test` target delegates to the shared test harness:
```makefile
# solvers/acoustics/makefile
test: acousticsMain
    @${MAKE} -C $(LIBP_TEST_DIR) --no-print-directory test-acoustics

# test/makefile
test-acoustics:
    @./testAcoustics.py
```

---

## What Tests Exist

All acoustics tests live in `/Users/nikolasborrel/github/libParanumal/test/testAcoustics.py`.

### Test Cases

| Test name | Element type | Dim | Polynomial degree | Mesh | Ranks | Reference norm |
|---|---|---|---|---|---|---|
| `testAcousticsTri` | Triangles (3) | 2D | 4 | BOX 10×10 | 1 | `10.1302322430996` |
| `testAcousticsQuad` | Quadrilaterals (4) | 2D | 4 | BOX 10×10 | 1 | `10.1299609797959` |
| `testAcousticsTet` | Tetrahedra (6) | 3D | 2 | BOX 10×10×10 | 1 | `31.6577046152384` |
| `testAcousticsHex` | Hexahedra (12) | 3D | 2 | BOX 10×10×10 | 1 | `31.6576028812776` |
| `testAcousticsTri_MPI` | Triangles (3) | 2D | 4 | BOX 10×10 | 4 | `10.1300558638317` |

All tests use:
- Initial condition: Gaussian pressure pulse `r = 1 + exp(-3*(x²+y²[+z²]))`, zero velocity
- Final time: `t = 1.0`
- Time integrator: DOPRI5
- CFL: 1.0
- Boundary flag: `-1` (default, periodic/no boundary in test; solver setups use flag `1` = wall)

The `testAcousticsTri_MPI` test additionally checks that VTU output files are produced (`OUTPUT TO FILE=TRUE`), then cleans up all `.vtu` files from the test directory.

---

## Test Mechanism

### How Tests Work

The test harness (`/Users/nikolasborrel/github/libParanumal/test/test.py`) works as follows:

1. A `setup.rc` configuration file is written to the test directory via `writeSetup("setup", settings)`.
2. The solver binary is invoked via `mpirun --oversubscribe -np {ranks} {acousticsBin} setup.rc`.
3. The **last line of stdout** is parsed for the pattern `"Solution norm = "` followed by a floating-point value.
4. The parsed norm is compared against a hardcoded `referenceNorm` with tolerance `TOL = 1.0e-5`.
5. If the norm matches, the test prints `PASS`; otherwise `FAIL` with expected vs. observed values.
6. On failure, the full stdout and stderr are printed, and the failing `.rc` file is saved for reproducibility.

### Binary Location
```python
acousticsBin = solverDir + "/acoustics/acousticsMain"
```
The binary must be built before running tests (`make acousticsMain` or `make` in `solvers/acoustics/`).

### Solution Norm Computation
The `Solution norm` printed at the end of `acousticsRun.cpp` is the mass-matrix-weighted L2 norm of the full state vector `q` at the final time:
```cpp
mesh.MassMatrixApply(o_q, o_Mq);
dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm));
printf("Solution norm = %17.15lg\n", norm2);
```
This is not a convergence error against an exact solution — it is the energy norm of the numerical solution itself. The test therefore checks that the solver reproduces the same numerical answer bitwise to `1e-5` relative tolerance, not that it is mathematically correct.

---

## CI Pipeline

**Config:** `/Users/nikolasborrel/github/libParanumal/.github/workflows/build.yml`

**Trigger:** Push or pull request to `master` branch; also `workflow_dispatch`.

**Runner:** `ubuntu-latest` (GitHub-hosted, CPU only)

**Steps:**
1. Checkout with submodules (`occa` is a submodule)
2. Install dependencies: `libopenmpi-dev openmpi-bin libopenblas-serial-dev`
3. Build: `make -j $(nproc) verbose=true` with `LIBP_COVERAGE=1`
4. Test: `make test` (runs all solvers including acoustics) then uploads coverage to Codecov via `gcov`

**Thread model in CI:** The CI does not have a GPU, so all tests run in `Serial` mode (OCCA falls back to CPU). The `device` variable in the test Python scripts defaults to `"Serial"` when no argument is provided.

**Coverage:** `LIBP_COVERAGE=1` adds `--coverage -fprofile-abs-path` to `LIBP_CXXFLAGS`. Coverage collection is via `gcov` reported to Codecov.

---

## Test Data Files

| File | Purpose |
|---|---|
| `/Users/nikolasborrel/github/libParanumal/solvers/acoustics/data/acousticsGaussian2D.h` | BCs and ICs for 2D tests. Defines `acousticsDirichletConditions2D` and `acousticsInitialConditions2D` macros. |
| `/Users/nikolasborrel/github/libParanumal/solvers/acoustics/data/acousticsGaussian3D.h` | BCs and ICs for 3D tests. Defines `acousticsDirichletConditions3D` and `acousticsInitialConditions3D` macros. |
| `/Users/nikolasborrel/github/libParanumal/test/squareTri.msh` | Gmsh triangle mesh (not used by acoustics tests — acoustics uses BOX generator) |
| `/Users/nikolasborrel/github/libParanumal/test/squareQuad.msh` | Gmsh quad mesh (not used by acoustics tests) |
| `/Users/nikolasborrel/github/libParanumal/test/cubeTet.msh` | Gmsh tet mesh (not used by acoustics tests) |
| `/Users/nikolasborrel/github/libParanumal/test/cubeHex.msh` | Gmsh hex mesh (not used by acoustics tests) |

The acoustics tests always use the built-in BOX mesh generator (`MESH FILE = BOX`), not external mesh files.

---

## Setup Files (Manual / Interactive Use)

The `setups/` directory contains reference `.rc` files for manual runs, not automated tests:

| File | Element type | Dim | Default integrator | Final time |
|---|---|---|---|---|
| `setups/setupTri2D.rc` | Triangles | 2D | DOPRI5 | 10 |
| `setups/setupQuad2D.rc` | Quadrilaterals | 2D | DOPRI5 | 10 |
| `setups/setupTet3D.rc` | Tetrahedra | 3D | DOPRI5 | 10 |
| `setups/setupHex3D.rc` | Hexahedra | 3D | DOPRI5 | 10 |

These files set `OUTPUT TO FILE=TRUE` and `THREAD MODEL=CUDA`, targeting interactive GPU runs with t_final=10. They are not consumed by `make test`.

---

## What Is NOT Tested

- **Convergence order** — No hp-refinement study. The tests verify a fixed norm value, not p-th order convergence rate.
- **Exact solution comparison** — No manufactured solution or analytic reference. The Gaussian pulse has no closed-form solution at t>0.
- **Neumann (absorbing) boundary conditions** — The surface kernel contains a comment: `//should also add the Neumann BC here, but need uxM, uyM, vxM, and vyM somehow`. This code path is absent.
- **Performance / throughput** — No timing assertions.
- **GPU execution in CI** — CI runs in Serial (CPU) mode; CUDA/HIP/OpenCL paths are not exercised automatically.
- **VTK output correctness** — `testAcousticsTri_MPI` only confirms VTU files are created and deleted; file content is not validated.
- **High polynomial degree** — Tet and Hex tests use degree 2; only Tri and Quad tests use degree 4.
- **Single-process 3D** — The MPI test covers only 2D Tri with 4 ranks; no MPI test for 3D elements.
- **Different time integrators** — All tests use DOPRI5 only; AB3 and LSERK4 paths are untested in CI.
- **Non-box (unstructured) meshes** — All tests use BOX mesh generator; external Gmsh meshes in `test/` are not used for acoustics.

---

## Test Coverage Gaps (Priority)

**High:**
- No convergence test: implement p-refinement or h-refinement study with Gaussian pulse or manufactured solution to verify spectral/optimal convergence rates for each element type.
- Neumann BC path is unimplemented and untested (see comment in `acousticsSurfaceTri2D.okl` line 118–119 and `acousticsSurfaceHex3D.okl`).

**Medium:**
- AB3 and LSERK4 time integrators have no test coverage.
- 3D MPI test is absent.
- GPU execution is never exercised in CI.

**Low:**
- High-degree (p≥4) tests for 3D elements.
- External mesh file tests for acoustics.

---

*Testing analysis: 2026-04-22*
