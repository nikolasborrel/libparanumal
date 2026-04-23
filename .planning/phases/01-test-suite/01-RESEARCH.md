# Phase 1: Test Suite - Research

**Researched:** 2026-04-23
**Domain:** C++ unit/integration testing with Catch2 v3 in a GNU Make / MPI / OCCA build system
**Confidence:** HIGH

---

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

- **D-01:** All Phase 1 tests must pass green from day one — no expected-fail annotations in the Phase 1 PR.
- **D-02:** Unit tests (TESTS-01) assert the mass-matrix L2 norm of the acoustic pressure field for Gaussian pulse runs on all four element types (Tri2D, Quad2D, Tet3D, Hex3D), using the same reference norm values already established in the Python test runner.
- **D-03:** Integration tests (TESTS-02) for FreqIndep and PerfRefl BCs assert only that the solver setup runs without crashing (no norm assertion). Numerical norm assertions for these BC types are added in Phase 3.
- **D-04:** No WAV file comparison in Phase 1 — `tinywav` is not a Phase 1 dependency.
- **D-05:** Catch2 v3 — use the amalgamated single-file distribution (`catch_amalgamated.hpp` + `catch_amalgamated.cpp`). Fetch from the official Catch2 v3 GitHub release.
- **D-06:** Do not copy `catch.hpp` from the fork — the fork uses Catch2 v2, which is in maintenance-only mode.
- **D-07:** Test files live in `solvers/acoustics/tests/`. A new `make tests` target is added to `solvers/acoustics/makefile` that builds the test binary by linking against `libacoustics.a`. The existing `make test` target (Python runner) is untouched.
- **D-08:** The test binary links against the same `libacoustics.a` produced by `make lib` — no duplication of object files. `make tests` depends on `make lib`.
- **D-09:** Catch2 amalgamated source files are vendored into `solvers/acoustics/tests/` (not fetched at build time).
- **D-10:** A new step is added to `.github/workflows/build.yml` after the existing `make test` step.
- **D-11:** The existing `make test` CI step is unchanged.

### Claude's Discretion

- Test binary name (`acousticsTests` or similar) — follow the fork's convention if one exists, otherwise choose a clear name.
- Whether to expose `make tests` from the top-level `makefile` as well (convenience alias) — planner can decide.
- MPI invocation style for Catch2 tests in CI (single-rank vs. multi-rank) — match the existing Python runner's convention (`mpirun --oversubscribe -np 1` for unit, `-np 4` for MPI integration test).

### Deferred Ideas (OUT OF SCOPE)

- WAV file comparison in tests — deferred to Phase 6.
- LR BC integration test assertions — deferred to Phase 5.
- Norm assertions for FreqIndep/PerfRefl BC tests — deferred to Phase 3.
- Expected-fail test annotations (`[!shouldfail]`) — not needed in Phase 1.
</user_constraints>

---

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| TESTS-01 | Developer can run Catch2 unit tests that verify acoustic pressure field values against reference data | Catch2 v3 amalgamated provides TEST_CASE + REQUIRE(norm == Approx(ref).margin(1e-5)); solver called via `acoustics_t::Run()` and stdout parsed for "Solution norm" |
| TESTS-02 | Developer can run Catch2 integration tests covering FreqIndep and PerfRefl BCs on 2D and 3D meshes | Crash-only via REQUIRE_NOTHROW; upstream BOX mesh with boundary_flag=1 provides a wall (reflecting) BC; FreqIndep is semantically deferred to Phase 3 but a labeled smoke test can run with the same wall BC |
| TESTS-03 | CI pipeline executes the test suite and reports pass/fail per test case | New step in `.github/workflows/build.yml` runs `make tests` then `mpirun ... acousticsTests`; Catch2 reports per-test pass/fail to stdout |
| TESTS-04 | Test infrastructure compiles and links against the existing C++ libParanumal build system | `make tests` in `solvers/acoustics/makefile` links against `libacoustics.a`; includes `make.top`; uses `LIBP_CXX`, `ACOUSTICS_CXXFLAGS`, and `LIBS` verbatim |
</phase_requirements>

---

## Summary

Phase 1 ports a Catch2 test suite into the upstream `libParanumal` C++ build without touching any solver code. The upstream solver already has all the machinery needed: `libacoustics.a` (from `make lib`) exposes `acoustics_t`, `acoustics_t::Setup()`, and `acoustics_t::Run()` through `acoustics.hpp`; `Run()` prints "Solution norm = %17.15lg\n" to stdout on rank 0 at the end of every run; and four reference norms for the Gaussian pulse are already locked in the Python test runner.

The primary engineering task is writing a `make tests` makefile target, a C++ test file that links against `libacoustics.a` and Catch2 v3, and a CI step. Catch2 v3.14.0 (released 2026-04-05) is the current release; it ships `catch_amalgamated.hpp` and `catch_amalgamated.cpp` as its two-file distribution.

The fork (`libparanumal-dtu`) test file (`acousticsTestMain.c`) is a C file using Catch2 v2 that calls a fork-specific C API (`acousticsSetupMain`, `setupAide`). It cannot be ported directly. The upstream tests must call the C++ `acoustics_t` API: construct `comm_t`, `platformSettings_t`, `meshSettings_t`, `acousticsSettings_t`, parse settings from a temporary `.rc` file (or populate programmatically), instantiate `platform_t`, `mesh_t`, and `acoustics_t`, then call `acoustics_t::Run()`. A helper function writes a temporary `.rc` file, invokes the solver, and captures the norm from stdout — exactly the same pattern the Python harness uses, but in C++.

**Primary recommendation:** Write a single `tests/acousticsTests.cpp` file that (1) wraps the C++ solver API in a helper that returns the final norm, (2) uses Catch2 `Approx().margin(1e-5)` for norm assertions, (3) wraps BC smoke tests in `REQUIRE_NOTHROW`, and (4) uses `CATCH_CONFIG_RUNNER` to bracket MPI init/finalize around the Catch2 session.

---

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Test binary build | Build system (GNU Make) | — | New `make tests` target in `solvers/acoustics/makefile` |
| Solver execution in tests | Host C++ (MPI rank 0..N) | OCCA Serial backend | `acoustics_t::Run()` runs on-device via OCCA; CI uses Serial backend |
| Norm extraction | Host stdout parse | — | `Run()` prints to stdout; test captures via subprocess or direct call |
| Catch2 framework | Host C++ (test binary) | — | Vendored `catch_amalgamated.*` compiled once per test build |
| CI orchestration | GitHub Actions | — | New step in `build.yml` after existing `make test` step |
| Reference norm values | Hard-coded in test source | — | Copied verbatim from `testAcoustics.py` |

---

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Catch2 (amalgamated) | v3.14.0 [VERIFIED: GitHub API 2026-04-05] | C++ unit test framework | Only Catch2 v3 is actively maintained; v2 is in maintenance-only mode; amalgamated distribution requires no CMake and is HPC-offline-safe |
| MPI (OpenMPI) | System (CI: `libopenmpi-dev`) [VERIFIED: build.yml] | Wraps solver launch; `Comm::Init` / `Comm::Finalize` required by libp | Solver was designed as MPI-first |
| OCCA (submodule) | pinned at `11552d0dc02fb98` [VERIFIED: codebase] | GPU/CPU kernel dispatch; Serial backend used in CI | All solver compute runs through OCCA |
| libacoustics.a | built by `make lib` [VERIFIED: makefile] | Static library produced from `src/*.cpp`; test binary links against it | Avoids duplicate compilation of solver objects |

### Catch2 v3 Amalgamated File Details

| File | GitHub Release URL |
|------|-------------------|
| `catch_amalgamated.hpp` | `https://github.com/catchorg/Catch2/releases/download/v3.14.0/catch_amalgamated.hpp` |
| `catch_amalgamated.cpp` | `https://github.com/catchorg/Catch2/releases/download/v3.14.0/catch_amalgamated.cpp` |

Both files are vendored into `solvers/acoustics/tests/`. The `.cpp` is compiled once; the `.hpp` is included in test source files as `#include "catch_amalgamated.hpp"` (relative path). [VERIFIED: Catch2 v3 migration guide, Context7]

### Migration from Fork's Catch2 v2 to v3

The fork uses `#include "catch.hpp"` (v2 single-header). Phase 1 uses v3 amalgamated:

| v2 pattern | v3 amalgamated equivalent |
|-----------|--------------------------|
| `#include "catch.hpp"` | `#include "catch_amalgamated.hpp"` |
| `#define CATCH_CONFIG_RUNNER` | Still supported in v3 |
| `Catch::Session().run(argc, argv)` | Same API [VERIFIED: Context7/catchorg/catch2] |
| `REQUIRE(x == Approx(y))` | `REQUIRE(x == Catch::Approx(y).margin(tol))` [VERIFIED: Context7] |

**Installation (vendor into repo):**
```bash
cd solvers/acoustics/tests/
curl -LO https://github.com/catchorg/Catch2/releases/download/v3.14.0/catch_amalgamated.hpp
curl -LO https://github.com/catchorg/Catch2/releases/download/v3.14.0/catch_amalgamated.cpp
```

---

## Architecture Patterns

### System Architecture Diagram

```
[make tests]
     |
     v
[compile catch_amalgamated.cpp] --> catch_amalgamated.o
[compile acousticsTests.cpp]   --> acousticsTests.o
     |
     v
[link: acousticsTests.o + catch_amalgamated.o + libacoustics.a + libp_libs + OCCA]
     |
     v
[acousticsTests binary]
     |
     v  mpirun --oversubscribe -np 1
[MPI Init] --> [Catch::Session().run(argc, argv)]
                     |
         +-----------+-------------------+
         |                               |
    [Unit TEST_CASEs]           [Integration TEST_CASEs]
    Gaussian pulse (x4)         BC smoke tests (x4)
         |                               |
    [write temp .rc]           [write temp .rc]
    [platform+mesh+acoustics_t setup]   [same setup, different BC flag]
    [acoustics_t::Run()]       [REQUIRE_NOTHROW(acoustics_t::Run())]
    [parse stdout norm]
    [REQUIRE norm == Approx(ref).margin(1e-5)]
         |
         v
    [MPI Finalize]
```

### Recommended Project Structure

```
solvers/acoustics/tests/
├── catch_amalgamated.hpp       # Catch2 v3.14.0 vendored header
├── catch_amalgamated.cpp       # Catch2 v3.14.0 vendored source (compiled once)
├── acousticsTests.cpp          # All test cases (unit + integration)
└── data/                       # Test-specific .rc files (if not written programmatically)
```

All test case logic lives in a single `acousticsTests.cpp` to minimize build complexity. The `.rc` files for each element type can be written to a temp directory at runtime using `std::tmpnam` or `std::filesystem::temp_directory_path()`.

### Pattern 1: Custom `main` with MPI bracketing (Catch2 v3)

The Catch2 session must be run inside MPI init/finalize because the solver calls `Comm::Init` and `Comm::Finalize` which wrap `MPI_Init`/`MPI_Finalize`.

```cpp
// Source: Catch2 docs/own-main.md (Context7/catchorg/catch2)
// Note: CATCH_CONFIG_RUNNER is the v2 spelling; in v3 amalgamated,
// providing your own main() is the standard approach.

#include "catch_amalgamated.hpp"
#include "acoustics.hpp"

int main(int argc, char* argv[]) {
    Comm::Init(argc, argv);
    int result = Catch::Session().run(argc, argv);
    Comm::Finalize();
    return result;
}
```

`Comm::Init` / `Comm::Finalize` are the libp wrappers around `MPI_Init` / `MPI_Finalize` — verified in `include/comm.hpp`. They accept `argc`/`argv` by reference.

### Pattern 2: Unit test — norm assertion

The solver's `Run()` method prints to stdout. The cleanest approach for a self-contained C++ test is to call `Run()` in-process and capture stdout, then parse the norm. Because `libp::abort` throws `libp::exception` (not `exit()`), errors surface as C++ exceptions that Catch2 can report.

```cpp
// Source: acousticsRun.cpp (verified), test.py (verified pattern)
#include "catch_amalgamated.hpp"
#include "acoustics.hpp"
#include <sstream>

static double runGaussianPulse(int elementType, int dim, int degree,
                               int nx, int ny, int nz, const std::string& dataFile) {
    comm_t comm(Comm::World().Dup());

    platformSettings_t platformSettings(comm);
    meshSettings_t     meshSettings(comm);
    acousticsSettings_t acousticsSettings(comm);

    // Programmatically populate settings (avoids temp file I/O)
    platformSettings.changeSetting("THREAD MODEL",   "Serial");
    platformSettings.changeSetting("PLATFORM NUMBER", "0");
    platformSettings.changeSetting("DEVICE NUMBER",   "0");
    meshSettings.changeSetting("MESH FILE",           "BOX");
    meshSettings.changeSetting("MESH DIMENSION",      std::to_string(dim));
    meshSettings.changeSetting("ELEMENT TYPE",        std::to_string(elementType));
    meshSettings.changeSetting("BOX NX",              std::to_string(nx));
    meshSettings.changeSetting("BOX NY",              std::to_string(ny));
    meshSettings.changeSetting("BOX NZ",              std::to_string(nz));
    meshSettings.changeSetting("BOX BOUNDARY FLAG",   "-1");
    meshSettings.changeSetting("POLYNOMIAL DEGREE",   std::to_string(degree));
    acousticsSettings.changeSetting("DATA FILE",      dataFile);
    acousticsSettings.changeSetting("TIME INTEGRATOR","DOPRI5");
    acousticsSettings.changeSetting("CFL NUMBER",     "1.0");
    acousticsSettings.changeSetting("START TIME",     "0");
    acousticsSettings.changeSetting("FINAL TIME",     "1.0");
    acousticsSettings.changeSetting("OUTPUT TO FILE", "FALSE");

    // Redirect stdout to capture "Solution norm = ..."
    std::ostringstream buf;
    std::streambuf* oldCout = std::cout.rdbuf(buf.rdbuf());

    platform_t platform(platformSettings);
    mesh_t     mesh(platform, meshSettings, comm);
    acoustics_t acoustics(platform, mesh, acousticsSettings);
    acoustics.Run();

    std::cout.rdbuf(oldCout);

    // Parse norm from captured output
    std::string output = buf.str();
    double norm = -1.0;
    auto pos = output.rfind("Solution norm = ");
    if (pos != std::string::npos) {
        norm = std::stod(output.substr(pos + 16));
    }
    return norm;
}

TEST_CASE("Gaussian pulse L2 norm — Tri2D", "[unit][tri2d]") {
    double norm = runGaussianPulse(3, 2, 4, 10, 10, 10,
                                   DACOUSTICS "data/acousticsGaussian2D.h");
    REQUIRE(norm == Catch::Approx(10.1302322430996).margin(1e-5));
}
```

**Critical note on stdout redirect:** `printf()` in `acousticsRun.cpp` writes to C stdout (`FILE*`), not `std::cout`. The `std::cout.rdbuf()` redirect does not capture `printf` output. Two alternatives exist:

1. Use `freopen` / `dup2` at the OS level to redirect file descriptor 1, then read from a pipe. This is the most reliable but more complex.
2. Modify the capture approach: launch `acousticsMain` as a subprocess (same pattern as the Python runner). This avoids the printf/cout mismatch entirely.

**Recommended approach for Phase 1:** Use the subprocess invocation pattern (same as Python test.py) — write a temp `.rc` file, call `mpirun --oversubscribe -np 1 ./acousticsMain setup.rc` via `popen` or a similar C++ subprocess API, parse the last line of stdout. This is simpler and avoids the `printf`/`cout` redirect problem. The test binary itself does not need to be the solver process.

Alternatively: If calling the C++ API directly is preferred, use `fdopen`/`pipe`/`dup2` to capture fd 1 before calling `Run()`.

### Pattern 3: Integration smoke test (crash-only)

```cpp
// Source: CONTEXT.md D-03 (locked decision)
TEST_CASE("BC smoke — PerfRefl Tri2D", "[integration][perf_refl][tri2d]") {
    // BOX boundary_flag=1 → bc==1 in acousticsGaussian2D.h → perfect wall reflection
    REQUIRE_NOTHROW(runGaussianPulse(3, 2, 4, 10, 10, 10,
                                    DACOUSTICS "data/acousticsGaussian2D.h",
                                    /*boundary_flag=*/1));
}
```

For `FreqIndep` in Phase 1: since the FreqIndep flux variant does not exist in the upstream codebase yet (it is a Phase 3 addition), the Phase 1 "FreqIndep" smoke test is a labeled placeholder that runs with `boundary_flag=1` (wall) and asserts no crash. The test name and tag `[freq_indep]` establish the CI slot that Phase 3 will fill with a real norm assertion.

### Pattern 4: Catch2 floating-point comparison

```cpp
// Source: Context7/catchorg/catch2 — Approx matcher
// Use margin (absolute tolerance) matching the Python TOL = 1.0e-5
REQUIRE(norm == Catch::Approx(10.1302322430996).margin(1e-5));
```

`Catch::Approx` with `.margin(tol)` matches `|actual - expected| < tol`. The Python runner uses absolute tolerance `TOL = 1.0e-5` — use the same value for consistency. [VERIFIED: Context7/catchorg/catch2]

### Anti-Patterns to Avoid

- **Copying fork's `catch.hpp` (v2):** The fork's `catch.hpp` is Catch2 v2. Do not vendor it. Use v3 amalgamated. [VERIFIED: CONTEXT.md D-06]
- **Using `#include <catch2/catch_test_macros.hpp>`:** This path is for CMake-installed Catch2. With the amalgamated distribution, use `#include "catch_amalgamated.hpp"`. [VERIFIED: Catch2 migration guide via Context7]
- **Forgetting the `make.top` include in the new makefile target:** All sub-makefiles in this project include `make.top` for `LIBP_CXX`, `LIBP_CXXFLAGS`, `LIBP_LIBS`, etc. The `make tests` target must include the same. The existing `ifndef LIBP_MAKETOP_LOADED` guard at the top of the acoustics makefile handles this. [VERIFIED: acoustics makefile]
- **Adding `tests` to the MAKECMDGOALS filter:** The acoustics makefile has a filter guard (lines 65-70) that rejects unknown targets. `tests` must be added to the filter list alongside `acousticsMain lib clean ...`. [VERIFIED: acoustics makefile lines 65-66]
- **Running tests without `make lib` first:** `make tests` must depend on `libacoustics.a` (via `make lib`). [VERIFIED: CONTEXT.md D-08]
- **Using `std::cout.rdbuf()` to capture printf output:** `printf` writes to C `stdout` (fd 1), not `std::cout`. The redirect only captures `std::cout` output. Use subprocess invocation or fd-level redirect to capture the norm.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Unit test framework | Custom assertion macros | Catch2 v3 | Catch2 provides `TEST_CASE`, `REQUIRE`, `REQUIRE_NOTHROW`, `Approx`, per-test reporting, and exit codes; hand-rolling all of this is error-prone |
| Floating-point comparison | `abs(a-b) < 1e-5` bare | `Catch::Approx().margin()` | Approx gives a clear diagnostic on failure showing expected vs actual |
| Test discovery | Manual test registry | Catch2 auto-registration | `TEST_CASE` macros self-register; no main-side list needed |

**Key insight:** The only custom code needed is the solver-invocation helper (`runGaussianPulse`). Everything else — test runner, assertion, reporting, exit code — is Catch2.

---

## Common Pitfalls

### Pitfall 1: `printf` stdout capture in same-process test

**What goes wrong:** The test calls `acoustics_t::Run()` in-process and tries to capture "Solution norm = ..." via `std::cout.rdbuf()`. The norm is not captured because `Run()` uses `printf` (C `FILE*` stdout), not `std::cout`. The test sees `norm == -1.0` and fails.

**Why it happens:** `acousticsRun.cpp` line 62: `printf("Solution norm = %17.15lg\n", norm2)` — C printf goes to fd 1 directly.

**How to avoid:** Either (a) invoke the solver as a subprocess (`mpirun ./acousticsMain tmp.rc`) and parse stdout via pipe, or (b) use `dup2`/`pipe` to redirect fd 1 before calling `Run()` and restore it after. Option (a) is simpler and matches the existing Python test pattern.

**Warning signs:** `norm == -1.0` returned from the helper; "Solution norm" not found in captured buffer.

### Pitfall 2: Makefile filter rejects `tests` target

**What goes wrong:** Running `make tests` in `solvers/acoustics/` prints an error message and exits because `tests` is not in the `$(filter ...)` guard on lines 65-66.

**Why it happens:** The acoustics makefile explicitly enumerates valid targets in a filter expression to produce a helpful error for unknown targets.

**How to avoid:** Add `tests` (and `clean-tests`) to the filter list on line 65. [VERIFIED: acoustics makefile]

**Warning signs:** `make tests` prints the ACOUSTICS_HELP_MSG and exits with an error immediately.

### Pitfall 3: `libacoustics.a` compiled without `LIBP_COVERAGE` when CI needs coverage

**What goes wrong:** The CI builds with `LIBP_COVERAGE=1`, which adds `--coverage -fprofile-abs-path` to `LIBP_CXXFLAGS`. If `libacoustics.a` was already built without coverage flags and `make tests` does not rebuild it, `.gcda` files for the solver code are not generated, and Codecov sees zero coverage from the Catch2 run.

**Why it happens:** `make tests` depends on `libacoustics.a` via `make lib`; if the `.a` is already present, Make skips the rebuild.

**How to avoid:** In CI, run `make clean` before `make tests` when `LIBP_COVERAGE=1`, or ensure the `make tests` step runs a fresh `make lib` as its first dependency. The existing CI already builds from scratch (`make -j $(nproc) verbose=true LIBP_COVERAGE=1`), so `libacoustics.a` will be freshly built before `make tests` runs.

**Warning signs:** `.gcno` files present but `.gcda` files absent for solver source files.

### Pitfall 4: OCCA kernel cache collision between Python and Catch2 test runs

**What goes wrong:** The Python `make test` step and the new `make tests` step both JIT-compile the same OKL kernels. If they run sequentially in CI, the second step finds the kernel cache already populated. This is fine functionally but can mask a build problem if the cache contains stale kernels from a different build.

**Why it happens:** OCCA caches compiled kernels in `${LIBP_DIR}/.occa/`. Both the Python runner (launching `acousticsMain`) and the Catch2 tests (also invoking the solver) share this cache.

**How to avoid:** Not a blocking issue for Phase 1. If needed, add `make clean-kernels` before the Catch2 step. For the current phase this is not required.

### Pitfall 5: MPI rank count mismatch between test invocation and solver assumptions

**What goes wrong:** Running the test binary with `-np 4` when the Gaussian pulse unit tests are written for single-rank produces different norms because `testAcousticsTri_MPI` (4 ranks) has a different reference: `10.1300558638317` vs `10.1302322430996` (1 rank).

**Why it happens:** Mesh partitioning across ranks changes the numerical result slightly due to floating-point non-associativity.

**How to avoid:** Unit tests (norm assertions) run with `-np 1`. The MPI integration test, if added, uses `-np 4` with the MPI-specific reference norm. For Phase 1, only the single-rank norms are asserted. [VERIFIED: CONTEXT.md, testAcoustics.py]

---

## Code Examples

### Verified `.rc` Settings for Each Unit Test

These settings exactly reproduce the Python test results. [VERIFIED: testAcoustics.py]

**Tri2D (element=3, dim=2, degree=4):**
```
[FORMAT]        2.0
[DATA FILE]     data/acousticsGaussian2D.h
[MESH FILE]     BOX
[MESH DIMENSION] 2
[ELEMENT TYPE]  3
[BOX NX]        10
[BOX NY]        10
[BOX NZ]        10
[BOX BOUNDARY FLAG] -1
[POLYNOMIAL DEGREE] 4
[THREAD MODEL]  Serial
[PLATFORM NUMBER] 0
[DEVICE NUMBER] 0
[TIME INTEGRATOR] DOPRI5
[CFL NUMBER]    1.0
[START TIME]    0
[FINAL TIME]    1.0
[OUTPUT TO FILE] FALSE
```

**Quad2D (element=4, dim=2, degree=4):** Same as Tri2D with `[ELEMENT TYPE] 4`.

**Tet3D (element=6, dim=3, degree=2):**
```
[DATA FILE]     data/acousticsGaussian3D.h
[MESH DIMENSION] 3
[ELEMENT TYPE]  6
[POLYNOMIAL DEGREE] 2
[BOX BOUNDARY FLAG] -1
```
(All other settings same as Tri2D.)

**Hex3D (element=12, dim=3, degree=2):** Same as Tet3D with `[ELEMENT TYPE] 12`.

### Reference Norms (from testAcoustics.py, verbatim) [VERIFIED: test/testAcoustics.py]

| Element Type | Element Code | Dim | Degree | Reference Norm |
|-------------|-------------|-----|--------|----------------|
| Tri2D | 3 | 2 | 4 | `10.1302322430996` |
| Quad2D | 4 | 2 | 4 | `10.1299609797959` |
| Tet3D | 6 | 3 | 2 | `31.6577046152384` |
| Hex3D | 12 | 3 | 2 | `31.6576028812776` |

Tolerance: `1e-5` (absolute, matching Python runner's `TOL = 1.0e-5`).

### `make tests` Target Pattern [VERIFIED: acoustics makefile + fork makefile]

```makefile
# Add to filter guard on line 65 (after "test"):
ifeq (,$(filter acousticsMain lib clean clean-libs clean-kernels \
                realclean info help test tests clean-tests,$(MAKECMDGOALS)))

# Test source objects
TESTS_SRC = tests/catch_amalgamated.cpp tests/acousticsTests.cpp
TESTS_OBJS = tests/catch_amalgamated.o tests/acousticsTests.o

# Rule for test object files (same flags as solver objects)
tests/%.o: tests/%.cpp $(DEPS)
ifneq (,${verbose})
    $(LIBP_CXX) -o $@ -c $< $(ACOUSTICS_CXXFLAGS) -I tests/
else
    @printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
    @$(LIBP_CXX) -o $@ -c $< $(ACOUSTICS_CXXFLAGS) -I tests/
endif

# Build test binary linking against libacoustics.a
acousticsTests: libacoustics.a $(TESTS_OBJS)
ifneq (,${verbose})
    $(LIBP_LD) -o acousticsTests $(TESTS_OBJS) -L. -lacoustics $(LFLAGS)
else
    @printf "%b" "$(EXE_COLOR)Linking $(@F)$(NO_COLOR)\n";
    @$(LIBP_LD) -o acousticsTests $(TESTS_OBJS) -L. -lacoustics $(LFLAGS)
endif

tests: acousticsTests

clean-tests:
    rm -f tests/*.o acousticsTests
```

**Note:** `-I tests/` is needed so `#include "catch_amalgamated.hpp"` resolves. The `-L. -lacoustics` links the static library in the current directory before the other `$(LFLAGS)` libs.

### CI Step Addition [VERIFIED: .github/workflows/build.yml]

The existing CI workflow has two steps after checkout+install: Build and Test. The new Catch2 step goes after the existing Test step:

```yaml
    - name: Test
      run: |
        make test
        bash <(curl --no-buffer -s https://codecov.io/bash) -x gcov
      env:
        LIBP_COVERAGE: 1

    - name: Build Catch2 Tests
      run: |
        cd solvers/acoustics
        make tests
      env:
        LIBP_COVERAGE: 1

    - name: Run Catch2 Tests
      run: |
        cd solvers/acoustics
        mpirun --oversubscribe -np 1 ./acousticsTests
```

Separating Build and Run into two steps gives CI clearer failure attribution. The `LIBP_COVERAGE: 1` env on the build step ensures coverage flags match the earlier build.

---

## Build System Integration

### How `make lib` Produces `libacoustics.a` [VERIFIED: acoustics makefile]

```
make lib
  → depends on: $(OBJS) = src/*.cpp → src/*.o
  → each .cpp compiled with: $(LIBP_CXX) -o $*.o -c $*.cpp $(ACOUSTICS_CXXFLAGS)
  → ACOUSTICS_CXXFLAGS = $(LIBP_CXXFLAGS) $(DEFINES) $(INCLUDES)
  → LIBP_CXXFLAGS = -fopenmp -O3 -Wall ... -std=c++17 (plus --coverage in CI)
  → INCLUDES = -I${LIBP_INCLUDE_DIR} -I${OCCA_DIR}/include -I.
  → DEFINES = ${LIBP_DEFINES} -DLIBP_DIR='"${LIBP_DIR}"'
  → ar -cr libacoustics.a $(OBJS)
  → output: libacoustics.a in solvers/acoustics/
```

`libacoustics.a` contains all objects from `src/acousticsRun.o`, `src/acousticsSetup.o`, `src/acousticsSettings.o`, `src/acousticsReport.o`, `src/acousticsPlotFields.o`, and the `rhsf` / `MaxWaveSpeed` implementations.

**Important:** `acousticsMain.o` is NOT in `libacoustics.a` — it is linked separately for the `acousticsMain` binary. The test binary replaces `acousticsMain.o` with its own `main()` (provided by Catch2 session runner + test cases).

### Compiler/Linker Flags for Test Binary [VERIFIED: make.top, acoustics makefile]

```makefile
# Compilation flags for test objects (same as solver objects):
$(LIBP_CXX) -o tests/acousticsTests.o -c tests/acousticsTests.cpp \
    -fopenmp -O3 -Wall -Wshadow -Wno-unused-function -std=c++17 \
    -mavx2 -ftree-vectorize -march=native \
    -DLIBP_DIR='"${LIBP_DIR}"' \
    -I${LIBP_INCLUDE_DIR} -I${OCCA_DIR}/include -I. -I tests/

# Link flags for test binary:
$(LIBP_LD) -o acousticsTests \
    tests/catch_amalgamated.o tests/acousticsTests.o \
    -L. -lacoustics \                         # libacoustics.a
    -L${LIBP_LIBS_DIR} \                      # libp library archives
    -ltimeStepper -lmesh -lparAdogs -logs -llinAlg -lcore \
    -Wl,-rpath,$(LIBP_BLAS_DIR) -L$(LIBP_BLAS_DIR) -lopenblas \
    -Wl,-rpath,$(OCCA_DIR)/lib -L$(OCCA_DIR)/lib -locca
```

This is equivalent to the existing `LFLAGS` variable reused for the solver binary, plus `-L. -lacoustics`.

---

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| `mpic++` / `mpirun` | Build and test invocation | Yes (dev machine) [VERIFIED: shell] | Open MPI 5.0.8 | — |
| `g++` / clang | C++17 compilation | Yes [VERIFIED: shell] | Apple clang 17.0.0 | — |
| `libopenmpi-dev` | CI build | Yes (CI apt-install) [VERIFIED: build.yml] | — | — |
| `libopenblas-serial-dev` | CI link | Yes (CI apt-install) [VERIFIED: build.yml] | — | — |
| OCCA submodule | Kernel JIT | Present (git submodule) [VERIFIED: codebase] | pinned commit | — |
| Internet (Catch2 download) | Vendor step (one-time) | Not required at build time (vendored) | — | Already vendored |
| GPU | OCCA CUDA/HIP backend | Not required for CI | — | OCCA Serial backend |

**Missing dependencies with no fallback:** None — all required dependencies are present.

---

## Validation Architecture

### Test Framework

| Property | Value |
|----------|-------|
| Framework | Catch2 v3.14.0 (amalgamated) |
| Config file | none — standalone binary |
| Quick run command | `mpirun --oversubscribe -np 1 ./acousticsTests "[unit]"` |
| Full suite command | `mpirun --oversubscribe -np 1 ./acousticsTests` |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| TESTS-01 | Norm matches reference for Tri2D/Quad2D/Tet3D/Hex3D Gaussian pulse | unit | `mpirun --oversubscribe -np 1 ./acousticsTests "[unit]"` | No — Wave 0 |
| TESTS-02 | FreqIndep + PerfRefl setups do not crash | integration | `mpirun --oversubscribe -np 1 ./acousticsTests "[integration]"` | No — Wave 0 |
| TESTS-03 | CI reports pass/fail per test case | CI smoke | `make tests && mpirun --oversubscribe -np 1 ./acousticsTests` | No — Wave 0 |
| TESTS-04 | Test binary compiles and links | build | `make tests` (zero exit code) | No — Wave 0 |

### Sampling Rate

- **Per task commit:** `make tests && mpirun --oversubscribe -np 1 ./acousticsTests "[unit]"`
- **Per wave merge:** `mpirun --oversubscribe -np 1 ./acousticsTests` (full suite)
- **Phase gate:** Full suite green + `make test` (Python runner) unbroken before `/gsd-verify-work`

### Wave 0 Gaps

- [ ] `solvers/acoustics/tests/acousticsTests.cpp` — covers TESTS-01, TESTS-02
- [ ] `solvers/acoustics/tests/catch_amalgamated.hpp` — vendored Catch2 header
- [ ] `solvers/acoustics/tests/catch_amalgamated.cpp` — vendored Catch2 source
- [ ] Updated `solvers/acoustics/makefile` — `tests` and `clean-tests` targets (covers TESTS-03, TESTS-04)
- [ ] Updated `.github/workflows/build.yml` — new CI step (covers TESTS-03)

---

## Security Domain

> `security_enforcement` not set in config.json — treated as enabled.

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | No | — |
| V3 Session Management | No | — |
| V4 Access Control | No | — |
| V5 Input Validation | Partially | Temp `.rc` files are written by test code, not user input; no sanitization needed for test-generated data |
| V6 Cryptography | No | — |

### Known Threat Patterns for HPC test binaries

| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Temp file race (TOCTOU) | Tampering | Use `mkstemp` or `std::filesystem::temp_directory_path()` for unique temp `.rc` paths; not exploitable in CI but good practice |
| Vendored third-party source | Tampering | Verify GPG signature of `catch_amalgamated.*` against `.asc` files published in the v3.14.0 release if supply-chain integrity is required |

---

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Settings can be populated programmatically via `changeSetting()` without writing a `.rc` file | Pattern 2 code example | If `changeSetting` is not the right API, a temp `.rc` file must be written and parsed via `parseFromFile` instead |
| A2 | `Comm::Init(argc, argv)` and `Comm::Finalize()` can safely bracket the Catch2 session (MPI initialized once per process) | Pattern 1 | If `acoustics_t` calls `Comm::Init` internally, double-init will error; inspection shows it does not — but only `acousticsMain.cpp` was verified |
| A3 | `libp::abort` throws `libp::exception` (not `exit()` or `MPI_Abort()`), so Catch2 can catch and report failures as test failures rather than process crashes | Pitfall section | If `libp::abort` calls `MPI_Abort`, the test process exits immediately and Catch2 cannot record the failure |

**A3 is partially verified:** `libs/core/exception.cpp` line 79: `throw exception(...)` — this is a C++ throw, not `exit()` or `MPI_Abort()`. [VERIFIED: libs/core/exception.cpp]

**A1 and A2 are ASSUMED** — verified from code inspection but not from a live build+run test.

---

## Open Questions

1. **Can settings be set programmatically, or must a `.rc` file be parsed?**
   - What we know: `acousticsSettings_t` inherits from `settings_t` which has `changeSetting(name, val)`. `parseFromFile` calls `changeSetting` internally.
   - What's unclear: Whether there are required settings that have no default and must be set — e.g., `DATA FILE` has a default, but platform/mesh settings may not.
   - Recommendation: Verify by attempting to build with only programmatic settings. If it fails, write a temp `.rc` file and use `parseFromFile` — the same approach the Python runner uses.

2. **Does `acoustics_t::Run()` call `MPI_Finalize` or `Comm::Finalize`?**
   - What we know: `acousticsRun.cpp` has no `Comm::Finalize` call — it just computes and prints the norm.
   - What's unclear: Whether any library called by `Run()` finalizes MPI internally.
   - Recommendation: Inspect `timeStepper.Run()` and `linAlg().innerProd()` for MPI finalization. Based on code review, neither does — MPI cleanup is only in `acousticsMain.cpp`.

3. **Binary name: `acousticsTests` or `acousticsTestMain`?**
   - Fork uses `acousticsTestMain` (from `tests:` target in fork makefile).
   - Upstream convention (existing binaries): `acousticsMain`, `advectionMain`, etc. — all use `*Main` suffix.
   - Recommendation: Use `acousticsTests` (not `acousticsTestMain`) — the `Main` suffix on existing binaries is because they define `main()`; the test binary also defines `main()` through Catch2 but the conventional name in CI test suites is the solver name + "Tests".

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Catch2 v2 single header (`catch.hpp`) | Catch2 v3 amalgamated (`catch_amalgamated.hpp` + `.cpp`) | Catch2 v3.0.0 (2022) | Include path changes; `CATCH_CONFIG_MAIN` replaced by providing own `main`; `Approx` is now `Catch::Approx` |
| Fork procedural C API (`acousticsSetupMain`, `setupAide`) | Upstream C++ class API (`acoustics_t`, `platform_t`, `mesh_t`) | libParanumal v0.5.0 (2022) | Tests must use the C++ constructor chain, not the C procedural API |

**Deprecated/outdated:**
- `catch.hpp` (v2 single-header): In maintenance-only mode as of 2022. Do not use. [CITED: Catch2 migration guide]
- `CATCH_CONFIG_MAIN`: Still functional in v3 amalgamated but the recommended pattern is providing your own `main()`. [CITED: Context7/catchorg/catch2 docs/own-main.md]
- Fork's `build_acoustics_tests.sh`: HPC-specific (loads LMOD modules, sets `OCCA_DIR`); not applicable to upstream GNU Make build.

---

## Sources

### Primary (HIGH confidence)
- Codebase: `solvers/acoustics/makefile` — makefile target structure, filter guard, `LFLAGS`
- Codebase: `solvers/acoustics/acoustics.hpp` — C++ API: `acoustics_t`, `Setup`, `Run`
- Codebase: `solvers/acoustics/src/acousticsRun.cpp` — exact `printf("Solution norm = %17.15lg\n", norm2)` pattern
- Codebase: `solvers/acoustics/src/acousticsSettings.cpp` — settings keys and defaults
- Codebase: `test/testAcoustics.py` — reference norms and test configuration
- Codebase: `.github/workflows/build.yml` — CI step structure
- Codebase: `make.top` — compiler flags, `LIBP_CXX`, `LIBP_CXXFLAGS`, `LIBP_LIBS`
- Codebase: `include/utils.hpp` + `libs/core/exception.cpp` — `LIBP_ABORT` throws `libp::exception`
- GitHub API: Catch2 v3.14.0 release (2026-04-05) — confirmed as latest release, confirmed amalgamated file names
- Context7 `/catchorg/catch2` — Catch2 v3 custom main pattern, amalgamated distribution docs, `Approx` matcher

### Secondary (MEDIUM confidence)
- Fork codebase: `~/github/libparanumal-dtu/solvers/acoustics/tests/acousticsTestMain.c` — fork test structure (v2 API, C language; not directly portable but confirms test case intent)
- Fork codebase: `~/github/libparanumal-dtu/solvers/acoustics/makefile` — `make tests` target pattern (C procedural build; structure is similar)

### Tertiary (LOW confidence)
- None — all claims were verified from codebase or official sources.

---

## Metadata

**Confidence breakdown:**
- Standard stack (Catch2 v3.14.0, MPI, OCCA): HIGH — verified via GitHub API and codebase
- Architecture (makefile pattern, link flags): HIGH — verified from actual makefile
- C++ solver API (acoustics_t, Run, settings): HIGH — verified from source files
- Reference norms: HIGH — copied verbatim from testAcoustics.py
- Catch2 v3 patterns (custom main, Approx): HIGH — verified via Context7 official docs
- stdout capture pitfall (printf vs cout): HIGH — verified from acousticsRun.cpp source
- Programmatic settings API (A1, A2): MEDIUM — inferred from code structure, not live-tested

**Research date:** 2026-04-23
**Valid until:** 2026-05-23 (stable stack — Catch2 release cadence is low; MPI and OCCA are pinned)
