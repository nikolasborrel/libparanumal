# Phase 1: Test Suite - Pattern Map

**Mapped:** 2026-04-23
**Files analyzed:** 5 (3 new, 2 modified)
**Analogs found:** 5 / 5

---

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|---|---|---|---|---|
| `solvers/acoustics/tests/acousticsTests.cpp` | test | request-response (subprocess invoke) | `test/testAcoustics.py` + `solvers/acoustics/acousticsMain.cpp` | role-match (Python→C++, same invoke pattern) |
| `solvers/acoustics/tests/catch_amalgamated.hpp` | config (vendored) | — | `~/github/libparanumal-dtu/solvers/acoustics/tests/catch.hpp` | version-upgrade (v2→v3) |
| `solvers/acoustics/tests/catch_amalgamated.cpp` | config (vendored) | — | none in upstream | no analog (new) |
| `solvers/acoustics/makefile` | config (build) | batch | `solvers/acoustics/makefile` (self — adding targets) | exact (modify in place) |
| `.github/workflows/build.yml` | config (CI) | batch | `.github/workflows/build.yml` (self — adding step) | exact (modify in place) |

---

## Pattern Assignments

### `solvers/acoustics/tests/acousticsTests.cpp` (test, request-response)

**Primary analog:** `test/testAcoustics.py` — subprocess invoke pattern, reference norms, tolerance
**Secondary analog:** `solvers/acoustics/acousticsMain.cpp` — C++ MPI init/finalize wrapper, solver construction chain

---

**MPI init/finalize wrapper pattern**
Source: `solvers/acoustics/acousticsMain.cpp` lines 29–68

```cpp
int main(int argc, char **argv){
  Comm::Init(argc, argv);

  { /*Scope so everything is destructed before MPI_Finalize */
    comm_t comm(Comm::World().Dup());

    platformSettings_t platformSettings(comm);
    meshSettings_t meshSettings(comm);
    acousticsSettings_t acousticsSettings(comm);

    acousticsSettings.parseFromFile(platformSettings, meshSettings, argv[1]);

    platform_t platform(platformSettings);
    mesh_t mesh(platform, meshSettings, comm);
    acoustics_t acoustics(platform, mesh, acousticsSettings);

    acoustics.Run();
  }

  Comm::Finalize();
  return LIBP_SUCCESS;
}
```

For the test binary, replace `parseFromFile(...)` with programmatic `changeSetting(...)` calls (see settings pattern below). Replace the bare `acoustics.Run()` call with the subprocess pattern or fd-redirect pattern described in RESEARCH.md pitfall 1.

**Catch2 custom main with MPI bracketing**
Source: RESEARCH.md Pattern 1 (verified from Catch2 docs + `acousticsMain.cpp` init pattern)

```cpp
#include "catch_amalgamated.hpp"
#include "acoustics.hpp"

int main(int argc, char* argv[]) {
    Comm::Init(argc, argv);
    int result = Catch::Session().run(argc, argv);
    Comm::Finalize();
    return result;
}
```

`Catch::Session().run(argc, argv)` replaces `acousticsMain`'s solver call. The MPI init/finalize bookends are identical to `acousticsMain.cpp` lines 32 and 66.

**Solver construction + settings pattern**
Source: `solvers/acoustics/acousticsMain.cpp` lines 37–62; settings keys from `test/testAcoustics.py` lines 34–56

```cpp
comm_t comm(Comm::World().Dup());

platformSettings_t platformSettings(comm);
meshSettings_t     meshSettings(comm);
acousticsSettings_t acousticsSettings(comm);

// Programmatic population — avoids temp file I/O (assumption A1 in RESEARCH.md)
platformSettings.changeSetting("THREAD MODEL",    "Serial");
platformSettings.changeSetting("PLATFORM NUMBER", "0");
platformSettings.changeSetting("DEVICE NUMBER",   "0");
meshSettings.changeSetting("MESH FILE",           "BOX");
meshSettings.changeSetting("MESH DIMENSION",      std::to_string(dim));
meshSettings.changeSetting("ELEMENT TYPE",        std::to_string(elementType));
meshSettings.changeSetting("BOX NX",              std::to_string(nx));
meshSettings.changeSetting("BOX NY",              std::to_string(ny));
meshSettings.changeSetting("BOX NZ",              std::to_string(nz));
meshSettings.changeSetting("BOX BOUNDARY FLAG",   std::to_string(boundaryFlag));
meshSettings.changeSetting("POLYNOMIAL DEGREE",   std::to_string(degree));
acousticsSettings.changeSetting("DATA FILE",      dataFile);
acousticsSettings.changeSetting("TIME INTEGRATOR","DOPRI5");
acousticsSettings.changeSetting("CFL NUMBER",     "1.0");
acousticsSettings.changeSetting("START TIME",     "0");
acousticsSettings.changeSetting("FINAL TIME",     "1.0");
acousticsSettings.changeSetting("OUTPUT TO FILE", "FALSE");

platform_t  platform(platformSettings);
mesh_t      mesh(platform, meshSettings, comm);
acoustics_t acoustics(platform, mesh, acousticsSettings);
acoustics.Run();
```

**Stdout norm output format (what to parse)**
Source: `solvers/acoustics/src/acousticsRun.cpp` line 63; `test/test.py` lines 126–130

`acousticsRun.cpp` line 63:
```c
printf("Solution norm = %17.15lg\n", norm2);
```

`test/test.py` parsing (lines 126–130):
```python
output = run.stdout.decode().splitlines()[-1]
if "Solution norm = " in output:
    norm = float(output.split()[3])
```

Critical: `printf` writes to C `FILE*` fd 1, not `std::cout`. `std::cout.rdbuf()` redirect does NOT capture it. Use subprocess invocation (write temp `.rc`, call `mpirun ./acousticsMain tmp.rc` via `popen`) or use `dup2`/`pipe` to redirect fd 1. The subprocess approach exactly mirrors `test/test.py` lines 101–160 (`subprocess.run(["mpirun", "--oversubscribe", "-np", str(ranks), cmd, inputRC], stdout=subprocess.PIPE, ...)`).

**Subprocess invocation helper pattern (recommended)**
Source: `test/test.py` lines 101–160 — write `.rc` then subprocess; adapt to C++ using `popen` or `FILE* pipe = popen("mpirun ...", "r")`.

```python
# Python pattern to adapt to C++:
run = subprocess.run(
    ["mpirun", "--oversubscribe", "-np", str(ranks), cmd, inputRC],
    stdout=subprocess.PIPE, stderr=subprocess.PIPE
)
output = run.stdout.decode().splitlines()[-1]
norm = float(output.split()[3])   # field index 3: "Solution norm = <value>"
```

In C++:
1. Write a temp `.rc` file to `std::filesystem::temp_directory_path()` using the same `[KEY]\nvalue\n` format as `test/test.py`'s `writeSetup` (lines 91–99).
2. Invoke `popen("mpirun --oversubscribe -np 1 " + binaryPath + " " + rcPath, "r")` and read stdout line by line.
3. Find the line containing `"Solution norm = "` and parse the double from field index 3.
4. Return the norm; the TEST_CASE asserts with `Catch::Approx(...).margin(1e-5)`.

**Unit test norm assertions**
Source: `test/testAcoustics.py` lines 61–81; tolerance from `test/test.py` line 63 (`TOL = 1.0e-5`)

```cpp
TEST_CASE("Gaussian pulse L2 norm - Tri2D", "[unit][tri2d]") {
    double norm = runGaussianPulse(/*element=*/3, /*dim=*/2, /*degree=*/4,
                                   /*nx=*/10, /*ny=*/10, /*nz=*/10,
                                   DACOUSTICS "data/acousticsGaussian2D.h",
                                   /*boundary_flag=*/-1);
    REQUIRE(norm == Catch::Approx(10.1302322430996).margin(1e-5));
}

TEST_CASE("Gaussian pulse L2 norm - Quad2D", "[unit][quad2d]") {
    double norm = runGaussianPulse(4, 2, 4, 10, 10, 10,
                                   DACOUSTICS "data/acousticsGaussian2D.h", -1);
    REQUIRE(norm == Catch::Approx(10.1299609797959).margin(1e-5));
}

TEST_CASE("Gaussian pulse L2 norm - Tet3D", "[unit][tet3d]") {
    double norm = runGaussianPulse(6, 3, 2, 10, 10, 10,
                                   DACOUSTICS "data/acousticsGaussian3D.h", -1);
    REQUIRE(norm == Catch::Approx(31.6577046152384).margin(1e-5));
}

TEST_CASE("Gaussian pulse L2 norm - Hex3D", "[unit][hex3d]") {
    double norm = runGaussianPulse(12, 3, 2, 10, 10, 10,
                                   DACOUSTICS "data/acousticsGaussian3D.h", -1);
    REQUIRE(norm == Catch::Approx(31.6576028812776).margin(1e-5));
}
```

`DACOUSTICS` macro (defined in `solvers/acoustics/acoustics.hpp` line 37) expands to the absolute path of `solvers/acoustics/`, giving a portable path to data files without hardcoding.

**Integration smoke tests (crash-only)**
Source: CONTEXT.md D-03; fork pattern from `libparanumal-dtu/solvers/acoustics/tests/acousticsTestMain.c` lines 96–120 (adapted from WAV-compare to no-assert)

```cpp
TEST_CASE("BC smoke - PerfRefl Tri2D", "[integration][perf_refl][tri2d]") {
    // boundary_flag=1 → bc==1 in acousticsGaussian2D.h → wall reflection
    REQUIRE_NOTHROW(runGaussianPulse(3, 2, 4, 10, 10, 10,
                                     DACOUSTICS "data/acousticsGaussian2D.h",
                                     /*boundary_flag=*/1));
}

TEST_CASE("BC smoke - FreqIndep Tri2D", "[integration][freq_indep][tri2d]") {
    // Phase 1 placeholder: runs with wall BC; norm assertion added in Phase 3
    REQUIRE_NOTHROW(runGaussianPulse(3, 2, 4, 10, 10, 10,
                                     DACOUSTICS "data/acousticsGaussian2D.h",
                                     /*boundary_flag=*/1));
}
```

Fork TEST_CASE names to adapt (from `acousticsTestMain.c` lines 96, 120):
- `"Studio with freq. indep. boundaries"` → `"BC smoke - FreqIndep Tri2D"` (upstream BOX geometry, Phase 1)
- `"Studio with perf. refl. boundaries"` → `"BC smoke - PerfRefl Tri2D"`

**Include pattern for acoustics.hpp**
Source: `solvers/acoustics/acousticsMain.cpp` line 27; `solvers/acoustics/acoustics.hpp` lines 27–91

```cpp
#include "catch_amalgamated.hpp"   // relative — vendored in same tests/ dir
#include "acoustics.hpp"           // relative — acoustics solver: acoustics_t, comm_t, platform_t, etc.
#include <cstdio>                  // popen/pclose for subprocess norm capture
#include <filesystem>              // std::filesystem::temp_directory_path
#include <string>
#include <stdexcept>
```

`acoustics.hpp` already pulls in `core.hpp`, `platform.hpp`, `mesh.hpp`, `solver.hpp`, `timeStepper.hpp`, `linAlg.hpp` and declares `using namespace libp;` (line 39) — no extra includes needed for `comm_t`, `platform_t`, `mesh_t`, `acoustics_t`.

---

### `solvers/acoustics/tests/catch_amalgamated.hpp` and `catch_amalgamated.cpp` (vendored)

**Analog:** `~/github/libparanumal-dtu/solvers/acoustics/tests/catch.hpp` (Catch2 v2 single-header)

These files are not written by hand — they are downloaded from the official release and vendored verbatim:

```bash
cd solvers/acoustics/tests/
curl -LO https://github.com/catchorg/Catch2/releases/download/v3.14.0/catch_amalgamated.hpp
curl -LO https://github.com/catchorg/Catch2/releases/download/v3.14.0/catch_amalgamated.cpp
```

The fork's `catch.hpp` is Catch2 v2 — do not copy it (CONTEXT.md D-06). The v3 amalgamated distribution replaces it.

---

### `solvers/acoustics/makefile` (config, batch) — adding `tests` and `clean-tests` targets

**Analog:** `solvers/acoustics/makefile` itself (the existing file; targets are added to it)
**Secondary analog:** `~/github/libparanumal-dtu/solvers/acoustics/makefile` lines 138–143 (fork `tests` target structure)

**Filter guard — line 65–66 (must be modified)**
Source: `solvers/acoustics/makefile` lines 65–66

Current:
```makefile
ifeq (,$(filter acousticsMain lib clean clean-libs clean-kernels \
                realclean info help test,$(MAKECMDGOALS)))
```

Replace with (add `tests clean-tests` to the list):
```makefile
ifeq (,$(filter acousticsMain lib clean clean-libs clean-kernels \
                realclean info help test tests clean-tests,$(MAKECMDGOALS)))
```

**Help message block — lines 27–63 (must be updated)**
Add `make tests` and `make clean-tests` entries to the `ACOUSTICS_HELP_MSG` define block, following the existing format:

```makefile
   make tests
   make clean-tests

...

make tests
   Build and run the Catch2 unit and integration test suite.
make clean-tests
   Clean the test binary and test object files.
```

**New object compilation rule for `tests/` directory**
Source pattern: `solvers/acoustics/makefile` lines 141–147 (existing `%.o` rule); same flags, different directory prefix

```makefile
# Test source files
TESTS_OBJS = tests/catch_amalgamated.o tests/acousticsTests.o

# Compile test object files (same flags as solver objects, plus -I tests/ for catch header)
tests/%.o: tests/%.cpp $(DEPS) | libp_libs
ifneq (,${verbose})
	$(LIBP_CXX) -o $@ -c $< $(ACOUSTICS_CXXFLAGS) -I tests/
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(LIBP_CXX) -o $@ -c $< $(ACOUSTICS_CXXFLAGS) -I tests/
endif
```

**`acousticsTests` link target and `tests` phony**
Source pattern: `solvers/acoustics/makefile` lines 124–130 (`acousticsMain` link rule); fork makefile lines 138–139

```makefile
acousticsTests: lib $(TESTS_OBJS)
ifneq (,${verbose})
	$(LIBP_LD) -o acousticsTests $(TESTS_OBJS) -L. -lacoustics $(LFLAGS)
else
	@printf "%b" "$(EXE_COLOR)Linking $(@F)$(NO_COLOR)\n";
	@$(LIBP_LD) -o acousticsTests $(TESTS_OBJS) -L. -lacoustics $(LFLAGS)
endif

tests: acousticsTests
```

Key differences from `acousticsMain` link rule:
- Dependency is `lib` (not `libp_libs`) — `lib` builds `libacoustics.a` which already depends on `libp_libs`
- Link order: `$(TESTS_OBJS) -L. -lacoustics $(LFLAGS)` — `-lacoustics` before `$(LFLAGS)` (which contains `-ltimeStepper -lmesh ...`)
- `acousticsMain.o` is replaced by `tests/catch_amalgamated.o tests/acousticsTests.o`

**`clean-tests` target**
Source pattern: `solvers/acoustics/makefile` line 151 (existing `clean` target)

```makefile
clean-tests:
	rm -f tests/*.o acousticsTests
```

**`.PHONY` update — line 110–111**
Source: `solvers/acoustics/makefile` lines 110–111

```makefile
.PHONY: all lib libp_libs clean clean-libs \
        clean-kernels realclean help info tests clean-tests
```

---

### `.github/workflows/build.yml` (config, CI) — adding two new steps

**Analog:** `.github/workflows/build.yml` itself (the existing file; new steps are appended after the existing `Test` step)

**Existing `Test` step structure (lines 27–30) — do not modify**
Source: `.github/workflows/build.yml` lines 27–30

```yaml
    - name: Test
      run: |
        make test
        bash <(curl --no-buffer -s https://codecov.io/bash) -x gcov
      env:
        LIBP_COVERAGE: 1
```

**New steps to append immediately after `Test`**
Source pattern: `.github/workflows/build.yml` lines 22–30 (Build and Test step structure); RESEARCH.md CI Step Addition pattern

```yaml
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

Rationale for two separate steps (Build vs Run): separate failure attribution — a link failure and a test assertion failure show up as distinct CI steps in the PR check. The `LIBP_COVERAGE: 1` env on the build step ensures the test objects are compiled with `--coverage -fprofile-abs-path` (from `make.top` lines 73–75), matching the earlier `libacoustics.a` build.

The `cd solvers/acoustics` is needed because `.github/workflows/build.yml` runs from `LIBP_DIR` (repo root) and the `make tests` target is in the solver subdirectory. This matches the existing `test/` directory convention where `make test` is run from the repo root (line 27: `make test`).

---

## Shared Patterns

### MPI Init/Finalize Bookend
**Source:** `solvers/acoustics/acousticsMain.cpp` lines 32, 66
**Apply to:** `acousticsTests.cpp` `main()` function

```cpp
Comm::Init(argc, argv);
// ... Catch::Session().run(argc, argv) ...
Comm::Finalize();
```

`Comm::Init` and `Comm::Finalize` are the libp wrappers around `MPI_Init`/`MPI_Finalize` (verified in `include/comm.hpp`). They must bracket the entire Catch2 session because the solver uses MPI internally.

### DACOUSTICS Macro for Portable Paths
**Source:** `solvers/acoustics/acoustics.hpp` line 37
**Apply to:** All data file path strings in `acousticsTests.cpp`

```cpp
#define DACOUSTICS LIBP_DIR"/solvers/acoustics/"
```

Usage in test source:
```cpp
DACOUSTICS "data/acousticsGaussian2D.h"   // → absolute path at compile time
DACOUSTICS "data/acousticsGaussian3D.h"
```

### Verbose-guarded printf build output
**Source:** `solvers/acoustics/makefile` lines 125–130, 133–138, 143–147
**Apply to:** All new makefile rules (`tests/%.o`, `acousticsTests` link rule)

```makefile
ifneq (,${verbose})
	$(LIBP_CXX) -o $@ ...
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(LIBP_CXX) -o $@ ...
endif
```

The `$(OBJ_COLOR)` / `$(EXE_COLOR)` / `$(NO_COLOR)` variables are defined in `make.top` lines 77–83 and are available to all sub-makefiles.

### ACOUSTICS_CXXFLAGS Reuse
**Source:** `solvers/acoustics/makefile` lines 80–99
**Apply to:** `tests/%.o` compilation rule

```makefile
ACOUSTICS_CXXFLAGS=${LIBP_CXXFLAGS} ${DEFINES} ${INCLUDES}
```

Test objects must use the same flags as solver objects so coverage, optimization, and include paths match. Add `-I tests/` to the compilation of test objects only (so `#include "catch_amalgamated.hpp"` resolves by relative path from the `tests/` directory).

### LFLAGS Reuse for Linking
**Source:** `solvers/acoustics/makefile` line 99
**Apply to:** `acousticsTests` link rule

```makefile
LFLAGS=${ACOUSTICS_CXXFLAGS} ${LIBS}
# where LIBS=-L${LIBP_LIBS_DIR} -ltimeStepper -lmesh -lparAdogs -logs -llinAlg -lcore \
#            -Wl,-rpath,$(LIBP_BLAS_DIR) -L$(LIBP_BLAS_DIR) -lopenblas \
#            -Wl,-rpath,$(OCCA_DIR)/lib -L$(OCCA_DIR)/lib -locca
```

The test binary link command prepends `-L. -lacoustics` before `$(LFLAGS)`:
```makefile
$(LIBP_LD) -o acousticsTests $(TESTS_OBJS) -L. -lacoustics $(LFLAGS)
```

---

## No Analog Found

| File | Role | Data Flow | Reason |
|---|---|---|---|
| `solvers/acoustics/tests/catch_amalgamated.cpp` | config (vendored) | — | No compiled test framework source exists anywhere in the upstream repo; file is downloaded from Catch2 v3.14.0 release |

---

## Key Pitfalls (from RESEARCH.md — planner must include mitigations in task actions)

| Pitfall | Mitigation | Source |
|---|---|---|
| `printf` stdout not captured by `std::cout.rdbuf()` | Use subprocess invoke (`popen`) or `dup2`/`pipe` fd redirect | RESEARCH.md Pitfall 1; `acousticsRun.cpp` line 63 |
| `make tests` rejected by filter guard | Add `tests clean-tests` to `$(filter ...)` on makefile line 65 | RESEARCH.md Pitfall 2; makefile lines 65–66 |
| Unit test norm with `-np 4` gives wrong reference | Run unit tests with `-np 1`; 4-rank reference is `10.1300558638317` (different value) | RESEARCH.md Pitfall 5; `testAcoustics.py` line 83 |
| `libacoustics.a` built without coverage flags | CI already runs fresh `make lib` before `make tests` — no extra action needed for Phase 1 | RESEARCH.md Pitfall 3 |

---

## Metadata

**Analog search scope:** `solvers/acoustics/`, `test/`, `.github/workflows/`, `make.top`, `~/github/libparanumal-dtu/solvers/acoustics/`
**Files read:** 10
**Pattern extraction date:** 2026-04-23
