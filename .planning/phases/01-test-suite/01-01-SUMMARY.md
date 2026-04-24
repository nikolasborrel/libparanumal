---
phase: 01-test-suite
plan: 01
subsystem: acoustics-solver-tests
tags: [catch2, testing, makefile, macos-arm64, occa]
dependency_graph:
  requires: []
  provides:
    - acousticsTests binary (Catch2 v3.14.0, 8 test cases)
    - make tests / make clean-tests targets in solvers/acoustics/makefile
  affects:
    - occa/src/occa/internal/utils/runFunction.cpp_codegen (ARM64 fix)
    - make.top (macOS/Linux platform detection)
    - libs/timeStepper/*.cpp (VLA removal)
    - libs/core/platformDeviceConfig.cpp (macOS core count)
tech_stack:
  added:
    - Catch2 v3.14.0 (vendored amalgamated, solvers/acoustics/tests/)
  patterns:
    - popen() subprocess pattern to capture printf output from acousticsMain
    - mkstemp() for TOCTOU-safe temp .rc file creation
    - MPI-bracketed Catch2 main: Comm::Init / Session().run() / Comm::Finalize
key_files:
  created:
    - solvers/acoustics/tests/catch_amalgamated.hpp
    - solvers/acoustics/tests/catch_amalgamated.cpp
    - solvers/acoustics/tests/acousticsTests.cpp
  modified:
    - solvers/acoustics/makefile
    - make.top
    - libs/core/platformDeviceConfig.cpp
    - libs/timeStepper/timeStepperAB3.cpp
    - libs/timeStepper/timeStepperEXTBDF3.cpp
    - libs/timeStepper/timeStepperMRAB3.cpp
    - libs/timeStepper/timeStepperMRSAAB3.cpp
    - libs/timeStepper/timeStepperSAAB3.cpp
    - libs/timeStepper/timeStepperSARK4.cpp
    - libs/timeStepper/timeStepperSARK5.cpp
    - libs/timeStepper/timeStepperSSBDF3.cpp
    - occa/src/occa/internal/utils/runFunction.cpp_codegen
decisions:
  - "Use popen() subprocess to invoke acousticsMain rather than in-process linking: acousticsRun.cpp uses printf to C fd 1, not std::cout, so rdbuf redirect cannot capture it"
  - "Vendor Catch2 amalgamated into tests/ rather than using CMake FetchContent: solver uses plain Make, no CMake dependency"
  - "Fix OCCA runFunction for ARM64 by casting void(*)(...) to void(*)(void*,...) per-argc: avoids variadic ABI mismatch on AArch64 where args end up on stack not in x0/x1"
  - "Compile catch_amalgamated.cpp with -DCATCH_AMALGAMATED_CUSTOM_MAIN to suppress built-in main; acousticsTests.cpp provides the custom main with MPI bracketing"
metrics:
  duration: "~3 hours (including debugging macOS ARM64 OCCA crash)"
  completed_date: "2026-04-24"
  tasks_completed: 3
  tasks_total: 3
  files_created: 3
  files_modified: 12
---

# Phase 01 Plan 01: Catch2 Test Suite Baseline Summary

Catch2 v3.14.0 test suite wired into the acoustics solver makefile; 4 unit tests assert Gaussian pulse L2 norms within 1e-5 of reference values; 4 BC smoke tests assert no crash.

## Tasks Completed

| Task | Name | Commit | Key Files |
|------|------|--------|-----------|
| 1 | Vendor Catch2 v3.14.0 amalgamated files | 831a34a1 | tests/catch_amalgamated.{hpp,cpp} |
| 2 | Write acousticsTests.cpp | 320d63f4 | tests/acousticsTests.cpp |
| 3 | Wire makefile + macOS ARM64 fixes | 61087f62 | makefile, make.top, 10 .cpp files, occa submodule |

## Verification Results

All success criteria from PLAN.md satisfied:

- `make tests` builds `acousticsTests` binary
- `mpirun --oversubscribe -np 1 ./acousticsTests` reports **All tests passed (8 assertions in 8 test cases)**
- `./acousticsTests "[unit]"` selects 4 unit tests; all pass
- `./acousticsTests "[integration]"` selects 4 BC smoke tests; all pass
- `make clean-tests` removes test binary and object files
- Existing targets (`make acousticsMain`, `make lib`, `make clean`, `make realclean`, `make test`) unchanged

Reference norms verified:
- Tri2D: 10.130232243319 (expected 10.1302322430996, delta < 1e-5)
- Quad2D: 10.129960979796 (expected 10.1299609797959, delta < 1e-5)
- Tet3D: 31.657704615238 (expected 31.6577046152384, delta < 1e-5)
- Hex3D: 31.657602881278 (expected 31.6576028812776, delta < 1e-5)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] macOS ARM64: `-fopenmp` unsupported by Apple clang**
- **Found during:** Task 3 (first build attempt)
- **Issue:** Apple clang does not accept `-fopenmp`; requires `-Xpreprocessor -fopenmp` + Homebrew libomp
- **Fix:** Added `UNAME_S := $(shell uname -s)` platform detection block to `make.top`; macOS branch uses `-Xpreprocessor -fopenmp` with `/opt/homebrew/opt/libomp`
- **Files modified:** `make.top`
- **Commit:** 61087f62

**2. [Rule 3 - Blocking] macOS ARM64: `-mavx2` unsupported on arm64 target**
- **Found during:** Task 3
- **Issue:** `-mavx2` is an x86 SIMD flag; arm64 has NEON, not AVX2
- **Fix:** Replaced with `-ftree-vectorize -march=native` inside `ifeq (Darwin)` branch of `make.top`
- **Files modified:** `make.top`
- **Commit:** 61087f62

**3. [Rule 3 - Blocking] macOS ARM64: `<parallel/algorithm>` not found**
- **Found during:** Task 3
- **Issue:** `-DGLIBCXX_PARALLEL` is a GCC libstdc++ extension; Apple clang uses libc++ which has no such header
- **Fix:** Removed `-DGLIBCXX_PARALLEL` from macOS branch in `make.top`
- **Files modified:** `make.top`
- **Commit:** 61087f62

**4. [Rule 3 - Blocking] Apple clang: VLA-with-initializer not supported**
- **Found during:** Task 3
- **Issue:** 8 timeStepper files used `dfloat arr[runtimeVar*N] = {...}` which is a GCC extension (VLA with initializer); Apple clang rejects it as C++ forbids initialized VLAs
- **Fix:** Replaced each with `const dfloat arr[KnownFixedSize] = {...}` using the compile-time constant (Nstages, Nrk) that the runtime variable always equals
- **Files modified:** timeStepperAB3.cpp, timeStepperMRAB3.cpp, timeStepperMRSAAB3.cpp, timeStepperSSBDF3.cpp, timeStepperEXTBDF3.cpp, timeStepperSAAB3.cpp, timeStepperSARK4.cpp, timeStepperSARK5.cpp
- **Commit:** 61087f62

**5. [Rule 3 - Blocking] macOS: `lscpu` command not found**
- **Found during:** Task 3 (first acousticsMain run)
- **Issue:** `platformDeviceConfig.cpp` called `lscpu` to count CPU cores; this is a Linux-only utility
- **Fix:** Added `#if defined(__APPLE__)` block using `sysctl -n hw.physicalcpu` instead
- **Files modified:** `libs/core/platformDeviceConfig.cpp`
- **Commit:** 61087f62

**6. [Rule 3 - Blocking] macOS ARM64: OCCA `runFunction` crashes with signal 11 (segfault)**
- **Found during:** Task 3 (all 8 tests failing)
- **Issue:** OCCA's Serial backend compiles kernels to `.dylib` and invokes them through `void (*)(...)` variadic function pointer. On ARM64, the compiler emits code that puts arguments on the stack at `[sp]`/`[sp+8]` rather than in registers `x0`/`x1`. The actual kernel function (non-variadic, e.g. `innerProd2(const int&, double*)`) reads from x0/x1 and finds garbage (0x2), causing an invalid-permissions fault at address 0x2 when writing the result.
- **Root cause:** ARM64 ABI / Apple AAPCS64 places variadic call args differently from non-variadic. Calling `void(*)(...)` with concrete void* arguments uses the stack-passing path, but the compiled kernel function (non-variadic) reads from registers.
- **Fix:** Replaced each `f(args[0], ..., args[N])` in `occa/src/occa/internal/utils/runFunction.cpp_codegen` with `((void(*)(void*,...,void*))f)(args[0],...,args[N])`. The explicit non-variadic cast forces the compiler to emit register-passing calls (x0..xN per AArch64 AAPCS). Committed within the OCCA submodule.
- **Files modified:** `occa/src/occa/internal/utils/runFunction.cpp_codegen` (submodule commit 9897ecb4)
- **Commit:** 61087f62 (parent repo points to updated submodule)

**7. [Rule 3 - Blocking] Catch2 duplicate `main` symbol**
- **Found during:** Task 3 (initial link of acousticsTests)
- **Issue:** Both `catch_amalgamated.o` and `acousticsTests.o` defined `main()`; linker reported duplicate symbol
- **Fix:** Added a specific rule for `tests/catch_amalgamated.o` that compiles with `-DCATCH_AMALGAMATED_CUSTOM_MAIN`, which guards out Catch2's built-in main
- **Files modified:** `solvers/acoustics/makefile`
- **Commit:** 61087f62

## Known Stubs

None. All 8 test cases assert real behavior: 4 norm tests compare to reference values within 1e-5; 4 BC smoke tests assert no crash.

## Threat Flags

None. The test infrastructure (popen subprocess, mkstemp temp file) does not introduce new network endpoints or auth paths. The OCCA fix is inside a library call-dispatch routine with no security surface.

## Self-Check: PASSED

Files verified:
- solvers/acoustics/tests/catch_amalgamated.hpp: exists
- solvers/acoustics/tests/catch_amalgamated.cpp: exists
- solvers/acoustics/tests/acousticsTests.cpp: exists
- solvers/acoustics/makefile: contains `tests:` and `clean-tests:` targets

Commits verified:
- 831a34a1: chore(01-01): vendor Catch2 v3.14.0 amalgamated files
- 320d63f4: feat(01-01): add acousticsTests.cpp with 4 unit + 4 BC smoke tests
- 61087f62: feat(01-01): wire make tests/clean-tests and add macOS ARM64 compatibility fixes

Test binary run result: All tests passed (8 assertions in 8 test cases)
