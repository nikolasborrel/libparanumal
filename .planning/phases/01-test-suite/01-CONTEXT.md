# Phase 1: Test Suite - Context

**Gathered:** 2026-04-23
**Status:** Ready for planning

<domain>
## Phase Boundary

Port Catch2 unit and integration tests from the `libparanumal-dtu` fork into the upstream `libParanumal` C++ build system, establishing a verifiable correctness baseline before any code changes land. The PR adds a new test binary and CI step — it does not modify existing solver code, build rules, or the existing Python test runner.

</domain>

<decisions>
## Implementation Decisions

### Test Assertions
- **D-01:** All Phase 1 tests must pass green from day one — no expected-fail annotations in the Phase 1 PR.
- **D-02:** Unit tests (TESTS-01) assert the mass-matrix L2 norm of the acoustic pressure field for Gaussian pulse runs on all four element types (Tri2D, Quad2D, Tet3D, Hex3D), using the same reference norm values already established in the Python test runner (`test/testAcoustics.py`).
- **D-03:** Integration tests (TESTS-02) for FreqIndep and PerfRefl BCs assert only that the solver setup runs without crashing (no norm assertion). Numerical norm assertions for these BC types are added in Phase 3 when the flux variants (`upwindBC`, `central`) land.
- **D-04:** No WAV file comparison in Phase 1 — `tinywav` is not a Phase 1 dependency. WAV-based assertions are introduced in Phase 6 when HDF5/WAV output is implemented.

### Test Framework
- **D-05:** Catch2 v3 — use the amalgamated single-file distribution (`catch_amalgamated.hpp` + `catch_amalgamated.cpp`). Fetch from the official Catch2 v3 GitHub release. The `.cpp` file is compiled once as part of the test build; it is not header-only.
- **D-06:** Do not copy `catch.hpp` from the fork — the fork uses Catch2 v2, which is in maintenance-only mode.

### Build Integration
- **D-07:** Test files live in `solvers/acoustics/tests/`. A new `make tests` target is added to `solvers/acoustics/makefile` that builds the test binary by linking against `libacoustics.a`. The existing `make test` target (Python runner) is untouched.
- **D-08:** The test binary links against the same `libacoustics.a` produced by `make lib` — no duplication of object files. `make tests` depends on `make lib`.
- **D-09:** Catch2 amalgamated source files are vendored into `solvers/acoustics/tests/` (not fetched at build time) so the build works offline and on HPC systems without internet access.

### CI Integration
- **D-10:** A new step is added to `.github/workflows/build.yml` after the existing `make test` step: build and run the Catch2 suite via `make tests` in `solvers/acoustics/`. Both the Python runner and the Catch2 suite appear as distinct steps in the PR check.
- **D-11:** The existing `make test` CI step is unchanged — the Python runner continues to run and serves as the baseline regression check.

### Claude's Discretion
- Test binary name (`acousticsTests` or similar) — follow the fork's convention if one exists, otherwise choose a clear name.
- Whether to expose `make tests` from the top-level `makefile` as well (convenience alias) — planner can decide.
- MPI invocation style for Catch2 tests in CI (single-rank vs. multi-rank) — match the existing Python runner's convention (`mpirun --oversubscribe -np 1` for unit, `-np 4` for MPI integration test).

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Existing Test Infrastructure
- `test/testAcoustics.py` — Reference norm values for all four element types; CI invocation pattern; Python test harness structure
- `.github/workflows/build.yml` — CI pipeline steps; where to add the new Catch2 step; existing `make test` invocation
- `solvers/acoustics/makefile` — Existing `make test`, `make lib`, `make acousticsMain` targets; how to add `make tests` without breaking upstream rules

### Fork Test Reference
- `.claude/ACOUSTICS_DIFF_ANALYSIS.md` — Full divergence analysis; fork test structure overview
- Fork path `~/github/libparanumal-dtu/solvers/acoustics/tests/acousticsTestMain.c` — Fork test cases (TEST_CASE names, BC configurations, setup files used); adapt to C++ API
- Fork path `~/github/libparanumal-dtu/solvers/acoustics/tests/data/ref/` — Reference WAV files and setup configurations (context only — not used in Phase 1 assertions)

### Build System
- `make.top` — Compiler flags, library paths, `LIBP_*` variables used by all sub-makefiles; new `make tests` must include this
- `.planning/codebase/STACK.md` — Technology stack overview; Catch2 v3 fit within GNU Make build system
- `.planning/codebase/TESTING.md` — Full analysis of current test patterns, gaps, and what is/isn't tested

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `libacoustics.a` (from `make lib`): Test binary links against this — no need to compile solver objects separately
- `solvers/acoustics/acoustics.hpp`: C++ class header; tests instantiate `acoustics_t` through this interface
- `solvers/acoustics/data/acousticsGaussian2D.h`, `acousticsGaussian3D.h`: IC/BC macros used by existing tests; Catch2 tests use the same data files

### Established Patterns
- Python test harness writes a `.rc` file to a temp directory, invokes the binary via `mpirun`, and parses stdout — Catch2 tests follow a similar pattern but call the C++ API directly rather than launching a subprocess
- All existing tests use BOX mesh generator, DOPRI5 integrator, CFL 1.0, t_final 1.0 — Catch2 unit tests use the same defaults for consistency with reference norms
- `LIBP_TEST_DIR` make variable points to `test/` — new `solvers/acoustics/tests/` dir is independent of this

### Integration Points
- `solvers/acoustics/makefile`: Add `make tests` and `make clean-tests` targets; link against `libacoustics.a`
- `.github/workflows/build.yml`: Add step after existing `make test` step
- The C++ solver entry is `acousticsRun()` in `src/acousticsRun.cpp` — unit tests that check norm values will call this (or the equivalent `acoustics_t` run method) and capture the printed norm

</code_context>

<specifics>
## Specific Ideas

- Catch2 v3 amalgamated files should be fetched from the official GitHub release (not built from source) and vendored into `solvers/acoustics/tests/` — this keeps the build self-contained for HPC environments
- The fork's test case structure (`TEST_CASE("Studio with freq. indep. boundaries", "[studio][freq_indep]")`) can be adapted for Phase 1 BC smoke tests, minus the WAV comparison
- Phase 1 unit test norm reference values come directly from `test/testAcoustics.py` (e.g., Tri2D: `10.1302322430996`, Quad2D: `10.1299609797959`, Tet3D: `31.6577046152384`, Hex3D: `31.6576028812776`)

</specifics>

<deferred>
## Deferred Ideas

- WAV file comparison in tests — deferred to Phase 6 when `tinywav` and WAV output are implemented
- LR BC integration test assertions — deferred to Phase 5
- Norm assertions for FreqIndep/PerfRefl BC tests — deferred to Phase 3 when flux variants land
- Expected-fail test annotations (`[!shouldfail]`) — not needed in Phase 1; relevant for later phases

</deferred>

---

*Phase: 01-test-suite*
*Context gathered: 2026-04-23*
