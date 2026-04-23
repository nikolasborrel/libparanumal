# Phase 1: Test Suite - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-04-23
**Phase:** 01-test-suite
**Areas discussed:** Test assertion strategy, Catch2 version, Build integration, CI integration

---

## Test Assertion Strategy

| Option | Description | Selected |
|--------|-------------|----------|
| Norm-based only — all pass | Unit tests: Gaussian pulse L2 norm. BC integration tests: solver runs without crash. No WAV, no tinywav. All tests green from day one. | ✓ |
| Port fork structure + expected-fail shells | Bring exact fork test structure with WAV assertions; mark BC/WAV tests as [!shouldfail] until Phase 3–6 land. | |
| Skeleton only, assertions added per phase | Phase 1 just sets up Catch2 build infra, actual assertions added per phase. | |

**User's choice:** Norm-based only — all pass

**Follow-up — BC tests specifically:**

| Option | Description | Selected |
|--------|-------------|----------|
| Setup + run without crash | Tests call solver setup with FreqIndep/PerfRefl BC flags and verify it exits cleanly. Norm assertions added in Phase 3. | ✓ |
| Norm with current upstream values | Run with existing upwind flux and record those norms; not physically meaningful for BCs. | |

**User's choice:** Setup + run without crash

---

## Catch2 Version

| Option | Description | Selected |
|--------|-------------|----------|
| v2 single-header from fork | Copy catch.hpp — already proven, zero build dependency, zero API changes. | |
| v3 amalgamated (catch_amalgamated.hpp) | Current release, two-file approach (header + cpp), future-proof. | ✓ |

**User's choice:** v3 amalgamated
**Notes:** User explicitly chose v3 after reviewing options out-of-band ("Let's use latest Catch2 v3 then"). v3 is not purely header-only — `catch_amalgamated.cpp` must be compiled once and linked.

---

## Build Integration

| Option | Description | Selected |
|--------|-------------|----------|
| solvers/acoustics/tests/ + make tests | Mirror fork layout. New make tests target links against libacoustics.a. make test (Python) unchanged. | ✓ |
| top-level test/catch2/ + make test-catch2 | Place Catch2 tests alongside Python runner in test/. | |

**User's choice:** solvers/acoustics/tests/ + make tests

---

## CI Integration

| Option | Description | Selected |
|--------|-------------|----------|
| New step in build.yml after make test | Dedicated CI step runs make tests in solvers/acoustics/. Both suites visible as separate steps. | ✓ |
| Separate workflow file | Standalone .github/workflows/test-acoustics-catch2.yml. | |

**User's choice:** New step in build.yml after make test

---

## Claude's Discretion

- Test binary name
- Whether to expose `make tests` from the top-level makefile as a convenience alias
- MPI invocation style for Catch2 tests in CI

## Deferred Ideas

- WAV file comparison — Phase 6
- LR BC assertions — Phase 5
- FreqIndep/PerfRefl norm assertions — Phase 3
- Expected-fail annotations — not needed until Phase 3+
