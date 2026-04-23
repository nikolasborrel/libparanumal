---
phase: 1
slug: test-suite
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-04-23
---

# Phase 1 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Catch2 v3.14.0 (amalgamated) |
| **Config file** | none — standalone binary |
| **Quick run command** | `mpirun --oversubscribe -np 1 ./acousticsTests "[unit]"` |
| **Full suite command** | `mpirun --oversubscribe -np 1 ./acousticsTests` |
| **Estimated runtime** | ~60 seconds (4 solver runs × ~15s each) |

---

## Sampling Rate

- **After every task commit:** Run `cd solvers/acoustics && make tests && mpirun --oversubscribe -np 1 ./acousticsTests "[unit]"`
- **After every plan wave:** Run `cd solvers/acoustics && mpirun --oversubscribe -np 1 ./acousticsTests` (full suite)
- **Before `/gsd-verify-work`:** Full suite must be green + `make test` (Python runner) unbroken
- **Max feedback latency:** ~60 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 1-01-01 | 01 | 0 | TESTS-04 | — | N/A | build | `cd solvers/acoustics && make tests` (zero exit code) | ❌ W0 | ⬜ pending |
| 1-01-02 | 01 | 1 | TESTS-01 | — | N/A | unit | `mpirun --oversubscribe -np 1 ./acousticsTests "[unit]"` | ❌ W0 | ⬜ pending |
| 1-01-03 | 01 | 1 | TESTS-02 | — | N/A | integration | `mpirun --oversubscribe -np 1 ./acousticsTests "[integration]"` | ❌ W0 | ⬜ pending |
| 1-01-04 | 01 | 2 | TESTS-03 | — | N/A | CI smoke | CI pass/fail per test case in PR check | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `solvers/acoustics/tests/catch_amalgamated.hpp` — vendored Catch2 v3.14.0 header (download from GitHub release)
- [ ] `solvers/acoustics/tests/catch_amalgamated.cpp` — vendored Catch2 v3.14.0 source (compiled once)
- [ ] `solvers/acoustics/tests/acousticsTests.cpp` — all test cases (unit + integration); covers TESTS-01, TESTS-02
- [ ] Updated `solvers/acoustics/makefile` — `tests` and `clean-tests` targets; covers TESTS-04
- [ ] Updated `.github/workflows/build.yml` — new Catch2 build + run CI steps; covers TESTS-03

*All phase test infrastructure is new — Wave 0 creates it from scratch.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Catch2 per-test pass/fail output readable in CI log | TESTS-03 | GitHub Actions log formatting requires visual inspection | After CI run, open the "Run Catch2 Tests" step log and confirm individual test names with ✓/✗ are visible |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 60s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
