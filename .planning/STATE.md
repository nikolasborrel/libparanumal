# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-04-22)

**Core value:** Every merged PR leaves the acoustics solver more capable and correct — correctness fixes first, then new features, all independently reviewable
**Current focus:** Phase 1 — Test Suite

## Current Position

Phase: 1 of 7 (Test Suite)
Plan: 0 of TBD in current phase
Status: Ready to plan
Last activity: 2026-04-22 — Roadmap created; all 7 phases defined with 27 requirements mapped

Progress: [░░░░░░░░░░] 0%

## Performance Metrics

**Velocity:**
- Total plans completed: 0
- Average duration: -
- Total execution time: 0 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| - | - | - | - |

**Recent Trend:**
- Last 5 plans: (none yet)
- Trend: -

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Roadmap init: ER BCs excluded — not working correctly in multi-GPU setup; not ready for upstream
- Roadmap init: Stacked PR strategy — keeps changes small and independently reviewable
- Roadmap init: Tests-first ordering — establishes correctness baseline before any code changes land
- Roadmap init: HDF5/WAV before GRF — GRF output needs compact I/O to be useful at scale
- Roadmap init: Fork C code is reference only — all new code goes into C++ `solver_t` class hierarchy

### Pending Todos

None yet.

### Blockers/Concerns

- Phase 1 tests may initially fail for features not yet landed (Phases 3–7); those tests should be annotated as expected-fail in CI until the relevant phase merges
- Phase 2 kernel fixes must be diffed carefully against `ACOUSTICS_DIFF_ANALYSIS.md` section 6 — the `@global` annotation is already present in the upstream Hex3D kernel; verify the other three surface kernels
- Phase 5 (LR BCs) introduces implicit/explicit time-stepping; interaction with the existing `timeStepper_t` hierarchy needs design review during planning

## Deferred Items

| Category | Item | Status | Deferred At |
|----------|------|--------|-------------|
| *(none)* | | | |

## Session Continuity

Last session: 2026-04-22
Stopped at: Roadmap created — ready to plan Phase 1
Resume file: None
