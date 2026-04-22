# libParanumal Acoustics Feature Merge

## What This Is

A structured merge of research-specific acoustics features from the `libparanumal-dtu` fork into the upstream `libParanumal` DG solver library. The work ports physically realistic boundary conditions, receiver networks, stochastic initial conditions, and efficient I/O formats into the upstream C++ framework via a sequence of focused, self-contained PRs — each reviewed and merged by the upstream maintainers before the next begins.

## Core Value

Every merged PR leaves the acoustics solver more capable and correct than before — correctness fixes first, then new features, all compatible with the upstream C++ framework and independently reviewable.

## Requirements

### Validated

- ✓ Working acoustics solver with standard upwind flux — existing
- ✓ C++ OOP framework (`solver_t` class hierarchy, settings classes, OGS) — existing
- ✓ GPU kernel dispatch via OCCA — existing
- ✓ VTU output and convergence reporting — existing
- ✓ Multi-GPU via MPI halo exchange — existing

### Active

- [ ] Test suite ported from fork (Catch2 unit + integration tests) — baseline before any changes
- [ ] OKL kernel correctness fixes (`@global` annotations, `DT` naming, `@barrier` placement)
- [ ] `upwindBC` and `central` flux variants in surface kernels (enables physically parametric BCs)
- [ ] Receiver system: GPU-accelerated interpolation at arbitrary microphone positions
- [ ] Locally Reacting (LR) boundary conditions
- [ ] HDF5/XDMF/WAV output for compact, large-scale simulation data
- [ ] GRF initial conditions: Gaussian Random Fields for stochastic room acoustics

### Out of Scope

- Extended Reaction (ER) BCs — not working correctly in multi-GPU setup; excluded from merge
- ML dataset generation pipelines (`simulationSetups/deeponet/`) — fork-specific tooling, not upstream material
- C procedural port — target is C++ OOP only; fork C code is reference only

## Context

The fork (`libparanumal-dtu`) diverged from `libParanumal` in mid-2019. The original continued evolving its framework (C++ OOP, settings classes, v0.5.0 release in 2022), while the fork added room acoustics research features in C (procedural). Merging requires adapting fork code to the C++ class hierarchy.

Key technical context:
- Fork code is C procedural; target is C++ OOP (`acoustics_t` inheriting `solver_t`)
- OKL kernels are the most directly portable artifacts — shared files have diverged in annotations and naming
- Each PR targets `main` and must be independently reviewable without depending on later PRs
- Full file-level divergence analysis: `.claude/ACOUSTICS_DIFF_ANALYSIS.md`
- Codebase map: `.planning/codebase/` (ARCHITECTURE, CONCERNS, CONVENTIONS, STACK, STRUCTURE, TESTING, INTEGRATIONS)

## Constraints

- **Framework**: Must target libParanumal v0.5.0 C++ API (settings classes, OGS, platform abstraction, `solver_t` hierarchy)
- **PR strategy**: Each phase = one stacked PR; changes must be self-contained and not break existing solvers
- **Compatibility**: Must not break existing tests or other solvers in libParanumal
- **External deps**: HDF5 linkage required for output phase; `tinywav` is a single-file library (easy to add)

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| ER BCs excluded | Not working correctly in multi-GPU setup; not ready for upstream | — Pending |
| Stacked PR strategy | Keeps changes small and independently reviewable; reduces upstream review burden | — Pending |
| Tests-first ordering | Establishes correctness baseline before any code changes land | — Pending |
| HDF5/WAV before GRF | GRF output needs compact I/O to be useful at scale | — Pending |

## Evolution

This document evolves at phase transitions and milestone boundaries.

**After each phase transition** (via `/gsd-transition`):
1. Requirements invalidated? → Move to Out of Scope with reason
2. Requirements validated? → Move to Validated with phase reference
3. New requirements emerged? → Add to Active
4. Decisions to log? → Add to Key Decisions
5. "What This Is" still accurate? → Update if drifted

**After each milestone** (via `/gsd-complete-milestone`):
1. Full review of all sections
2. Core Value check — still the right priority?
3. Audit Out of Scope — reasons still valid?
4. Update Context with current state

---
*Last updated: 2026-04-22 after initialization*