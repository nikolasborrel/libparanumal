---
name: Acoustics merge project
description: Ongoing work to merge libparanumal-dtu fork changes into libParanumal via stacked PRs
type: project
---

Merging research-specific acoustics features from the `libparanumal-dtu` fork into `libParanumal` via multiple stacked PRs on branch `nikolasborrel/acoustics_tests`.

**Why:** The fork added physically realistic boundary conditions (ER, LR), receiver networks, stochastic ICs (GRF), and efficient I/O (HDF5/XDMF/WAV) for room acoustics simulation. The original repo has since modernized its framework (C++ OOP, settings classes, `DT` operator naming, `@global` annotations).

**How to apply:** Detailed file-level divergence analysis is in `.claude/ACOUSTICS_DIFF_ANALYSIS.md`. Prioritize OKL kernel fixes first (shared code, correctness), then new features in order of integration effort.

Planned PR order (rough):
1. Unit and integration tests (ported from fork) — baseline to verify correctness throughout the merge
2. OKL kernel correctness fixes (`@global` annotations, `DT` naming, `@barrier`)
3. `upwindBC` / `central` flux variants in surface kernels
4. Receiver system
5. Locally Reacting (LR) BCs
6. GRF initial conditions
7. HDF5/XDMF/WAV output

Excluded: Extended Reaction (ER) BCs — research code, not working as expected, will not be merged.
