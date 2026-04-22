# Roadmap: libParanumal Acoustics Feature Merge

## Overview

A sequence of seven self-contained PRs that merge room acoustics research features from the `libparanumal-dtu` fork into the upstream `libParanumal` C++ framework. Each PR is independently reviewable and leaves the solver more capable than before. The test suite lands first to establish a correctness baseline; subsequent phases add kernel fixes, new flux variants, the receiver system, locally reacting boundary conditions, compact I/O formats, and Gaussian Random Field initial conditions.

## Phases

**Phase Numbering:**
- Integer phases (1–7): One PR each, strictly sequential
- Decimal phases: Reserved for urgent insertions via `/gsd-insert-phase`

- [ ] **Phase 1: Test Suite** - Port Catch2 unit and integration tests from the fork as a correctness baseline
- [ ] **Phase 2: Kernel Correctness Fixes** - Apply `@global` annotations, `DT` naming, and `@barrier` placement to OKL kernels
- [ ] **Phase 3: Flux Variants** - Add `upwindBC` and `central` flux functions to all four surface kernels
- [ ] **Phase 4: Receiver System** - GPU-accelerated interpolation at arbitrary receiver positions, multi-GPU capable
- [ ] **Phase 5: Locally Reacting BCs** - Frequency-dependent wall impedance BCs with implicit/explicit time-stepping
- [ ] **Phase 6: HDF5/XDMF/WAV Output** - Compact output formats for large-scale simulation data
- [ ] **Phase 7: GRF Initial Conditions** - Gaussian Random Field ICs on sparse grids with configurable random source positions

## Phase Details

### Phase 1: Test Suite
**Goal**: Developers can run a Catch2 test suite that establishes a verifiable correctness baseline for the acoustics solver before any code changes land
**Depends on**: Nothing (first phase)
**Requirements**: TESTS-01, TESTS-02, TESTS-03, TESTS-04
**Success Criteria** (what must be TRUE):
  1. `make test` (or equivalent) in `solvers/acoustics/` compiles and links the Catch2 suite against the existing C++ libParanumal build system without modification to upstream build rules
  2. Developer can run unit tests that compare acoustic pressure field values against reference data from the fork and see pass/fail per test case on stdout
  3. Developer can run integration tests covering FreqIndep and PerfRefl BCs on 2D and 3D meshes and see per-test pass/fail results
  4. CI pipeline executes the full test suite and reports pass/fail in the PR check; tests that depend on later phases are annotated as expected-fail
**Plans**: TBD

### Phase 2: Kernel Correctness Fixes
**Goal**: All acoustics OKL kernels use correct GPU annotations and naming consistent with the rest of libParanumal v0.5.0, eliminating latent GPU memory-access bugs and naming inconsistencies inherited from the pre-v0.5.0 era
**Depends on**: Phase 1
**Requirements**: KERN-01, KERN-02, KERN-03
**Success Criteria** (what must be TRUE):
  1. All four surface kernels pass `@global` annotation review: helper function pointer arguments carry the `@global` qualifier, matching the pattern already applied in the upstream `acousticsSurfaceHex3D.okl`
  2. All four volume kernels use the `DT` name for the transpose differentiation matrix, consistent with the libParanumal v0.5.0 uniform operator naming (`D` → `DT` rename, commit `6d3ec5a3`)
  3. Volume kernels contain `@barrier("local")` in all locations required for correct GPU thread synchronisation (as audited against the upstream Hex3D fix)
  4. Phase 1 test suite passes without regression after all kernel changes are applied
**Plans**: TBD

### Phase 3: Flux Variants
**Goal**: Surface kernels expose `upwindBC` and `central` flux functions alongside the existing upwind flux, enabling physically parametric boundary conditions and making physical constants (density, speed of sound) available as kernel parameters
**Depends on**: Phase 2
**Requirements**: FLUX-01, FLUX-02, FLUX-03, FLUX-04
**Success Criteria** (what must be TRUE):
  1. Each of the four surface kernel files (Tri2D, Tet3D, Quad2D, Hex3D) defines an `upwindBC` flux function that accepts a wall velocity term `vn` and can be selected at runtime
  2. Each surface kernel file defines a `central` flux function for interior faces
  3. Physical constants `p_AcConstant`, `p_rho`, and `p_c` are defined as kernel parameters accessible within the surface kernels
  4. Existing solver runs (Gaussian pulse, all four element types) produce the same L2 norm as before the flux additions, confirming no regression to the default upwind path
**Plans**: TBD

### Phase 4: Receiver System
**Goal**: Researchers can place N microphones at arbitrary positions in the domain and record acoustic pressure time series on-device without CPU round-trips, correctly across multi-GPU configurations
**Depends on**: Phase 3
**Requirements**: RECV-01, RECV-02, RECV-03, RECV-04
**Success Criteria** (what must be TRUE):
  1. A researcher adds receiver positions to the `.rc` setup file and the solver reads them without error, with the count and positions echoed to stdout during setup
  2. The solver identifies the containing mesh element for each receiver during setup (element search runs on host); element indices are verified correct for a known mesh and receiver layout
  3. Acoustic pressure is interpolated at all receiver positions via an OKL kernel at each RK stage; recorded time series match reference values from the fork's test data to within discretisation tolerance
  4. Receiver interpolation produces correct results when the simulation is run across multiple MPI ranks / GPUs, with each rank handling only receivers in its partition
**Plans**: TBD

### Phase 5: Locally Reacting BCs
**Goal**: Researchers can configure frequency-dependent, locally reacting wall impedance boundary conditions via the setup file; the solver handles the implicit/explicit time-stepping required for numerical stability and passes integration tests on canonical room geometries
**Depends on**: Phase 4
**Requirements**: BC-01, BC-02, BC-03, BC-04
**Success Criteria** (what must be TRUE):
  1. A researcher sets `[BC TYPE] LR` and provides impedance coefficients in the `.rc` file; the solver starts without error and echoes the configured BC parameters during setup
  2. LR BC time-stepping is numerically stable for the studio room and cylindrical mesh integration test cases at the polynomial degrees tested in the fork
  3. Integration tests for LR BCs on studio room and cylindrical mesh configurations pass (pressure time series match fork reference data to within tolerance)
  4. Setting constant (frequency-independent) LR coefficients reproduces the FreqIndep BC behaviour, confirmed by matching the FreqIndep integration test reference output
**Plans**: TBD

### Phase 6: HDF5/XDMF/WAV Output
**Goal**: Simulation output can be written as compact HDF5 wave field snapshots, XDMF metadata for visualisation at actual Gauss-Lobatto nodes, and WAV audio files for receiver time series — all selectable via the setup file, reducing output size by approximately 25x versus VTU
**Depends on**: Phase 5
**Requirements**: IO-01, IO-02, IO-03, IO-04
**Success Criteria** (what must be TRUE):
  1. A simulation configured with `[OUTPUT FORMAT] HDF5` writes a single `.h5` file containing all time steps as float32 data; file size is substantially smaller than the equivalent VTU output for the same run
  2. An XDMF `.xdmf` metadata file is written alongside the HDF5 file; opening it in ParaView displays the wave field at actual Gauss-Lobatto node positions without visible oversampling artifacts
  3. Receiver time series are written as `.wav` audio files via `tinywav`; files are readable by standard audio software and contain the expected number of samples
  4. Setting `[OUTPUT FORMAT] VTU` in the `.rc` file preserves existing VTU behaviour with no regression; HDF5/XDMF/WAV output is additive, not a replacement
**Plans**: TBD

### Phase 7: GRF Initial Conditions
**Goal**: Researchers can initialise the acoustic wave field with a Gaussian Random Field on a sparse uniform grid, with source positions optionally sampled at random from an input mesh — enabling stochastic room acoustics simulations and ML dataset generation
**Depends on**: Phase 6
**Requirements**: IC-01, IC-02, IC-03, IC-04
**Success Criteria** (what must be TRUE):
  1. Setting `[IC TYPE] GRF` in the `.rc` file initialises the wave field with a Gaussian Random Field; the solver starts without error and the initial pressure field is non-trivially distributed across the domain
  2. The GRF sparse grid resolution `ppw` (points per wavelength) is read from the setup file; changing `ppw` visibly changes the spatial frequency content of the initial condition
  3. Providing an input mesh file for source positions causes source locations to be drawn from that mesh; two runs with different random seeds produce different but statistically similar initial conditions
  4. GRF values are correctly interpolated from the sparse IC grid to the simulation quadrature points during setup, confirmed by comparing the resulting `o_q` to a reference snapshot from the fork's test data
**Plans**: TBD

## Progress

**Execution Order:** 1 → 2 → 3 → 4 → 5 → 6 → 7

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Test Suite | 0/TBD | Not started | - |
| 2. Kernel Correctness Fixes | 0/TBD | Not started | - |
| 3. Flux Variants | 0/TBD | Not started | - |
| 4. Receiver System | 0/TBD | Not started | - |
| 5. Locally Reacting BCs | 0/TBD | Not started | - |
| 6. HDF5/XDMF/WAV Output | 0/TBD | Not started | - |
| 7. GRF Initial Conditions | 0/TBD | Not started | - |

---
*Roadmap created: 2026-04-22*
*Last updated: 2026-04-22 after initial creation*
