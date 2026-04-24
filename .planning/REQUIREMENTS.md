# Requirements: libParanumal Acoustics Feature Merge

**Defined:** 2026-04-22
**Core Value:** Every merged PR leaves the acoustics solver more capable and correct — correctness fixes first, then new features, all independently reviewable.

## v1 Requirements

### Tests

- [x] **TESTS-01
**: Developer can run Catch2 unit tests that verify acoustic pressure field values against reference data from the fork
- [x] **TESTS-02
**: Developer can run Catch2 integration tests covering FreqIndep and PerfRefl BCs on 2D and 3D meshes
- [ ] **TESTS-03**: CI pipeline executes the test suite and reports pass/fail per test case
- [x] **TESTS-04
**: Test infrastructure compiles and links against the existing C++ libParanumal build system

### Kernel Fixes

- [ ] **KERN-01**: Acoustic surface kernels use `@global` pointer annotations on all helper function arguments for correct GPU memory access
- [ ] **KERN-02**: Volume kernels use `DT` (transpose differentiation matrix) naming consistently across all element types (Tet3D, Hex3D, Tri2D, Quad2D)
- [ ] **KERN-03**: Volume kernels have `@barrier("local")` placed correctly for GPU thread synchronization

### Flux Variants

- [ ] **FLUX-01**: Surface kernels provide an `upwindBC` flux function with wall velocity term `vn` for all element types
- [ ] **FLUX-02**: Surface kernels provide a `central` flux function for interior faces for all element types
- [ ] **FLUX-03**: Physical constants `p_AcConstant`, `p_rho`, `p_c` (density, speed of sound) are accessible as kernel parameters in surface kernels
- [ ] **FLUX-04**: Existing upwind flux behavior is preserved (no regression for existing setups)

### Receiver System

- [ ] **RECV-01**: Researcher can specify N arbitrary receiver positions in the setup `.rc` file
- [ ] **RECV-02**: Solver finds the mesh element containing each receiver position during setup
- [ ] **RECV-03**: Acoustic pressure is interpolated at all receiver positions on-device via an OKL kernel at each RK stage
- [ ] **RECV-04**: Receiver system functions correctly in multi-GPU / multi-MPI-rank configurations

### Boundary Conditions

- [ ] **BC-01**: Researcher can configure Locally Reacting (LR) frequency-dependent wall impedance BCs via setup file
- [ ] **BC-02**: LR BCs use implicit/explicit time-stepping for numerical stability
- [ ] **BC-03**: LR BCs pass integration tests for studio room and cylindrical mesh configurations
- [ ] **BC-04**: FreqIndep (frequency-independent) wall impedance BCs work as a degenerate case of LR with constant coefficients

### I/O Formats

- [ ] **IO-01**: Simulation writes wave field snapshots in compact HDF5 format (float32, all time steps in one dataset)
- [ ] **IO-02**: Simulation writes XDMF metadata files enabling visualization of HDF5 data at actual Gauss-Lobatto nodes
- [ ] **IO-03**: Receiver time series are written as WAV audio files via `tinywav`
- [ ] **IO-04**: Output format is selectable via setup `.rc` file parameter

### Initial Conditions

- [ ] **IC-01**: Researcher can initialize the wave field with a Gaussian Random Field (GRF) on a sparse uniform grid
- [ ] **IC-02**: GRF sparse mesh resolution (`ppw`) is configurable as a setup parameter
- [ ] **IC-03**: Source positions can be randomly sampled from an input mesh file
- [ ] **IC-04**: GRF is interpolated from the sparse IC grid to simulation quadrature points during setup

## v2 Requirements

(None — scope is fully defined for this milestone)

## Out of Scope

| Feature | Reason |
|---------|--------|
| Extended Reaction (ER) BCs | Not working correctly in multi-GPU setup; not ready for upstream |
| ML dataset generation pipelines (`simulationSetups/deeponet/`) | Fork-specific tooling; not upstream material |
| FreqIndep BCs as a separate PR | Implicit in FLUX-01/02 (upwindBC enables wall velocity term) and BC-01–03 (LR as the general case); no standalone implementation needed |
| Literal C procedural port | Fork `.c` files are reference only; all upstream code is adapted into the C++ `solver_t` class hierarchy |
| VTU output removal | Fork removed VTU; upstream preserves it for backwards compatibility — HDF5/XDMF are additions, not replacements |
| Curvilinear mesh support | Not mature in the general case; lacks documentation |
| macOS build fixes | Not documented or generalised enough for upstream |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| TESTS-01 | Phase 1 | Pending |
| TESTS-02 | Phase 1 | Pending |
| TESTS-03 | Phase 1 | Pending |
| TESTS-04 | Phase 1 | Pending |
| KERN-01 | Phase 2 | Pending |
| KERN-02 | Phase 2 | Pending |
| KERN-03 | Phase 2 | Pending |
| FLUX-01 | Phase 3 | Pending |
| FLUX-02 | Phase 3 | Pending |
| FLUX-03 | Phase 3 | Pending |
| FLUX-04 | Phase 3 | Pending |
| RECV-01 | Phase 4 | Pending |
| RECV-02 | Phase 4 | Pending |
| RECV-03 | Phase 4 | Pending |
| RECV-04 | Phase 4 | Pending |
| BC-01 | Phase 5 | Pending |
| BC-02 | Phase 5 | Pending |
| BC-03 | Phase 5 | Pending |
| BC-04 | Phase 5 | Pending |
| IO-01 | Phase 6 | Pending |
| IO-02 | Phase 6 | Pending |
| IO-03 | Phase 6 | Pending |
| IO-04 | Phase 6 | Pending |
| IC-01 | Phase 7 | Pending |
| IC-02 | Phase 7 | Pending |
| IC-03 | Phase 7 | Pending |
| IC-04 | Phase 7 | Pending |

**Coverage:**
- v1 requirements: 27 total
- Mapped to phases: 27
- Unmapped: 0 ✓

---
*Requirements defined: 2026-04-22*
*Last updated: 2026-04-22 after initial definition*
