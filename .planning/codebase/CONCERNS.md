# Codebase Concerns — solvers/acoustics/

**Analysis Date:** 2026-04-22
**Scope:** `solvers/acoustics/` within the libParanumal HPC finite-element library.

---

## Numerical Correctness Concerns

### Inconsistent upwind flux formulas across element types

**Files:**
- `solvers/acoustics/okl/acousticsSurfaceTri2D.okl` lines 43–46
- `solvers/acoustics/okl/acousticsSurfaceTet3D.okl` lines 47–52
- `solvers/acoustics/okl/acousticsSurfaceQuad2D.okl` lines 43–47
- `solvers/acoustics/okl/acousticsSurfaceHex3D.okl` lines 47–53

**Issue:** The `upwind()` flux helper uses two distinct algebraic conventions depending on the element type. Tri2D and Tet3D compute the *correction form* (numerical flux minus the interior trace flux):

```c
// Tri2D / Tet3D (correction form — matches DG weak form + LIFT*(F*-F^-))
*rflux = p_half*((ndotUP - ndotUM) - (rP - rM));
```

Quad2D and Hex3D compute the *total numerical flux form*:

```c
// Quad2D / Hex3D (total flux form — matches SBP strong form)
*rflux = p_half*(ndotUM + ndotUP - (rP - rM));
```

These conventions are intentionally different because the Tri/Tet kernels use the standard DG weak formulation while the Quad/Hex kernels use a summation-by-parts (SBP) strong form. The two forms are algebraically consistent with their respective volume operators. However, there is no comment in the code explaining this distinction, making it easy to mistake one form for a bug in the other. The test suite passes with different norms for Tri (10.1302) and Quad (10.1299), which is expected discretisation-error variation, not a sign of correctness.

**Risk:** Future contributors editing one surface kernel without awareness of the matching volume-kernel convention will silently corrupt the scheme.

**Fix approach:** Add a comment in each surface kernel identifying which DG formulation it implements (weak/strong) and cross-referencing the matching volume kernel.

---

### Incomplete Neumann boundary condition — Tri2D only

**File:** `solvers/acoustics/okl/acousticsSurfaceTri2D.okl` line 118

**Issue:** An explicit comment acknowledges that the Neumann BC is not implemented:

```c
//should also add the Neumann BC here, but need uxM, uyM, vxM, abd vyM somehow
```

This comment appears only in the Tri2D surface kernel. The other three surface kernels (Tet3D, Quad2D, Hex3D) have no equivalent comment and no Neumann implementation. The boundary-condition macros in the data files (`acousticsGaussian2D.h`, `acousticsGaussian3D.h`) only define conditions for `bc==1` (wall/reflect) and `bc==2` (outflow). Neumann conditions require gradient values that are not threaded into the surface kernels.

**Impact:** Outflow boundaries use a simple `*(rB) = -rM` anti-reflection hack rather than a true non-reflecting or Neumann condition. This produces non-physical reflections for any problem where the wave reaches the domain boundary.

**Fix approach:** For the simplex (Tri/Tet) kernels, compute and pass face-normal derivatives through an auxiliary gradient solve or embed them in the state vector. For Quad/Hex, the SBP structure may simplify this. Alternatively, implement a proper absorbing boundary via the Bayliss-Turkel or `*(rB) = rM` characteristic approach.

---

### `MaxWaveSpeed()` always returns 1.0 — no physical speed parameterisation

**File:** `solvers/acoustics/src/acousticsStep.cpp` lines 29–33

```cpp
dfloat acoustics_t::MaxWaveSpeed(){
  //wavespeed is constant 1 everywhere
  const dfloat vmax = 1.0;
  return vmax;
}
```

The solver assumes a unit wave speed everywhere in the domain. There is no mechanism to set a spatially varying or problem-specific speed of sound. The `data/*.h` files and initial condition macros do not expose a sound-speed parameter.

**Impact:** The solver cannot model problems where the sound speed differs from 1 (e.g., heterogeneous media, different physical units). Extending to non-unit wave speed requires changes to both the volume and surface kernels.

**Fix approach:** Add a sound-speed field `c` as a mesh-level array (analogous to BNS's viscosity parameter), and thread it into both the volume flux computations and the `MaxWaveSpeed()` kernel.

---

### Undocumented `JW` weighting in Hex3D volume kernel

**File:** `solvers/acoustics/okl/acousticsVolumeHex3D.okl` line 70

```c
// (1/J) \hat{div} (G*[F;G])
// questionable: why JW
```

The author left an unresolved question about why the Jacobian-weight `JW` is used in the SBP flux computation. This is correct for the SBP skew-symmetric entropy-stable formulation (where fluxes are pre-multiplied by the quadrature weight before differentiation), but the comment signals uncertainty. The matching `invJW` multiplication at the output stage (line 115) is the correct inverse operation. Both Quad2D and Hex3D share this pattern.

**Risk:** If the code is refactored without understanding the SBP context, removing `JW` and `invJW` will break the scheme for these element types.

**Fix approach:** Replace the `// questionable: why JW` comment with a short explanation of the SBP strong form and a reference to the relevant section of the paper/documentation.

---

### `time` variable implicitly captured in `acousticsSurfaceHex3D.okl` device function

**File:** `solvers/acoustics/okl/acousticsSurfaceHex3D.okl` lines 57–115

The `surfaceTerms()` device function calls `acousticsDirichletConditions3D(bc, time, ...)` (line 102) but `time` is not a parameter of `surfaceTerms()`. The function signature (lines 57–71) has no `time` argument. OCCA device functions can access kernel arguments by implicit capture, but this is an undocumented reliance on OCCA scoping semantics.

**Impact:** If OCCA's device function scoping rules change, or if `surfaceTerms()` is moved/reused outside the kernel scope, it will fail silently with a wrong (likely zero) `time` value. The equivalent Quad2D `surfaceTerms()` helper has the same pattern (line 93). The simplex kernels (Tri2D, Tet3D) do not use a helper function, so they avoid this issue.

**Fix approach:** Pass `time` explicitly as an argument to `surfaceTerms()` in both `acousticsSurfaceHex3D.okl` and `acousticsSurfaceQuad2D.okl`, matching the explicit parameter passing used in the Tri2D/Tet3D kernels.

---

## Technical Debt

### Three dead volume kernel variants in `acousticsVolumeTet3D.okl`

**File:** `solvers/acoustics/okl/acousticsVolumeTet3D.okl` lines 29–329

The file contains four `@kernel` definitions:
- `acousticsVolumeTet3D_v0` (lines 29–130): F/G/H shared array approach
- `acousticsVolumeTet3D_v1` (lines 134–207): separate rho/u/v/w shared arrays, `p_Nvol=1`
- `acousticsVolumeTet3D_v2` (lines 212–329): batched elements, hardcoded `#define p_Nvol 2`
- `acousticsVolumeTet3D` (lines 333–457): final variant, hardcoded `#define p_Nvol 1`, `#define p_NblockV 4`

The setup code in `acousticsSetup.cpp` (line 122) selects only `"acousticsVolumeTet3D"` by name. The three `_v0`, `_v1`, `_v2` variants are never selected by any host code path; they are dead code from a performance-tuning exploration.

**Impact:** OCCA JIT-compiles all kernels in the file at startup, increasing compilation time. The dead variants add maintenance noise and confusion.

**Fix approach:** Move the dead variants to a separate `acousticsVolumeTet3D_archive.okl` or delete them; retain only the production variant.

---

### Hardcoded tuning parameters inside OKL kernel bodies

**File:** `solvers/acoustics/okl/acousticsVolumeTet3D.okl` lines 339–340

```c
#define p_Nvol 1
#define p_NblockV 4
```

These `#define` statements sit inside the kernel function body rather than being defined by host-side `kernelInfo["defines/..."]` properties. This prevents tuning without modifying the OKL file. All other block-size parameters (`p_NblockV`, `p_NblockS`) are correctly injected from host code in `acousticsSetup.cpp` lines 96–99.

**Fix approach:** Remove the in-kernel `#define` statements and inject these values from `acousticsSetup.cpp` via `kernelInfo["defines/" "p_Nvol"]` and `kernelInfo["defines/" "p_NblockV"]`, consistent with the rest of the codebase.

---

### `CFL NUMBER` is set but never reported in settings output

**File:** `solvers/acoustics/src/acousticsSettings.cpp` lines 42–44 and 67–78

`CFL NUMBER` is registered as a setting (line 42) but `acousticsSettings_t::report()` does not call `reportSetting("CFL NUMBER")`. The CFL number directly affects solution stability but is invisible in the console output.

**Fix approach:** Add `reportSetting("CFL NUMBER");` inside the `report()` function.

---

## Missing Features Compared to Sibling Solvers

### No PML / absorbing layer support

**Comparison:** `solvers/bns/` has a full perfectly-matched layer (PML) infrastructure including:
- `bnsPmlSetup.cpp`, `bnsRelaxation*.okl` kernels, `rhsf_pml()` and `rhsf_MR_pml()` virtual overrides
- Settings: `PML PROFILE ORDER`, `PML SIGMAX MAX`, `PML SIGMAY MAX`, `PML SIGMAZ MAX`, `PML INTEGRATION`

The acoustics solver has no absorbing boundary or PML capability. The base `solver_t` class declares `rhsf_pml()` as a virtual method (`include/solver.hpp` lines 90–93), but `acoustics_t` never overrides it. Outflow boundaries use the reflection-suppression hack (`*(rB) = -rM`) rather than a proper absorbing condition.

**Impact:** Simulating unbounded domains (free-field acoustic propagation) requires a workaround with a large domain, wasting compute resources.

---

### No GPU-computed `MaxWaveSpeed` kernel

**Comparison:** `solvers/advection/` has per-element `maxWaveSpeedKernel` for all four element types (`advectionMaxWaveSpeed*.okl`). The wave speed is computed on the GPU using the actual solution state, enabling adaptive time stepping with DOPRI5.

Acoustics has a trivial CPU-side `MaxWaveSpeed()` that returns `1.0` unconditionally and takes no arguments. The method signature diverges from the advection solver's `MaxWaveSpeed(deviceMemory<dfloat>& o_Q, const dfloat T)` interface.

**Impact:** When DOPRI5 is used (the default time integrator per `acousticsSettings.cpp` line 39), adaptive step control relies on the error estimate from DOPRI5 itself, not wave-speed-based CFL monitoring. This is functional but not consistent with other solvers.

---

### No convergence study or error norm output

The `Report()` function (`acousticsReport.cpp`) and `Run()` (`acousticsRun.cpp`) only output the mass-matrix norm of `q`. There is no comparison against an analytical solution, no L2 error norm, and no convergence rate measurement.

The Gaussian initial condition in `acousticsGaussian2D.h` and `acousticsGaussian3D.h` does not have a closed-form exact solution at later times, which would be needed for convergence testing. The test suite (`test/testAcoustics.py`) only compares a reference norm at a fixed final time, not convergence order.

**Impact:** It is impossible to verify p-convergence or detect scheme regressions without adding an analytically solvable test case (e.g., a plane wave or mode with known exact solution).

---

## Performance Concerns

### Hex3D surface kernel uses direct global writes with multiple barrier points

**File:** `solvers/acoustics/okl/acousticsSurfaceHex3D.okl`

The Hex3D surface kernel calls `surfaceTerms()` which writes directly to global memory `rhsq` (not to shared memory). Three separate `@barrier()` synchronisation points are required between the six face passes. The Quad2D kernel correctly uses `@shared` arrays to accumulate all face contributions before a single write to global memory.

**Impact:** The Hex3D pattern generates more global memory traffic and more synchronisation overhead than the Quad2D approach. At high polynomial degree with many warps, this can be a significant bottleneck.

**Fix approach:** Refactor `acousticsSurfaceHex3D.okl` to accumulate flux contributions into `@shared` arrays (matching the Quad2D pattern) and perform a single global write pass after all faces are processed.

---

### Report function copies device-to-host on every output interval

**File:** `solvers/acoustics/src/acousticsReport.cpp` lines 34–54

The `Report()` function calls `mesh.MassMatrixApply(o_q, o_Mq)` (a GPU kernel), then `platform.linAlg().innerProd()` (a GPU reduction), then conditionally `o_q.copyTo(q)` for file output. The mass-matrix application and inner product are efficient. The `o_q.copyTo(q)` transfers the full solution array for VTU output, which is unavoidable but is gated on `OUTPUT TO FILE = TRUE`.

No concern with normal runs; the pattern is consistent with other solvers.

---

### PlotFields writes ASCII VTU on rank-local files (O(Nelements) fprintf calls)

**File:** `solvers/acoustics/src/acousticsPlotFields.cpp`

The VTU writer calls `fprintf` per node (one call per plot node per element). For high-order elements with many plot nodes, this produces millions of individual `fprintf` calls per output, which is slow. The same issue exists in all sibling solvers that share this VTU writing pattern; it is a library-wide concern, not acoustics-specific.

---

## Fragile Areas

### Data file is the only mechanism for boundary and initial conditions

**Files:** `solvers/acoustics/data/acousticsGaussian2D.h`, `solvers/acoustics/data/acousticsGaussian3D.h`

Both initial conditions and boundary conditions are defined as C preprocessor macros in a single header that is compiled into the OCCA kernel at JIT time. There is only one 2D data file and one 3D data file. Adding a new test case requires creating a new `.h` file and changing the `DATA FILE` setting.

The 2D outflow comment in `acousticsGaussian2D.h` line 29 says `/* wall 1, outflow 2 */` while the 3D file says `/* wall 1, inflow 2 */` — inconsistent labelling for `bc==2`.

**Impact:** The `bc==2` boundary condition in 2D (line 31–34 of `acousticsGaussian2D.h`) sets `*(rB) = -rM` (anti-phase pressure) to suppress reflection, while 3D uses the same formula with the comment calling it `inflow`. The physical meaning of the two boundary types is inconsistently documented, and both headers use the same numerical treatment for `bc==2`, so the 3D label is misleading.

**Safe modification:** When adding a new problem, copy the closest existing `data/*.h` file and register it via a new `DATA FILE` setting. Do not change `bc==1` or `bc==2` semantics without updating all four surface kernels simultaneously.

---

## Summary Table

| Area | Severity | File(s) |
|------|----------|---------|
| Neumann BC not implemented (acknowledged in comment) | High | `okl/acousticsSurfaceTri2D.okl:118` |
| `time` implicit capture in device function | Medium | `okl/acousticsSurfaceHex3D.okl:102`, `okl/acousticsSurfaceQuad2D.okl:93` |
| Inconsistent flux convention (no documentation) | Medium | all `okl/acousticsSurface*.okl` |
| Wave speed hardcoded to 1.0 | Medium | `src/acousticsStep.cpp:31` |
| Dead kernel variants (`_v0`, `_v1`, `_v2`) | Low | `okl/acousticsVolumeTet3D.okl` |
| Hardcoded `p_Nvol`/`p_NblockV` in kernel body | Low | `okl/acousticsVolumeTet3D.okl:339-340` |
| `CFL NUMBER` missing from settings report | Low | `src/acousticsSettings.cpp:67-78` |
| No PML / absorbing boundary | Medium | solver-level gap |
| No convergence / error norm output | Medium | `src/acousticsRun.cpp`, `src/acousticsReport.cpp` |
| Hex3D surface: global writes + multiple barriers | Low-Medium | `okl/acousticsSurfaceHex3D.okl` |
| Inconsistent `bc==2` label in 2D vs 3D data files | Low | `data/acousticsGaussian2D.h:29`, `data/acousticsGaussian3D.h:30` |
| Undocumented `JW` weighting question | Low | `okl/acousticsVolumeHex3D.okl:70` |

---

*Concerns audit: 2026-04-22*
