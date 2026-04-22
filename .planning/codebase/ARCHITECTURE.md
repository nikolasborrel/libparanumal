# Architecture: Acoustics Solver

**Analysis Date:** 2026-04-22

## Overview

The acoustics solver implements isotropic linear acoustics using a nodal **Discontinuous Galerkin (DG)** spatial discretization combined with a **method-of-lines** temporal approach. The PDE is cast as a first-order hyperbolic system for pressure-like density `r` and velocity components `(u, v[, w])`. The spatial residual is evaluated on a GPU via OCCA kernel pairs (volume + surface), and the resulting ODE `dq/dt = rhsf(q,t)` is integrated forward in time by a pluggable time-stepper.

The design follows the libParanumal solver pattern: every solver is a subclass of `solver_t` (which exposes `rhsf`, `Run`, `Report`) and delegates time advancement entirely to `timeStepper_t`.

---

## Governing Equations

The linearized acoustic system solved is (2D example):

```
dr/dt = -(du/dx + dv/dy)
du/dt = -dr/dx
dv/dt = -dr/dy
```

State vector `q = [r, u, v]` (2D) or `q = [r, u, v, w]` (3D). Wave speed is fixed at `c = 1` (see `MaxWaveSpeed()` in `src/acousticsStep.cpp`).

---

## Class Hierarchy

### `acousticsSettings_t` — `solvers/acoustics/acoustics.hpp`
- Inherits: `settings_t` (from `include/settings.hpp`)
- Parses a `.rc` run-configuration file; distributes keys to `platformSettings_t` and `meshSettings_t`
- Owns settings: `DATA FILE`, `TIME INTEGRATOR`, `CFL NUMBER`, `START TIME`, `FINAL TIME`, `OUTPUT INTERVAL`, `OUTPUT TO FILE`, `OUTPUT FILE NAME`

### `acoustics_t` — `solvers/acoustics/acoustics.hpp`
- Inherits: `solver_t` → `operator_t` (from `include/solver.hpp`)
- Central solver object; owns all field data and GPU kernels
- Key members:

| Member | Type | Role |
|---|---|---|
| `mesh` | `mesh_t` | Mesh, geometry, DG operators |
| `Nfields` | `int` | 3 (2D) or 4 (3D) |
| `timeStepper` | `timeStepper_t` | Wraps chosen RK/AB scheme |
| `traceHalo` | `ogs::halo_t` | MPI halo exchange for face traces |
| `q` / `o_q` | `memory<dfloat>` / `deviceMemory<dfloat>` | Host/device solution array |
| `o_Mq` | `deviceMemory<dfloat>` | Temporary for mass-matrix product (norms) |
| `volumeKernel` | `kernel_t` | OCCA volume DG kernel |
| `surfaceKernel` | `kernel_t` | OCCA surface DG + upwind flux kernel |
| `initialConditionKernel` | `kernel_t` | OCCA IC kernel |

### `solver_t` — `include/solver.hpp`
- Abstract base; defines the interface `rhsf`, `rhs_imex_f/g`, `rhs_subcycle_f`, `Run`, `Report`, `Operator`
- The acoustics solver only overrides `rhsf`, `Run`, and `Report`

### `timeStepper_t` — `include/timeStepper.hpp`
- Thin wrapper around a `shared_ptr<timeStepperBase_t>`
- Exposes `Run(solver, o_q, start, end)` and `SetTimeStep(dt)`
- Concrete implementations available: `ab3`, `lserk4`, `dopri5` (and semi-analytic / PML variants not used here)

### `mesh_t` — `include/mesh.hpp`
- Provides all mesh topology and DG geometric arrays on host and device
- Used arrays in acoustics: `o_vgeo` (volume geometry), `o_sgeo` (surface geometry), `o_D`/`o_DT` (derivative matrices), `o_LIFT` (surface-to-volume lift), `o_vmapM`/`o_vmapP` (face node maps), `o_EToB` (boundary tags), `o_x/y/z` (physical coordinates), `o_internalElementIds`/`o_haloElementIds`

---

## Setup Flow (`src/acousticsSetup.cpp`)

```
acoustics_t::Setup(platform, mesh, settings)
  ├── Set Nfields = (dim==3) ? 4 : 3
  ├── ogs::InitializeKernels()          // JIT kernel build for OGS
  ├── platform.linAlg().InitKernels()   // innerProd kernel
  ├── mesh.HaloTraceSetup(Nfields)      // MPI halo setup
  ├── timeStepper.Setup<AB3|LSERK4|DOPRI5>(...)
  ├── malloc q + o_q (local + halo entries)
  ├── malloc o_Mq; mesh.MassMatrixKernelSetup(Nfields)
  ├── Build kernelInfo (OCCA properties)
  │     ├── include DATA FILE header (BCs + ICs)
  │     ├── define p_Nfields, p_half, p_maxNodes
  │     ├── set p_NblockV, p_NblockS (thread-block sizes)
  │     └── define p_Lambda2 = 0.5
  ├── platform.buildKernel(acousticsVolume{suffix})
  ├── platform.buildKernel(acousticsSurface{suffix})
  └── platform.buildKernel(acousticsInitialCondition{2D|3D})
```

The `suffix` is one of `Tri2D`, `Quad2D`, `Tet3D`, `Hex3D`, selected at runtime from `mesh.elementType`.

---

## Run Flow (`src/acousticsRun.cpp`)

```
acoustics_t::Run()
  ├── initialConditionKernel(Nelements, t0, o_x, o_y, o_z, o_q)
  ├── Compute dt = cfl * hmin / (vmax*(N+1)^2)
  │     vmax = 1.0 (fixed constant in MaxWaveSpeed())
  ├── timeStepper.SetTimeStep(dt)
  └── timeStepper.Run(*this, o_q, startTime, finalTime)
        └── calls this->rhsf(o_q, o_rhs, t) at each stage
              └── (see RHS Evaluation below)
```

After time integration, the solution L2 norm is printed via `mesh.MassMatrixApply` + `linAlg().innerProd`.

---

## RHS Evaluation (`src/acousticsStep.cpp`)

`acoustics_t::rhsf(o_Q, o_RHS, T)` implements one full spatial residual evaluation:

```
1. traceHalo.ExchangeStart(o_Q, 1)       // begin async MPI halo exchange

2. volumeKernel(Nelements,
                mesh.o_vgeo, mesh.o_D[/DT],
                o_Q, o_RHS)              // volume DG on ALL elements (can overlap with MPI)

3. surfaceKernel(NinternalElements,      // internal elements (no halo dependency)
                 o_internalElementIds,
                 mesh.o_sgeo, mesh.o_LIFT,
                 mesh.o_vmapM, mesh.o_vmapP,
                 mesh.o_EToB, T,
                 mesh.o_x/y/z, o_Q, o_RHS)

4. traceHalo.ExchangeFinish(o_Q, 1)     // complete halo exchange

5. surfaceKernel(NhaloElements,          // halo elements (need completed halo data)
                 o_haloElementIds, ...)
```

This overlaps communication (step 1/4) with computation (step 2/3), a standard latency-hiding pattern in DG codes.

---

## GPU Kernel Dispatch

Kernels are compiled JIT by OCCA at first run and cached in `.occa/`. The dispatch signature for volume and surface kernels follows the element-type naming convention. Thread-block sizes are baked in at build time via `p_NblockV` (volume) and `p_NblockS` (surface) defines:

```
p_NblockV = max(1, blockMax / Np)
p_NblockS = max(1, blockMax / maxNodes)
blockMax   = 256 (CPU/HIP) or 512 (CUDA)
```

The OCCA `@outer` loop iterates over element blocks; `@inner` loops over nodes within an element.

---

## DG Spatial Discretization

### Volume kernel
Computes `D^T * F(q)` in reference coordinates for each element. For simplex elements (Tri2D, Tet3D), a dense differentiation matrix `D[Np×Np]` is applied. For tensor-product elements (Quad2D, Hex3D), a 1D differentiation matrix `DT[Nq×Nq]` is applied direction-by-direction using shared memory.

Physical fluxes for acoustics (2D):
```
F[0] = -u    G[0] = -v      (density eq)
F[1] = -r    G[1] = 0       (u-momentum eq)
F[2] = 0     G[2] = -r      (v-momentum eq)
```

### Surface kernel
Computes the upwind numerical flux at each face node and accumulates via LIFT matrix:

```
rhs += LIFT * (sJ/J) * flux_correction
```

The upwind (Roe) flux for the acoustic system:
```c
ndotUM = nx*uM + ny*vM
ndotUP = nx*uP + ny*vP
rflux  = 0.5 * (ndotUP - ndotUM - (rP - rM))
uflux  = 0.5 * nx * ((rP - rM) - (ndotUP - ndotUM))
vflux  = 0.5 * ny * ((rP - rM) - (ndotUP - ndotUM))
```

The `p_Lambda2 = 0.5` constant encodes the upwind parameter.

For Quad2D and Hex3D, the surface kernel processes each pair of opposite faces in separate `@outer`-indexed loops with `@barrier()` calls between them, avoiding race conditions when writing to global `rhsq`.

---

## Boundary Condition Handling

Boundary conditions are defined entirely in the **data header** included at kernel build time (e.g., `data/acousticsGaussian2D.h`). The surface kernel looks up the BC tag `EToB[face + Nfaces*element]`:

- `bc == 0`: interior face — use neighbor trace directly
- `bc > 0`: boundary face — call `acousticsDirichletConditions{2D|3D}(bc, t, x, y[,z], n..., qM..., &qP...)` which overwrites the exterior trace `qP`

The standard BCs defined in the Gaussian data files:
- `bc == 1`: **wall** (reflecting) — mirror velocity normal component: `uP = uM - 2*(n·uM)*n`
- `bc == 2`: **outflow** (absorbing) — flip density sign: `rP = -rM`, keep velocities

To add new boundary conditions: add a new `bc == N` branch inside the macro in the data header file. No C++ source changes required.

---

## Time-Stepping Schemes

All three schemes supported by the acoustics solver are explicit:

| Setting key | Class | Notes |
|---|---|---|
| `AB3` | `TimeStepper::ab3` | 3rd-order Adams-Bashforth, fixed step |
| `LSERK4` | `TimeStepper::lserk4` | Low-storage 4th-order Runge-Kutta, fixed step |
| `DOPRI5` | `TimeStepper::dopri5` | Dormand-Prince order 5(4), adaptive step (default) |

The time step for fixed-step schemes is initialized as:
```
dt = CFL * hmin / (vmax * (N+1)^2)
```
where `hmin = mesh.MinCharacteristicLength()` and `vmax = 1.0`.

`DOPRI5` uses this as the initial step guess and adapts using embedded error control with absolute tolerance `ATOL` and relative tolerance `RTOL`.

---

## Reporting and Output

`acoustics_t::Report(time, tstep)` — called by the time stepper at intervals controlled by `OUTPUT INTERVAL`:
1. Computes weighted L2 norm: `sqrt(q . M*q)` where M is the DG mass matrix
2. Prints to stdout: `time (tstep), norm`
3. If `OUTPUT TO FILE == TRUE`: copies `o_q` to host, writes VTU file `{name}_{rank}_{frame}.vtu`

VTU output (`src/acousticsPlotFields.cpp`) interpolates fields to visualization nodes using `mesh.PlotInterp` and writes pressure (`r`) and velocity (`u,v[,w]`) as separate `DataArray` entries in the UnstructuredGrid format (VTK cell type 5=triangle for 2D, 10=tetrahedron for 3D).

---

## Solver Parameterization

The solver is fully parameterized via `.rc` files and the choice of data header:

| Dimension of variation | Mechanism |
|---|---|
| Element type | `[ELEMENT TYPE]` + `[MESH DIMENSION]` in `.rc` → selects OKL kernel suffix |
| Polynomial degree | `[POLYNOMIAL DEGREE]` → `mesh.N`, sizes `Np`, `Nq` baked into kernel defines |
| Time integrator | `[TIME INTEGRATOR]` → `AB3`, `LSERK4`, or `DOPRI5` |
| IC + BC functions | `[DATA FILE]` → C header included into OKL kernels at JIT compile |
| Mesh | `[MESH FILE] BOX` with `[BOX NX/NY/NZ]`, or external mesh file |
| Hardware backend | `[THREAD MODEL]` → CUDA, HIP, OpenCL, OpenMP, Serial |

---

*Architecture analysis: 2026-04-22*
