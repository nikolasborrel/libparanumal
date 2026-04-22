# Codebase Structure: Acoustics Solver

**Analysis Date:** 2026-04-22

---

## Directory Layout

```
solvers/acoustics/
├── acoustics.hpp              # Solver class declaration + settings class
├── acousticsMain.cpp          # main() entry point
├── makefile                   # Build rules; produces acousticsMain or libacoustics.a
│
├── src/                       # C++ solver implementation files
│   ├── acousticsSetup.cpp     # acoustics_t::Setup() — mesh, kernels, time stepper
│   ├── acousticsRun.cpp       # acoustics_t::Run() — IC, dt selection, time loop
│   ├── acousticsStep.cpp      # acoustics_t::rhsf() + MaxWaveSpeed()
│   ├── acousticsReport.cpp    # acoustics_t::Report() — norm output + VTU dispatch
│   ├── acousticsPlotFields.cpp# acoustics_t::PlotFields() — VTU file writer
│   └── acousticsSettings.cpp  # acousticsSettings_t — setting definitions + parser
│
├── okl/                       # OCCA kernel files (JIT-compiled to GPU/CPU)
│   ├── acousticsVolumeTri2D.okl
│   ├── acousticsVolumeTet3D.okl
│   ├── acousticsVolumeQuad2D.okl
│   ├── acousticsVolumeHex3D.okl
│   ├── acousticsSurfaceTri2D.okl
│   ├── acousticsSurfaceTet3D.okl
│   ├── acousticsSurfaceQuad2D.okl
│   ├── acousticsSurfaceHex3D.okl
│   ├── acousticsInitialCondition2D.okl
│   └── acousticsInitialCondition3D.okl
│
├── setups/                    # Run-configuration files (.rc) for each element type
│   ├── setupTri2D.rc
│   ├── setupQuad2D.rc
│   ├── setupTet3D.rc
│   └── setupHex3D.rc
│
└── data/                      # C headers with IC and BC macro definitions
    ├── acousticsGaussian2D.h
    └── acousticsGaussian3D.h
```

---

## Header and Entry Point

### `acoustics.hpp`
The single solver header. Includes `core.hpp`, `platform.hpp`, `mesh.hpp`, `solver.hpp`, `timeStepper.hpp`, `linAlg.hpp`. Declares:
- `acousticsSettings_t` — key-value settings, parsed from `.rc` file
- `acoustics_t` — the solver class; inherits `solver_t`
- Macro `DACOUSTICS` = path to solver directory (used to locate OKL files at runtime)

### `acousticsMain.cpp`
Standard libParanumal entry point pattern:
1. `Comm::Init`
2. Parse settings from `argv[1]` (a `.rc` file)
3. Construct `platform_t`, `mesh_t`, `acoustics_t`
4. Call `acoustics.Run()`
5. `Comm::Finalize`

---

## Source Files (`src/`)

### `acousticsSettings.cpp`
Defines the settings registry for the acoustics solver. Each `newSetting()` call registers a key with default value and allowed values. The `parseFromFile` method reads all settings from a `.rc` file and routes keys to either `platformSettings`, `meshSettings`, or the solver's own settings.

**Registered settings:**

| Key | Default | Description |
|---|---|---|
| `DATA FILE` | `data/acousticsGaussian2D.h` | C header included into OKL kernels |
| `TIME INTEGRATOR` | `DOPRI5` | One of `AB3`, `DOPRI5`, `LSERK4` |
| `CFL NUMBER` | `1.0` | dt multiplier |
| `START TIME` | `0` | Integration start |
| `FINAL TIME` | `10` | Integration end |
| `OUTPUT INTERVAL` | `.1` | Time between Report() calls |
| `OUTPUT TO FILE` | `FALSE` | Write VTU files |
| `OUTPUT FILE NAME` | `acoustics` | VTU filename prefix |

### `acousticsSetup.cpp`
Implements `acoustics_t::Setup(platform, mesh, settings)`. Responsibilities:
- Allocates `q`, `o_q`, `o_Mq`
- Calls `mesh.HaloTraceSetup(Nfields)` to initialize MPI trace exchange
- Selects and constructs time stepper via `timeStepper.Setup<T>(...)`
- Builds OCCA `kernelInfo` (properties dict): copies `mesh.props`, adds includes (DATA FILE), and defines all `p_*` constants
- Calls `platform.buildKernel(...)` for `volumeKernel`, `surfaceKernel`, `initialConditionKernel`
- Kernel suffix (e.g. `Tri2D`) is chosen from `mesh.elementType`

### `acousticsRun.cpp`
Implements `acoustics_t::Run()`. Responsibilities:
- Fires `initialConditionKernel` to set `o_q` on device
- Computes time step `dt = cfl * hmin / (vmax*(N+1)^2)`
- Calls `timeStepper.Run(*this, o_q, startTime, finalTime)` which drives the time loop and periodically calls `Report()`
- After time loop: computes and prints solution L2 norm

### `acousticsStep.cpp`
Implements `acoustics_t::rhsf(o_Q, o_RHS, T)` and `acoustics_t::MaxWaveSpeed()`.

`MaxWaveSpeed()` returns `1.0` (constant wave speed, no material variation).

`rhsf` performs one full spatial residual evaluation with communication-computation overlap:
1. Start async halo exchange: `traceHalo.ExchangeStart(o_Q, 1)`
2. Launch `volumeKernel` on all elements (overlaps with MPI)
3. Launch `surfaceKernel` on internal elements
4. Finish halo exchange: `traceHalo.ExchangeFinish(o_Q, 1)`
5. Launch `surfaceKernel` on halo elements

### `acousticsReport.cpp`
Implements `acoustics_t::Report(time, tstep)`. Called by time stepper at output intervals. Computes `sqrt(q . M*q)` using `mesh.MassMatrixApply` + `platform.linAlg().innerProd`. If `OUTPUT TO FILE == TRUE`, copies `o_q` to host and calls `PlotFields`.

### `acousticsPlotFields.cpp`
Implements `acoustics_t::PlotFields(Q, fileName)`. Writes a per-MPI-rank VTU XML file. Fields written:
- `Density`: field `q[0]` (pressure-like `r`)
- `Velocity`: fields `q[1], q[2][, q[3]]` (velocity components)

Uses `mesh.PlotInterp` to interpolate from DG nodes to visualization nodes.

---

## OKL Kernel Files (`okl/`)

All kernel files use the [OCCA Kernel Language](https://libocca.org/) syntax. `@kernel` marks a GPU kernel entry point; `@outer`/`@inner` annotate the loop hierarchy mapped to thread blocks/threads. `@shared` allocates GPU shared memory; `@restrict` marks non-aliased pointer arguments.

### Volume Kernels

| File | Kernel | Element | Notes |
|---|---|---|---|
| `acousticsVolumeTri2D.okl` | `acousticsVolumeTri2D` | Triangle 2D | Dense `D[Np×Np]` per direction; constant geometry per element |
| `acousticsVolumeTet3D.okl` | `acousticsVolumeTet3D` | Tetrahedra 3D | Four variants (`_v0`–`_v2` + default); default uses blocked `p_NblockV=4` outer loop with shared-memory load of field variables |
| `acousticsVolumeQuad2D.okl` | `acousticsVolumeQuad2D` | Quad 2D | Tensor-product `DT[Nq×Nq]` applied along each direction; JW-weighted fluxes |
| `acousticsVolumeHex3D.okl` | `acousticsVolumeHex3D` | Hex 3D | Tensor-product `DT[Nq×Nq]` applied along r/s/t directions; 3D shared arrays `s_F[Nfields][Nq][Nq][Nq]` |

**Volume kernel physics** — all variants implement the same linear acoustic fluxes in reference coordinates:
```
Field 0 (density r):  F = -u, G = -v [, H = -w]
Field 1 (velocity u): F = -r, G = 0  [, H = 0 ]
Field 2 (velocity v): F = 0,  G = -r [, H = 0 ]
Field 3 (velocity w): F = 0,  G = 0  [, H = -r]  (3D only)
```

### Surface Kernels

| File | Kernel | Element | Notes |
|---|---|---|---|
| `acousticsSurfaceTri2D.okl` | `acousticsSurfaceTri2D` | Triangle 2D | Shared flux arrays `[NblockS][NfacesNfp]`; LIFT matrix applied in second inner loop |
| `acousticsSurfaceTet3D.okl` | `acousticsSurfaceTet3D` | Tetrahedra 3D | Same pattern as Tri2D extended to 4 fields |
| `acousticsSurfaceQuad2D.okl` | `acousticsSurfaceQuad2D` | Quad 2D | Processes face pairs (0&2, then 1&3) in separate inner loops; flux accumulated into shared `[NblockS][Nq][Nq]` arrays then written back |
| `acousticsSurfaceHex3D.okl` | `acousticsSurfaceHex3D` | Hex 3D | Processes face pairs (0&5, 1&3, 2&4) with `@barrier()` between passes; writes directly to global `rhsq` |

Each surface kernel file defines a local device function `upwind(nx, ny[, nz], rM, uM, ..., rP, uP, ..., *rflux, *uflux, ...)` implementing the exact Roe upwind flux for the acoustic system, then the BC lookup calls `acousticsDirichletConditions{2D|3D}(bc, ...)` before invoking `upwind`.

### Initial Condition Kernels

| File | Kernel | Calls |
|---|---|---|
| `acousticsInitialCondition2D.okl` | `acousticsInitialCondition2D` | `acousticsInitialConditions2D(t, x, y, &r, &u, &v)` |
| `acousticsInitialCondition3D.okl` | `acousticsInitialCondition3D` | `acousticsInitialConditions3D(t, x, y, z, &r, &u, &v, &w)` |

The actual IC values are provided by macros from the included data header; the kernel just sets `q[e*Np*Nfields + n + field*Np]` for each node.

---

## Setup / Data Files

### `.rc` Run Configuration Files (`setups/`)

Key-value text files parsed at startup. Each file maps to one element type:

| File | Element Type | Dim | Element code |
|---|---|---|---|
| `setupTri2D.rc` | Triangles | 2D | `[ELEMENT TYPE] 3` |
| `setupQuad2D.rc` | Quadrilaterals | 2D | `[ELEMENT TYPE] 4` |
| `setupTet3D.rc` | Tetrahedra | 3D | `[ELEMENT TYPE] 6` |
| `setupHex3D.rc` | Hexahedra | 3D | `[ELEMENT TYPE] 12` |

All four use `[MESH FILE] BOX` (built-in box mesh generator) and `[TIME INTEGRATOR] DOPRI5`. Box sizes differ: Tri2D/Quad2D use 10×10, Tet3D uses 2×2×2 (coarser due to 3D cost), Hex3D uses 5×5×5.

### Data Headers (`data/`)

C headers included directly into OKL kernels at JIT compile time via the `kernelInfo["includes"]` mechanism. Each defines two macros:

**`acousticsGaussian2D.h`**
```c
// Boundary conditions (bc index → ghost state)
#define acousticsDirichletConditions2D(bc, t, x, y, nx, ny, rM, uM, vM, rB, uB, vB) { ... }
// bc==1: wall (reflecting), bc==2: outflow (absorbing)

// Initial conditions
#define acousticsInitialConditions2D(t, x, y, r, u, v) {
  *(r) = 1.0 + exp(-3*(x*x+y*y));   // Gaussian pressure pulse
  *(u) = 0.0; *(v) = 0.0;
}
```

**`acousticsGaussian3D.h`** — same structure with 3D extension (`w`, `z`, `nz`).

To implement a new problem: create a new `.h` file in `data/` and point `[DATA FILE]` in the `.rc` setup to it.

---

## Shared Library Infrastructure (`include/`)

Files in `include/` that the acoustics solver directly uses:

| File | Purpose |
|---|---|
| `include/solver.hpp` | `solver_t` base class; defines `rhsf`, `Run`, `Report` interface |
| `include/timeStepper.hpp` | `timeStepper_t` wrapper + all concrete stepper classes (`ab3`, `lserk4`, `dopri5`, etc.) |
| `include/mesh.hpp` | `mesh_t`: topology, geometry, DG operators, MPI partitioning |
| `include/platform.hpp` | `platform_t`: OCCA device, kernel builder, memory allocator |
| `include/settings.hpp` | `settings_t`: key-value registry, `.rc` file parser |
| `include/linAlg.hpp` | GPU linear algebra kernels (`innerProd`, etc.) |
| `include/ogs.hpp` | Overlapping Gather-Scatter; provides halo exchange (`halo_t`) |
| `include/core.hpp` | Common types: `dfloat`, `dlong`, `hlong`, `memory<T>`, `deviceMemory<T>` |

---

## Naming Conventions

### Files
- `acoustics{Action}{ElementType}.{ext}` — for OKL kernels, e.g. `acousticsVolumeTri2D.okl`, `acousticsSurfaceHex3D.okl`
- `acoustics{Concern}.cpp` — for source files, e.g. `acousticsSetup.cpp`, `acousticsRun.cpp`
- `setup{ElementType}.rc` — for run configs, e.g. `setupQuad2D.rc`
- `acoustics{IC}{Dim}.h` — for data headers, e.g. `acousticsGaussian2D.h`

### Variables and Fields
- Field state array: `q` (host) / `o_q` (device); indexed as `q[e*Np*Nfields + n + field*Np]`
- Field indices: `0=r` (density/pressure), `1=u`, `2=v`, `3=w` (3D only)
- RHS array: `rhsq` / `o_RHS` / `o_rhs`
- Device arrays: prefixed `o_` (OCCA convention)
- Geometry arrays: `vgeo` (volume), `sgeo` (surface)
- Kernel defines: `p_` prefix (e.g. `p_Nfields`, `p_Np`, `p_NblockV`)

### Element Type Suffixes
| Suffix | Type | Dim |
|---|---|---|
| `Tri2D` | Triangle | 2 |
| `Quad2D` | Quadrilateral | 2 |
| `Tet3D` | Tetrahedron | 3 |
| `Hex3D` | Hexahedron | 3 |

---

## Where to Add New Code

**New problem (different IC/BC):**
- Create `data/acoustics{ProblemName}{2D|3D}.h` with `acousticsInitialConditions{2D|3D}` and `acousticsDirichletConditions{2D|3D}` macros
- Point `[DATA FILE]` in the setup `.rc` to the new header

**New setup configuration:**
- Copy an existing `setups/setup{Element}.rc` and modify parameters
- Run with `./acousticsMain setups/myNewSetup.rc`

**New time integrator:**
- Add `timeStepper.Setup<TimeStepper::newScheme>(...)` branch in `src/acousticsSetup.cpp`
- The new scheme class must already exist in `include/timeStepper.hpp`

**New output field:**
- Add field computation in `src/acousticsReport.cpp` or `src/acousticsPlotFields.cpp`

**New volume/surface kernel variant:**
- Create `okl/acousticsVolume{NewType}.okl` and `okl/acousticsSurface{NewType}.okl`
- Add a new `suffix` branch in `src/acousticsSetup.cpp`
- Add element type to `Mesh::ElementType` enum in `include/mesh.hpp`

---

*Structure analysis: 2026-04-22*
