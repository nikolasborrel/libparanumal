# External Integrations

**Analysis Date:** 2026-04-22

## How the Acoustics Solver Integrates with libParanumal

The acoustics solver lives at `solvers/acoustics/` and is one of eight solvers (acoustics, advection, bns, cns, elliptic, fokkerPlanck, gradient, ins) in the libParanumal framework. It is not a standalone library; it depends on the shared framework stack built in `libs/`.

### Inheritance Chain

```
libp::operator_t          (include/operator.hpp)
  └── libp::solver_t      (include/solver.hpp)
        └── acoustics_t   (solvers/acoustics/acoustics.hpp)
```

`solver_t` holds a `platform_t`, `settings_t`, and `comm_t`. `acoustics_t` adds a `mesh_t`, `timeStepper_t`, `ogs::halo_t`, device/host field arrays, and OCCA kernel handles.

### Lifecycle

1. `main()` in `acousticsMain.cpp` initialises MPI, constructs `platformSettings_t`, `meshSettings_t`, and `acousticsSettings_t`, then calls `acousticsSettings.parseFromFile()` which reads one flat `.rc` key-value file and routes each key to the appropriate settings object.
2. `platform_t platform(platformSettings)` — initialises OCCA device, `linAlg_t`, and the OCCA cache.
3. `mesh_t mesh(platform, meshSettings, comm)` — loads/generates mesh, sets up connectivity, quadrature data, and uploads all geometric factors to device memory.
4. `acoustics_t acoustics(platform, mesh, acousticsSettings)` — calls `acoustics_t::Setup(...)`, which JIT-compiles OKL kernels via `platform.buildKernel(...)`, sets up the halo exchange, and configures the time stepper.
5. `acoustics.Run()` — sets initial conditions on device, computes CFL timestep, then delegates to `timeStepper.Run(*this, o_q, startTime, finalTime)`.

---

## Key Library Dependencies

### OCCA (`occa/` submodule)

All GPU/accelerator compute goes through OCCA. The acoustics solver touches OCCA through the `platform_t` wrapper.

- **Kernel compilation:** `platform.buildKernel(fileName, kernelName, kernelInfo)` — compiles `.okl` files at runtime (JIT) and caches results in `${LIBP_DIR}/.occa/` (cleaned via `make clean-kernels`).
- **Device memory:** `platform.malloc<dfloat>(count)` → `deviceMemory<dfloat>` (wraps `occa::memory`). The convention is that device arrays are prefixed `o_` (e.g., `o_q`, `o_Mq`).
- **Kernel launches:** `volumeKernel(Nelements, mesh.o_vgeo, mesh.o_D, o_Q, o_RHS)` — OCCA kernel objects are callable.
- **Backend selection:** Set at runtime via `[THREAD MODEL]` key in `.rc` file. Supported values: `CUDA`, `HIP`, `OpenCL`, `OpenMP`, `Serial`.
- **Block size tuning:** In `acousticsSetup.cpp`, `blockMax = 256` for non-CUDA backends and `512` for CUDA, checked via `platform.device.mode() == "CUDA"`.

### OGS / Halo (`libs/ogs/`, `include/ogs.hpp`)

`ogs::halo_t traceHalo` handles inter-process exchange of DG face trace data.

- **Setup:** `traceHalo = mesh.HaloTraceSetup(Nfields)` — called in `acoustics_t::Setup`.
- **Usage in `rhsf` (`src/acousticsStep.cpp`):**
  ```
  traceHalo.ExchangeStart(o_Q, 1);    // start async MPI + device packing
  volumeKernel(...);                  // overlap with volume computation
  traceHalo.ExchangeFinish(o_Q, 1);  // complete exchange
  surfaceKernel(...);                 // then process halo surface fluxes
  ```
  The split `ExchangeStart`/`ExchangeFinish` pattern hides MPI latency behind the volume kernel.
- **Kernel initialisation:** `ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add)` is called explicitly in `Setup` to JIT-compile OGS exchange kernels before first use.

### linAlg (`libs/linAlg/`, `include/linAlg.hpp`)

Used for global inner-product reductions.

- **Setup:** `platform.linAlg().InitKernels({"innerProd"})` — in `acoustics_t::Setup`.
- **Usage:** `platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm)` — in `Run()` and `Report()` to compute the L2 norm `sqrt(q^T M q)`. The inner product performs an `MPI_Allreduce` internally across `mesh.comm`.

### mesh (`libs/mesh/`, `include/mesh.hpp`)

The mesh object is the central data store. The acoustics solver reads but never writes mesh data after construction.

Key mesh members consumed by acoustics:
- `mesh.Nelements`, `mesh.Np`, `mesh.Nfp`, `mesh.Nfaces`, `mesh.dim`, `mesh.elementType`, `mesh.rank`
- `mesh.o_vgeo` — volume geometric factors (Jacobians, metric terms)
- `mesh.o_D` — differentiation matrix
- `mesh.o_sgeo` — surface geometric factors (normals, surface Jacobians)
- `mesh.o_LIFT` — DG lifting matrix
- `mesh.o_vmapM`, `mesh.o_vmapP` — face-to-volume node maps (interior/exterior)
- `mesh.o_EToB` — element-to-boundary-condition type map
- `mesh.o_x`, `mesh.o_y`, `mesh.o_z` — node coordinates (used in BC evaluation and initial conditions)
- `mesh.o_internalElementIds`, `mesh.o_haloElementIds` — partition lists for overlap

Mesh also provides:
- `mesh.MassMatrixKernelSetup(Nfields)` — JIT-compiles mass matrix apply kernel.
- `mesh.MassMatrixApply(o_q, o_Mq)` — applies mass matrix on device.
- `mesh.MinCharacteristicLength()` — returns minimum element size for CFL computation.
- `mesh.HaloTraceSetup(Nfields)` — returns an `ogs::halo_t`.
- `mesh.PlotInterp(...)` and `mesh.plotNp` — CPU-side interpolation to plot nodes (used in VTU output).

### timeStepper (`libs/timeStepper/`, `include/timeStepper.hpp`)

The time integration loop is fully encapsulated. The solver implements `rhsf(o_q, o_rhs, time)` and the stepper calls it.

- **Setup pattern:** `timeStepper.Setup<TimeStepper::dopri5>(Nelements, totalHaloPairs, Np, Nfields, platform, comm)` — templated on the integrator type.
- **Available integrators (as configured in `acousticsSettings.cpp`):**
  - `TimeStepper::ab3` — Adams-Bashforth order 3
  - `TimeStepper::lserk4` — Low-Storage Explicit Runge-Kutta order 4
  - `TimeStepper::dopri5` — Dormand-Prince adaptive RK45 (default)
- **Run interface:** `timeStepper.Run(*this, o_q, startTime, finalTime)` — the stepper calls `solver.rhsf(o_q, o_rhs, t)` at each stage. `dopri5` performs adaptive step-size control via error estimate.
- **Note:** PML and multi-rate variants (`ab3_pml`, `mrab3`, `sark4_pml`, etc.) are available in the framework but are not used by the acoustics solver.

### parAdogs (`libs/parAdogs/`, `include/parAdogs.hpp`)

Parallel mesh partitioning and reordering. Used internally by `mesh_t` construction; the acoustics solver has no direct calls to parAdogs.

---

## Data Flow Between Solver and Framework Layers

```
.rc file (disk)
    │
    ▼
acousticsSettings_t::parseFromFile()
    │ routes keys to platformSettings / meshSettings / acousticsSettings
    ▼
platform_t  ─────────────────────────────────────────────────────┐
    │ OCCA device handle, linAlg kernels, kernel cache            │
    ▼                                                             │
mesh_t  ──────────────────────────────────────────────────────── │
    │ geometry, quadrature, connectivity → all on device (o_*)   │
    ▼                                                             │
acoustics_t::Setup()                                             │
    │                                                             │
    ├─ ogs::InitializeKernels(platform, ...)                      │
    ├─ linAlg().InitKernels({"innerProd"})                        │
    ├─ traceHalo = mesh.HaloTraceSetup(Nfields)  ─ OGS halo      │
    ├─ timeStepper.Setup<dopri5>(...)                             │
    ├─ o_q  = platform.malloc<dfloat>(...)  ─ solution on device │
    ├─ o_Mq = platform.malloc<dfloat>(...)  ─ M*q on device      │
    └─ volumeKernel / surfaceKernel / initialConditionKernel      │
         compiled via platform.buildKernel(...)  ◄────────────── ┘
    ▼
acoustics_t::Run()
    │
    ├─ initialConditionKernel(...)  ─ sets o_q on device
    ├─ mesh.MinCharacteristicLength() + MaxWaveSpeed() → dt
    ├─ timeStepper.SetTimeStep(dt)
    └─ timeStepper.Run(*this, o_q, t0, tf)
              │ calls solver.rhsf() each stage
              ▼
         acoustics_t::rhsf(o_Q, o_RHS, T)
              ├─ traceHalo.ExchangeStart(o_Q, 1)    ─ async MPI
              ├─ volumeKernel(...)                   ─ GPU volume integral
              ├─ surfaceKernel(NinternalElements,...) ─ GPU internal surface
              ├─ traceHalo.ExchangeFinish(o_Q, 1)   ─ complete halo
              └─ surfaceKernel(NhaloElements,...)    ─ GPU halo surface
    ▼
acoustics_t::Report(time, tstep)   [called by timeStepper at OUTPUT INTERVAL]
    ├─ mesh.MassMatrixApply(o_q, o_Mq)
    ├─ linAlg().innerProd(...)  → global L2 norm via MPI_Allreduce
    └─ if OUTPUT TO FILE: o_q.copyTo(q), acoustics_t::PlotFields(q, fname)
```

---

## External Interfaces

### MPI

- **Initialisation:** `Comm::Init(argc, argv)` / `Comm::Finalize()` called in `acousticsMain.cpp`.
- **Communicator:** `comm_t comm(Comm::World().Dup())` — a duplicate of `MPI_COMM_WORLD`.
- **Communication in solver:** exclusively through OGS halo exchange (`traceHalo.ExchangeStart/Finish`) and `linAlg().innerProd()` (which calls `MPI_Allreduce` internally). The solver itself never calls MPI directly.
- **Multi-rank testing:** `test/testAcoustics.py` runs `testAcousticsTri_MPI` with `ranks=4` using `mpirun`.

### File I/O (Mesh Input)

- **Built-in box mesh:** When `[MESH FILE]` is set to `BOX`, the mesh is generated internally from `[BOX NX/NY/NZ]` and `[MESH DIMENSION]` parameters. No external file needed.
- **External mesh files:** The `mesh_t` constructor also accepts `.msh` (Gmsh format) mesh files. Example meshes are in `test/` (`squareTri.msh`, `cubeHex.msh`, etc.).

### File I/O (Solution Output)

- **VTU files:** When `[OUTPUT TO FILE] = TRUE`, `acoustics_t::Report()` writes one `.vtu` file per MPI rank per output interval. Files are named `{OUTPUT FILE NAME}_{rank:04d}_{frame:04d}.vtu` and use the VTK UnstructuredGrid XML format (ASCII, version 0.1, BigEndian byte order). Written with C `FILE*` in `src/acousticsPlotFields.cpp`. Fields exported: `Density` (scalar), `Velocity` (vector, dim components).
- **Stdout:** Solution norm is printed to stdout from rank 0 only: `printf("%5.2f (%d), %5.2f (time, timestep, norm)\n", ...)`.

### Settings / Configuration

- **Format:** Plain-text key-value `.rc` files (version 2.0 format). Keys are uppercase in brackets, e.g. `[TIME INTEGRATOR]`. See `solvers/acoustics/setups/setupTri2D.rc` for a complete example.
- **Parsing:** `acousticsSettings_t::parseFromFile()` reads all keys and routes them to `platformSettings`, `meshSettings`, or `acousticsSettings` based on which object `hasSetting(name)`. Unknown keys cause `LIBP_FORCE_ABORT`.
- **Boundary/IC header:** `[DATA FILE]` key points to a C header (`.h`) that defines the macros `acousticsDirichletConditions2D/3D` and `acousticsInitialConditions2D/3D`. This file is included directly into OKL kernels at JIT-compile time via `kernelInfo["includes"] += dataFileName`. Default: `data/acousticsGaussian2D.h`.

### OCCA Kernel Cache

- Compiled kernels are cached in `${LIBP_DIR}/.occa/` (project-local) or `~/.occa/` (user-global, cleaned by `make clean-kernels`).

---

## What the Acoustics Solver Does NOT Use

- `libp::linearSolver` (`libs/linearSolver/`) — no implicit solves.
- `libp::parAlmond` (`libs/parAlmond/`) — no AMG preconditioner.
- `libp::initialGuess` (`include/initialGuess.hpp`) — not applicable (explicit time stepping).
- The IMEX, sub-cycle, PML, or multi-rate `solver_t` interfaces (`rhs_imex_f`, `rhsf_pml`, `rhsf_MR`, etc.) — only `rhsf` is implemented.

---

*Integration audit: 2026-04-22*
