# Coding Conventions — acoustics solver

**Analysis Date:** 2026-04-22
**Scope:** `solvers/acoustics/` within the libParanumal HPC finite element library

---

## Naming Patterns

### Files
- C++ source files: `acoustics{Action}.cpp` — e.g., `acousticsSetup.cpp`, `acousticsStep.cpp`, `acousticsRun.cpp`, `acousticsReport.cpp`, `acousticsPlotFields.cpp`
- C++ header: `acoustics.hpp` at solver root
- OKL kernel files: `acoustics{Operator}{ElementType}.okl` — e.g., `acousticsVolumeTri2D.okl`, `acousticsSurfaceHex3D.okl`, `acousticsInitialCondition2D.okl`
- Data/BC files: `acoustics{Case}{Dim}.h` — e.g., `acousticsGaussian2D.h`, `acousticsGaussian3D.h`
- Setup/config files: `setup{ElementType}.rc` — e.g., `setupTri2D.rc`, `setupHex3D.rc`
- Main entry point: `acousticsMain.cpp`

### Classes and Types
- Solver class: `acoustics_t` (suffix `_t` for all libParanumal types)
- Settings class: `acousticsSettings_t`
- Typedef aliases in `libp` namespace: `dfloat` (double or float), `dlong` (device-compatible long), `memory<T>` (host memory), `deviceMemory<T>` (device memory)

### Functions / Methods
- Setup entry: `Setup(platform_t&, mesh_t&, acousticsSettings_t&)`
- Time-integration RHS: `rhsf(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T)` — uppercase `Q`, `RHS`, `T` signal device/time arguments
- Other methods use PascalCase: `Run()`, `Report()`, `PlotFields()`, `MaxWaveSpeed()`
- Kernel builder calls follow the pattern: `platform.buildKernel(fileName, kernelName, kernelInfo)`

### Variables — Host Side
- Device memory objects are prefixed with `o_`: `o_q`, `o_Mq`, `o_Q`, `o_RHS`
- Host memory objects have no prefix: `q`, `scratch`, `Ix`, `Iy`
- Count variables: `Nfields`, `Nlocal`, `Nhalo`, `Nentries`
- Index bases: `qbase`, `gbase`, `base`, `id`, `idM`, `idP`
- Block-size tuning: `NblockV` (volume), `NblockS` (surface), `blockMax`
- Penalty parameter: `Lambda2`
- Time step: `dt`, `cfl`, `hmin`, `vmax`

### Settings Keys (`.rc` format)
Settings use ALL-CAPS bracket syntax: `[DATA FILE]`, `[TIME INTEGRATOR]`, `[POLYNOMIAL DEGREE]`, `[THREAD MODEL]`, `[BOX NX]`, etc. Keys are strings compared with `settings.compareSetting(...)` or retrieved via `settings.getSetting(...)`.

---

## Physics Variable Naming

The acoustic state vector `q` stores fields in field-index order (not interleaved). The layout within an element block is:

| Field index | 2D name | 3D name | Physical meaning      |
|-------------|---------|---------|----------------------|
| 0           | `r`     | `r`     | acoustic pressure (density perturbation) |
| 1           | `u`     | `u`     | x-velocity           |
| 2           | `v`     | `v`     | y-velocity           |
| 3           | —       | `w`     | z-velocity (3D only) |

In OKL kernels the physical fields are loaded as:
```c
const dfloat r = q[qbase + 0*p_Np];
const dfloat u = q[qbase + 1*p_Np];
const dfloat v = q[qbase + 2*p_Np];
// 3D only:
const dfloat w = q[qbase + 3*p_Np];
```

In the Tet3D alternative kernels (`acousticsVolumeTet3D_v1`, `_v2`) the same field is called `rho`/`s_rho` in shared-memory staging arrays (reflecting its physical interpretation as acoustic density), but the primary (`_v0` and the active `acousticsVolumeTet3D`) kernel uses `r`.

Flux variables: `f`, `g`, `h` for the physical-space flux components in the r-, s-, t-reference directions.

Riemann solver outputs: `rflux`, `uflux`, `vflux`, `wflux` (one per field). Normal-velocity dot-products: `ndotUM`, `ndotUP`.

Boundary condition ghost values use the `P` suffix (`rP`, `uP`, `vP`, `wP`); interior values use `M` suffix (`rM`, `uM`, `vM`, `wM`).

---

## Data Layout

### State Vector (q)
**Field-major AoS-within-element layout.** For element `e`, node `n`, field `f`:

```
index = e * p_Np * p_Nfields + f * p_Np + n
```

In code: `q[e*p_Np*p_Nfields + n + f*p_Np]` where the node offset `n` is the inner index and field stride is `p_Np`. This is an Array-of-Structures variant where within each element the nodes are contiguous per field, not per node.

For Quad/Hex tensor-product elements the node index `n` maps to `(k*p_Nq*p_Nq + j*p_Nq + i)` with i as the innermost (fastest-varying) dimension.

### Geometry Arrays
- `vgeo`: volume geometry. For simplex elements: `vgeo[e*p_Nvgeo + geoID]` (element-stride, constant per element). For tensor-product elements: `vgeo[e*p_Np*p_Nvgeo + node_index + p_Np*geoID]` (per-node varying factors with field-stride interleave).
- `sgeo`: surface geometry, `sgeo[sk*p_Nsgeo + geoID]` where `sk` is the surface node index.

---

## OKL Kernel Coding Patterns

OKL (OCCA Kernel Language) is C-like with OCCA parallelism annotations. Conventions used throughout the acoustics kernels:

### Parallelism Structure
```c
// Outer loop = one threadblock per element (simplex) or per element-batch (tensor)
for(dlong e=0;e<Nelements;++e;@outer(0)){
  @shared dfloat s_F[p_Nfields][p_Np];    // shared memory for flux staging

  // Inner loop = one thread per node
  for(int n=0;n<p_Np;++n;@inner(0)){
    ...
  }
}
```

For tensor-product elements (Quad2D, Hex3D) the inner loops are 2D or 3D matching the tensor structure:
```c
for(int j=0;j<p_Nq;++j;@inner(1)){
  for(int i=0;i<p_Nq;++i;@inner(0)){
    ...
  }
}
```

Surface kernels batch multiple elements per block using `p_NblockS`:
```c
for(dlong eo=0;eo<Nelements;eo+=p_NblockS;@outer(0)){
  @shared dfloat s_rflux[p_NblockS][p_NfacesNfp];
  @exclusive dlong r_e, element;        // register-only variable per thread
  for(int es=0;es<p_NblockS;++es;@inner(1)){
    for(int n=0;n<p_maxNodes;++n;@inner(0)){
      ...
    }
  }
}
```

### Qualifiers
- `@restrict` on all pointer arguments (aliasing never assumed)
- `const` on all read-only pointer arguments
- `@shared` for intra-block staging arrays
- `@exclusive` for per-thread register variables visible across `@inner` loops
- `@barrier()` between phases that write then read global memory from within the same `@outer` iteration (used in `acousticsSurfaceHex3D.okl` between face batches)

### Compile-time Constants (kernel defines)
All tuning and mesh parameters are injected as preprocessor defines through `kernelInfo["defines/p_XXX"]`. Acoustics-specific defines:

| Define | Source | Meaning |
|---|---|---|
| `p_Nfields` | `Nfields` (3 or 4) | Number of conserved fields |
| `p_half` | `0.5` | Riemann flux coefficient |
| `p_maxNodes` | `max(Np, Nfp*Nfaces)` | Inner loop bound for surface kernels |
| `p_NblockV` | computed | Volume kernel block size |
| `p_NblockS` | computed | Surface kernel block size |
| `p_Lambda2` | `0.5` | Penalty parameter |

Mesh geometry constants (`p_Np`, `p_Nq`, `p_Nfp`, `p_Nfaces`, `p_Nvgeo`, `p_Nsgeo`, `p_RXID`, `p_NXID`, etc.) come from `mesh.props` propagated into `kernelInfo`.

### Two-Phase Volume Kernel Pattern
All volume kernels follow a strict two-phase pattern with a synchronization barrier between phases:
1. **Phase 1** (`@inner`): Load `q`, compute physical-space flux components `f`, `g`, `h`, transform to reference-space and store in `@shared` arrays `s_F`, `s_G`, `s_H`.
2. **Phase 2** (`@inner`, same outer block): Perform differentiation matrix multiply using `s_F/G/H`, accumulate into local `rhsq0`, `rhsq1`, ..., write to global `rhsq`.

### Riemann Solver Helper Functions
Device-side `upwind(nx, ny, [nz], rM, uM, vM, [wM], rP, uP, vP, [wP], *rflux, *uflux, *vflux, [*wflux])` is defined as a plain (non-`@kernel`) device function at file scope in each surface OKL file. Results are returned via pointer arguments.

For Hex3D, surface terms are additionally factored into a `surfaceTerms(...)` helper that writes directly to global `rhsq`, requiring `@barrier()` between face batches.

### Boundary Condition Macros
BCs are defined in data header files (`data/acousticsGaussian2D.h`, `data/acousticsGaussian3D.h`) and `#include`d into kernels via `kernelInfo["includes"] += dataFileName`. They are implemented as C preprocessor macros:

```c
#define acousticsDirichletConditions2D(bc, t, x, y, nx, ny, rM, uM, vM, rB, uB, vB) \
{ if(bc==2){ ... } else if(bc==1){ ... } }
```

BC codes: `1` = reflective wall (velocity reflected, pressure unchanged), `2` = outflow/absorbing (pressure negated, velocity unchanged).

---

## C++ Code Style

### Error Handling
- `LIBP_ABORT(msg, condition)` macro for fatal errors with a condition guard (evaluated at runtime). Example: `LIBP_ABORT("Usage: ...", argc!=2)`.
- `LIBP_FORCE_ABORT(msg)` for unconditional fatal errors (used in settings parser for unknown keys).
- No exceptions; no `try/catch` blocks.

### Memory Management
- Host allocations: `memory<T> buf; buf.malloc(N);` — RAII wrapper, no raw `new/delete`.
- Device allocations: `deviceMemory<T> o_buf = platform.malloc<T>(host_buf);` — allocates and optionally copies from host.
- Scratch memory: `memory<dfloat> scratch(2*Nscratch);` — constructor-allocated, automatic lifetime.

### MPI Scope
`main()` wraps all work in a block scope `{ ... }` to force destructors before `Comm::Finalize()`. MPI communicator is duplicated: `comm_t comm(Comm::World().Dup())`.

### Output
`printf(...)` used for solution norms and progress (not `std::cout`). `std::cout` used in settings `report()`. VTK output via `fprintf(fp, ...)` to `.vtu` files, with file naming pattern `{name}_{rank:04d}_{frame:04d}.vtu`.

### Standards
- C++17 (`-std=c++17` in `make.top`)
- Compiled with `mpic++` (GNU or Intel toolchain, default GNU)
- Optimization: `-O3 -fopenmp -mavx2 -ftree-vectorize -march=native` (release)
- Debug: `-O0 -g -Wall -Wshadow`

### No Inline Comments on `//` Lines Inside Macros
Macro definitions in data headers use `\`-continuation lines without inline comments (comments would break the macro).

---

## Common Preprocessor Patterns

```cpp
// Kernel suffix selection based on element type
std::string suffix;
if(mesh.elementType==Mesh::TRIANGLES)     suffix = "Tri2D";
if(mesh.elementType==Mesh::QUADRILATERALS) suffix = "Quad2D";
if(mesh.elementType==Mesh::TETRAHEDRA)    suffix = "Tet3D";
if(mesh.elementType==Mesh::HEXAHEDRA)     suffix = "Hex3D";

// Kernel file naming
std::string fileName = oklFilePrefix + "acousticsVolume" + suffix + oklFileSuffix;
std::string kernelName = "acousticsVolume" + suffix;
```

```cpp
// CUDA vs. CPU block-size selection
int blockMax = 256;
if (platform.device.mode() == "CUDA") blockMax = 512;
```

```cpp
// Dimension-dependent field count
Nfields = (mesh.dim==3) ? 4 : 3;
```

The macro `DACOUSTICS` expands to the absolute path of the acoustics solver directory, injected via `-DLIBP_DIR='"..."'` compile define:
```cpp
#define DACOUSTICS LIBP_DIR"/solvers/acoustics/"
```

---

## Module / File Organization

Each logical operation has its own `.cpp` file under `src/`:
- `acousticsSetup.cpp` — `Setup()`: allocates state, builds kernels, registers settings
- `acousticsRun.cpp` — `Run()`: sets initial conditions, computes time step, runs time integrator
- `acousticsStep.cpp` — `rhsf()` and `MaxWaveSpeed()`: the per-step ODE right-hand side
- `acousticsReport.cpp` — `Report()`: periodic norm output and optional VTK write
- `acousticsPlotFields.cpp` — `PlotFields()`: VTK ASCII output
- `acousticsSettings.cpp` — `acousticsSettings_t`: setting defaults and file parser

---

*Convention analysis: 2026-04-22*
