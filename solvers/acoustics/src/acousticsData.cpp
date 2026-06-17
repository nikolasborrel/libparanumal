/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

// ===========================================================================
// ML data-generation spatial sampler.
//
// Builds a uniform points-per-wavelength (PPW) Cartesian grid over the mesh
// bounding box, jitters each sample point's coordinates (seeded, uniform in
// +/- DATA JITTER * dx per axis), locates the element containing each jittered
// point and builds a Lagrange interpolation operator there. During the run the
// pressure field is interpolated onto these (jittered) points on the data-gen
// temporal cadence and persisted with the jittered coordinates, so a consuming
// ML model trains on varied locations rather than a fixed grid.
//
// This path is orthogonal to the visualization output: it does not touch
// Report()/OUTPUT INTERVAL or the VTU/XDMF writers.
// ===========================================================================

#include "acoustics.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

// ---------------------------------------------------------------------------
// Point location: given a physical point, find the local element that contains
// it and the corresponding reference coordinates (r,s,t).
//
//  - Simplices (Tri2D, Tet3D): the geometric map is affine, so the reference
//    coordinates follow from a single linear solve against the element vertices
//    (standard libParanumal reference-vertex ordering).
//  - Tensor elements (Quad2D, Hex3D): the vertex map is (bi/tri)linear; invert
//    it with a few Newton iterations. Corner signs are read from the reference
//    nodes at the vertices so we do not hard-code a vertex ordering.
// ---------------------------------------------------------------------------

namespace {

// Reference-coordinate tolerance for "inside element" acceptance.
constexpr dfloat kRefTol = 1.0e-4;

bool insideReference(int dim, bool simplex,
                     dfloat r, dfloat s, dfloat t) {
  const dfloat lo = -1.0 - kRefTol;
  if (simplex) {
    if (dim == 2)
      return (r >= lo) && (s >= lo) && (r + s <= -1.0 + kRefTol);
    return (r >= lo) && (s >= lo) && (t >= lo) &&
           (r + s + t <= -1.0 + kRefTol);
  } else {
    const dfloat hi = 1.0 + kRefTol;
    if (dim == 2)
      return (r >= lo && r <= hi) && (s >= lo && s <= hi);
    return (r >= lo && r <= hi) && (s >= lo && s <= hi) && (t >= lo && t <= hi);
  }
}

// Affine inverse for simplices. Vertices are X_v = (EX,EY[,EZ]) in standard
// order; reference vertex 0 is at (-1,-1[,-1]) and vertex k (k>0) advances
// reference axis (k-1) to +1, so xi_d = (ref_d + 1)/2 and
//   X(xi) = X0 + sum_d xi_d * (X_d+1 - X0).
bool locateSimplex(mesh_t& mesh, dlong e,
                   dfloat px, dfloat py, dfloat pz,
                   dfloat& rOut, dfloat& sOut, dfloat& tOut) {
  const dlong v = e * mesh.Nverts;
  const dfloat x0 = mesh.EX[v + 0], y0 = mesh.EY[v + 0];

  if (mesh.dim == 2) {
    const dfloat a00 = mesh.EX[v + 1] - x0, a01 = mesh.EX[v + 2] - x0;
    const dfloat a10 = mesh.EY[v + 1] - y0, a11 = mesh.EY[v + 2] - y0;
    const dfloat det = a00 * a11 - a01 * a10;
    if (std::fabs(det) < 1e-30) return false;
    const dfloat bx = px - x0, by = py - y0;
    const dfloat xi0 = ( a11 * bx - a01 * by) / det;
    const dfloat xi1 = (-a10 * bx + a00 * by) / det;
    rOut = -1.0 + 2.0 * xi0;
    sOut = -1.0 + 2.0 * xi1;
    tOut = -1.0;
    return insideReference(2, true, rOut, sOut, tOut);
  }

  const dfloat z0 = mesh.EZ[v + 0];
  const dfloat J00 = mesh.EX[v + 1] - x0, J01 = mesh.EX[v + 2] - x0, J02 = mesh.EX[v + 3] - x0;
  const dfloat J10 = mesh.EY[v + 1] - y0, J11 = mesh.EY[v + 2] - y0, J12 = mesh.EY[v + 3] - y0;
  const dfloat J20 = mesh.EZ[v + 1] - z0, J21 = mesh.EZ[v + 2] - z0, J22 = mesh.EZ[v + 3] - z0;

  const dfloat det = J00 * (J11 * J22 - J12 * J21)
                   - J01 * (J10 * J22 - J12 * J20)
                   + J02 * (J10 * J21 - J11 * J20);
  if (std::fabs(det) < 1e-30) return false;

  const dfloat bx = px - x0, by = py - y0, bz = pz - z0;
  const dfloat xi0 = ( bx * (J11 * J22 - J12 * J21) - J01 * (by * J22 - bz * J12) + J02 * (by * J21 - bz * J11) ) / det;
  const dfloat xi1 = ( J00 * (by * J22 - bz * J12) - bx * (J10 * J22 - J12 * J20) + J02 * (J10 * bz - by * J20) ) / det;
  const dfloat xi2 = ( J00 * (J11 * bz - by * J21) - J01 * (J10 * bz - by * J20) + bx * (J10 * J21 - J11 * J20) ) / det;

  rOut = -1.0 + 2.0 * xi0;
  sOut = -1.0 + 2.0 * xi1;
  tOut = -1.0 + 2.0 * xi2;
  return insideReference(3, true, rOut, sOut, tOut);
}

// Newton inverse for tensor-product elements (Quad2D/Hex3D). `signR/S/T` are the
// reference coordinates (+/-1) of each vertex, used to form the multilinear
// corner shape functions N_v(xi) = prod_d 0.5*(1 + signD[v]*xi_d).
bool locateTensor(mesh_t& mesh, dlong e,
                  const std::vector<dfloat>& signR,
                  const std::vector<dfloat>& signS,
                  const std::vector<dfloat>& signT,
                  dfloat px, dfloat py, dfloat pz,
                  dfloat& rOut, dfloat& sOut, dfloat& tOut) {
  const int dim = mesh.dim;
  const int nv  = mesh.Nverts;
  const dlong vb = e * nv;

  dfloat r = 0.0, s = 0.0, t = 0.0;  // start at element centre
  for (int iter = 0; iter < 20; ++iter) {
    // residual F = X(xi) - p, and Jacobian dX/dxi
    dfloat Fx = -px, Fy = -py, Fz = -pz;
    dfloat Jx[3] = {0,0,0}, Jy[3] = {0,0,0}, Jz[3] = {0,0,0};
    for (int v = 0; v < nv; ++v) {
      const dfloat sr = signR[v], ss = signS[v], st = (dim == 3) ? signT[v] : 0.0;
      const dfloat ar = 0.5 * (1.0 + sr * r);
      const dfloat as = 0.5 * (1.0 + ss * s);
      const dfloat at = (dim == 3) ? 0.5 * (1.0 + st * t) : 1.0;
      const dfloat N  = ar * as * at;
      const dfloat Xv = mesh.EX[vb + v], Yv = mesh.EY[vb + v];
      const dfloat Zv = (dim == 3) ? mesh.EZ[vb + v] : 0.0;
      Fx += N * Xv; Fy += N * Yv; Fz += N * Zv;
      // dN/dr, dN/ds, dN/dt
      const dfloat dNr = 0.5 * sr * as * at;
      const dfloat dNs = 0.5 * ss * ar * at;
      Jx[0] += dNr * Xv; Jy[0] += dNr * Yv; Jz[0] += dNr * Zv;
      Jx[1] += dNs * Xv; Jy[1] += dNs * Yv; Jz[1] += dNs * Zv;
      if (dim == 3) {
        const dfloat dNt = 0.5 * st * ar * as;
        Jx[2] += dNt * Xv; Jy[2] += dNt * Yv; Jz[2] += dNt * Zv;
      }
    }

    dfloat dr, ds, dt = 0.0;
    if (dim == 2) {
      const dfloat det = Jx[0] * Jy[1] - Jx[1] * Jy[0];
      if (std::fabs(det) < 1e-30) return false;
      dr = ( Jy[1] * Fx - Jx[1] * Fy) / det;
      ds = (-Jy[0] * Fx + Jx[0] * Fy) / det;
    } else {
      const dfloat det = Jx[0]*(Jy[1]*Jz[2]-Jy[2]*Jz[1])
                       - Jx[1]*(Jy[0]*Jz[2]-Jy[2]*Jz[0])
                       + Jx[2]*(Jy[0]*Jz[1]-Jy[1]*Jz[0]);
      if (std::fabs(det) < 1e-30) return false;
      // solve J * dxi = F via Cramer's rule
      dr = ( Fx*(Jy[1]*Jz[2]-Jy[2]*Jz[1]) - Jx[1]*(Fy*Jz[2]-Fz*Jz[1]) + Jx[2]*(Fy*Jz[1]-Fz*Jy[1]) ) / det;
      ds = ( Jx[0]*(Fy*Jz[2]-Fz*Jz[1]) - Fx*(Jy[0]*Jz[2]-Jy[2]*Jz[0]) + Jx[2]*(Jy[0]*Fz-Fy*Jz[0]) ) / det;
      dt = ( Jx[0]*(Jy[1]*Fz-Fy*Jz[1]) - Jx[1]*(Jy[0]*Fz-Fy*Jz[0]) + Fx*(Jy[0]*Jz[1]-Jy[1]*Jz[0]) ) / det;
    }

    r -= dr; s -= ds; t -= dt;
    if (std::fabs(dr) + std::fabs(ds) + std::fabs(dt) < 1e-12) break;
  }

  rOut = r; sOut = s; tOut = (dim == 3) ? t : -1.0;
  return insideReference(dim, false, rOut, sOut, tOut);
}

// Build the InterpolationMatrix row (Np weights) for one reference coordinate.
void interpolationRow(mesh_t& mesh, dfloat r, dfloat s, dfloat t,
                      memory<dfloat>& row) {
  memory<dfloat> rOut(1), sOut(1), tOut(1);
  rOut[0] = r; sOut[0] = s; tOut[0] = t;
  switch (mesh.elementType) {
    case Mesh::TRIANGLES:
      mesh_t::InterpolationMatrixTri2D(mesh.N, mesh.r, mesh.s, rOut, sOut, row);
      break;
    case Mesh::QUADRILATERALS:
      mesh_t::InterpolationMatrixQuad2D(mesh.N, mesh.r, mesh.s, rOut, sOut, row);
      break;
    case Mesh::TETRAHEDRA:
      mesh_t::InterpolationMatrixTet3D(mesh.N, mesh.r, mesh.s, mesh.t, rOut, sOut, tOut, row);
      break;
    case Mesh::HEXAHEDRA:
      mesh_t::InterpolationMatrixHex3D(mesh.N, mesh.r, mesh.s, mesh.t, rOut, sOut, tOut, row);
      break;
    default:
      LIBP_FORCE_ABORT("Data-gen: unsupported element type");
  }
}

} // namespace

// ---------------------------------------------------------------------------
// SetupDataGen
// ---------------------------------------------------------------------------
void acoustics_t::SetupDataGen() {
  dataGenEnabled = false;
  NDataPoints = 0;
  NDataSamples = 0;
  dataSampleIdx = 0;

  // Data generation is enabled simply by referencing a data config from the main
  // setup via [DATA CONFIG FILE] — no separate on/off flag. Visualization output
  // is controlled independently in the main config (OUTPUT TO FILE / OUTPUT FORMAT).
  std::string dataConfig;
  if (settings.hasSetting("DATA CONFIG FILE"))
    settings.getSetting("DATA CONFIG FILE", dataConfig);
  if (dataConfig.empty()) return;
  dataGenEnabled = true;

  settings.getSetting("SPATIAL PPW DATA",  dataPPW);
  settings.getSetting("DATA JITTER",       dataJitter);
  settings.getSetting("DATA SEED",         dataSeed);

  // fmax is a shared source/mesh property (top-level [FMAX]), not a data-only knob.
  dataFmax = fmax;

  dfloat temporalPPW = 20.0;
  settings.getSetting("DATA TEMPORAL PPW", temporalPPW);

  dataSimID.clear();
  if (settings.hasSetting("DATA SIMULATION ID"))
    settings.getSetting("DATA SIMULATION ID", dataSimID);
  if (dataSimID.empty()) dataSimID = simulationID;

  LIBP_ABORT("Data-gen: FMAX must be > 0",              dataFmax <= 0.0);
  LIBP_ABORT("Data-gen: SPATIAL PPW DATA must be > 0",  dataPPW  <= 0.0);
  // Nyquist: at least 2 samples per period at fmax, otherwise the requested
  // configuration itself aliases regardless of the solver resolution.
  LIBP_ABORT("Data-gen: DATA TEMPORAL PPW must be >= 2 (Nyquist for fmax)",
             temporalPPW < 2.0);

  // Grid spacing from the spatial points-per-wavelength at the max frequency.
  dataDx = c / (dataFmax * dataPPW);
  // Temporal sampling interval from the temporal points-per-period.
  dataDt = 1.0 / (dataFmax * temporalPPW);

  // ---- global bounding box (reduce over ranks on local node coords) -------
  dfloat xmin =  1e300, ymin =  1e300, zmin =  1e300;
  dfloat xmax = -1e300, ymax = -1e300, zmax = -1e300;
  const dlong Nnode = mesh.Nelements * mesh.Np;
  for (dlong n = 0; n < Nnode; ++n) {
    xmin = std::min(xmin, mesh.x[n]); xmax = std::max(xmax, mesh.x[n]);
    ymin = std::min(ymin, mesh.y[n]); ymax = std::max(ymax, mesh.y[n]);
    if (mesh.dim == 3) { zmin = std::min(zmin, mesh.z[n]); zmax = std::max(zmax, mesh.z[n]); }
  }
  if (mesh.dim != 3) { zmin = 0.0; zmax = 0.0; }
  mesh.comm.Allreduce(xmin, Comm::Min); mesh.comm.Allreduce(xmax, Comm::Max);
  mesh.comm.Allreduce(ymin, Comm::Min); mesh.comm.Allreduce(ymax, Comm::Max);
  if (mesh.dim == 3) { mesh.comm.Allreduce(zmin, Comm::Min); mesh.comm.Allreduce(zmax, Comm::Max); }

  const int nx = std::max(1, (int)std::floor((xmax - xmin) / dataDx) + 1);
  const int ny = std::max(1, (int)std::floor((ymax - ymin) / dataDx) + 1);
  const int nz = (mesh.dim == 3) ? std::max(1, (int)std::floor((zmax - zmin) / dataDx) + 1) : 1;
  const long long nGrid = (long long)nx * ny * nz;

  Logf("Data-gen: PPW grid %d x %d x %d = %lld points, dx=%.4g m, "
       "dt=%.4g s (fmax=%.4g Hz, sppw=%.3g, tppw=%.3g, jitter=%.3g, seed=%d)\n",
       nx, ny, nz, nGrid, dataDx, dataDt, dataFmax, dataPPW, temporalPPW,
       dataJitter, dataSeed);
  if (nGrid > 5000000LL)
    Logf("Data-gen: WARNING large sample grid (%lld points) — check FMAX / SPATIAL PPW DATA\n", nGrid);

  // Corner signs for tensor-element point location (read from reference nodes).
  const bool simplex = (mesh.elementType == Mesh::TRIANGLES ||
                        mesh.elementType == Mesh::TETRAHEDRA);
  std::vector<dfloat> signR, signS, signT;
  if (!simplex) {
    signR.resize(mesh.Nverts); signS.resize(mesh.Nverts); signT.resize(mesh.Nverts);
    for (int v = 0; v < mesh.Nverts; ++v) {
      const int nodeId = mesh.vertexNodes[v];
      signR[v] = mesh.r[nodeId];
      signS[v] = mesh.s[nodeId];
      signT[v] = (mesh.dim == 3) ? mesh.t[nodeId] : 0.0;
    }
  }

  // Per-element vertex bounding boxes for a cheap point-location pre-filter.
  std::vector<dfloat> ebxmin(mesh.Nelements), ebxmax(mesh.Nelements);
  std::vector<dfloat> ebymin(mesh.Nelements), ebymax(mesh.Nelements);
  std::vector<dfloat> ebzmin(mesh.Nelements), ebzmax(mesh.Nelements);
  for (dlong e = 0; e < mesh.Nelements; ++e) {
    dfloat xlo= 1e300,ylo= 1e300,zlo= 1e300, xhi=-1e300,yhi=-1e300,zhi=-1e300;
    for (int v = 0; v < mesh.Nverts; ++v) {
      const dlong id = e * mesh.Nverts + v;
      xlo=std::min(xlo,mesh.EX[id]); xhi=std::max(xhi,mesh.EX[id]);
      ylo=std::min(ylo,mesh.EY[id]); yhi=std::max(yhi,mesh.EY[id]);
      if (mesh.dim==3){ zlo=std::min(zlo,mesh.EZ[id]); zhi=std::max(zhi,mesh.EZ[id]); }
    }
    const dfloat pad = 1e-6 * (std::fabs(xhi-xlo)+std::fabs(yhi-ylo)+1.0);
    ebxmin[e]=xlo-pad; ebxmax[e]=xhi+pad;
    ebymin[e]=ylo-pad; ebymax[e]=yhi+pad;
    ebzmin[e]=zlo-pad; ebzmax[e]=zhi+pad;
  }

  // Build a uniform PPW grid (optionally jittered) and DG interpolation operators.
  //   keepOutside=false → drop points outside the local domain (query/output grid).
  //   keepOutside=true  → keep all nx*ny*nz points (elem=-1, zero weights outside),
  //                       preserving the rectilinear shape (source/umesh grid).
  // Every rank generates the identical grid (same seed + traversal); for the query
  // grid each rank keeps only the points whose containing element is local.
  auto buildGrid = [&](dfloat ppw, dfloat jitterFrac, unsigned seed, bool keepOutside,
                       std::vector<dfloat>& outXYZ, std::vector<dlong>& outElem,
                       std::vector<dfloat>& outIP, int shape[3]) {
    const dfloat dx = c / (dataFmax * ppw);
    const int gnx = std::max(1, (int)std::floor((xmax - xmin) / dx) + 1);
    const int gny = std::max(1, (int)std::floor((ymax - ymin) / dx) + 1);
    const int gnz = (mesh.dim == 3) ? std::max(1, (int)std::floor((zmax - zmin) / dx) + 1) : 1;
    shape[0] = gnx; shape[1] = gny; shape[2] = gnz;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<dfloat> jit(-jitterFrac * dx, jitterFrac * dx);
    memory<dfloat> row;

    for (int k = 0; k < gnz; ++k) {
      for (int j = 0; j < gny; ++j) {
        for (int i = 0; i < gnx; ++i) {
          dfloat jx = 0.0, jy = 0.0, jz = 0.0;
          if (jitterFrac > 0.0) {       // jitterFrac==0 → exact rectilinear grid
            jx = jit(rng); jy = jit(rng);
            if (mesh.dim == 3) jz = jit(rng);
          }
          const dfloat px = xmin + i * dx + jx;
          const dfloat py = ymin + j * dx + jy;
          const dfloat pz = (mesh.dim == 3) ? (zmin + k * dx + jz) : 0.0;

          dlong  foundE = -1;
          dfloat rr = 0, ss = 0, tt = -1;
          for (dlong e = 0; e < mesh.Nelements; ++e) {
            if (px < ebxmin[e] || px > ebxmax[e] ||
                py < ebymin[e] || py > ebymax[e] ||
                (mesh.dim == 3 && (pz < ebzmin[e] || pz > ebzmax[e]))) continue;
            bool in = simplex
              ? locateSimplex(mesh, e, px, py, pz, rr, ss, tt)
              : locateTensor(mesh, e, signR, signS, signT, px, py, pz, rr, ss, tt);
            if (in) { foundE = e; break; }
          }

          if (foundE >= 0) {
            interpolationRow(mesh, rr, ss, tt, row);
            outXYZ.push_back(px); outXYZ.push_back(py); outXYZ.push_back(pz);
            outElem.push_back(foundE);
            for (int n = 0; n < mesh.Np; ++n) outIP.push_back(row[n]);
          } else if (keepOutside) {
            outXYZ.push_back(px); outXYZ.push_back(py); outXYZ.push_back(pz);
            outElem.push_back(-1);
            for (int n = 0; n < mesh.Np; ++n) outIP.push_back(0.0);
          }
        }
      }
    }
  };

  // ---- query / output grid (jittered; outside-domain points dropped) ------
  std::vector<dfloat> keepXYZ, keepIP;
  std::vector<dlong>  keepElem;
  int qshape[3];
  buildGrid(dataPPW, dataJitter, (unsigned)dataSeed, /*keepOutside=*/false,
            keepXYZ, keepElem, keepIP, qshape);

  NDataPoints = (dlong)keepElem.size();
  if (mesh.rank == 0)
    printf("Data-gen: located %lld output sample points inside the domain\n",
           (long long)NDataPoints);

  if (NDataPoints > 0) {
    dfloat startTime, finalTime;
    settings.getSetting("START TIME", startTime);
    settings.getSetting("FINAL TIME", finalTime);

    // The data grid is sampled every dataStride solver steps. The point of DATA
    // TEMPORAL PPW is to land as close as possible to the requested interval, so
    // round the ratio to NEAREST — except rounding up must never push the rate
    // below Nyquist for fmax (achievedDt <= 1/(2*fmax)), which would alias. So
    // cap the stride at the coarsest one that still satisfies Nyquist. If even
    // sampling every step (stride 1) is coarser than the Nyquist interval, the
    // solver dt itself cannot resolve fmax and we warn — aliasing is then
    // unavoidable without refining the mesh / lowering CFL or FMAX. NDataSamples
    // is sized from the achieved interval so sampling never overruns the buffer.
    const dfloat dtNyquist   = 1.0/(2.0*dataFmax);                  // coarsest non-aliasing interval
    const int    strideRound = std::max(1, (int)std::lround(dataDt   / dt));
    const int    strideNyq   = (int)std::floor(dtNyquist / dt);     // 0 if dt itself aliases
    dataStride = std::max(1, std::min(strideRound, std::max(1, strideNyq)));

    const dfloat achievedDt = dataStride * dt;
    NDataSamples = (dlong)((finalTime - startTime) / achievedDt) + 2;

    const dfloat configPPW   = 1.0/(dataFmax*dataDt);
    const dfloat achievedPPW = 1.0/(dataFmax*achievedDt);
    Logf("Data-gen temporal cadence: solver dt=%.4g s, requested dt=%.4g s "
         "(tppw=%.3g) -> every %d step(s) = %.4g s (tppw=%.3g)\n",
         dt, dataDt, configPPW, dataStride, achievedDt, achievedPPW);

    const bool rounded = std::fabs(achievedDt - dataDt) > 1e-9 * dataDt;
    if (strideNyq < 1)
      Logf("Data-gen ERROR: data dt=%.4g s could not be met (solver dt=%.4g s "
           "aliases fmax=%.4g Hz); using dt=%.4g s, tppw=%.3g < 2\n",
           dataDt, dt, dataFmax, achievedDt, achievedPPW);
    else if (rounded)
      Logf("WARNING: data dt=%.4g s could not be met, rounded to dt=%.4g s "
           "(tppw %.3g -> %.3g)\n", dataDt, achievedDt, configPPW, achievedPPW);

    // Pack host arrays + device buffers.
    dataXYZ.malloc(NDataPoints * 3);
    dataElements.malloc(NDataPoints);
    dataElementsIdx.malloc(NDataPoints);
    for (dlong p = 0; p < NDataPoints; ++p) {
      dataXYZ[p*3+0] = keepXYZ[p*3+0];
      dataXYZ[p*3+1] = keepXYZ[p*3+1];
      dataXYZ[p*3+2] = keepXYZ[p*3+2];
      dataElements[p]    = keepElem[p];
      dataElementsIdx[p] = p;   // identity (shared receiver kernel applies indirection)
    }

    memory<dfloat> ip(NDataPoints * mesh.Np);
    for (dlong n = 0; n < NDataPoints * mesh.Np; ++n) ip[n] = keepIP[n];

    o_dataIP          = platform.malloc<dfloat>(ip);
    o_dataElements    = platform.malloc<dlong>(dataElements);
    o_dataElementsIdx = platform.malloc<dlong>(dataElementsIdx);

    dataVals.malloc(NDataPoints * NDataSamples, 0.0);
    o_dataVals = platform.malloc<dfloat>(dataVals);
    dataTimes.clear();
    dataTimes.reserve(NDataSamples);

    // Reuse the receiver interpolation kernel (interpolates field 0 = pressure at
    // arbitrary points); a separate build pins p_NRecvSamples to NDataSamples.
    properties_t dataInfo = mesh.props;
    dataInfo["defines/p_blockSize"]    = 256;
    dataInfo["defines/p_Np"]           = mesh.Np;
    dataInfo["defines/p_Nfields"]      = Nfields;
    dataInfo["defines/p_NRecvSamples"] = (int)NDataSamples;
    dataKernel = platform.buildKernel(
        DACOUSTICS "okl/acousticsReceiverKernel.okl",
        "acousticsReceiverInterpolation",
        dataInfo);
  }

  // ---- rectilinear SOURCE grid (DeepONet branch input) --------------------
  // The initial condition sampled on a uniform PPW grid (kept full nx*ny*nz so
  // the rectilinear shape is preserved; points outside the domain are 0).
  settings.getSetting("DATA SOURCE PPW",    srcPPW);
  settings.getSetting("DATA SOURCE JITTER", srcJitter);
  if (srcPPW > 0.0) {
    std::vector<dfloat> sxyz, sip;
    std::vector<dlong>  selem;
    buildGrid(srcPPW, srcJitter, (unsigned)(dataSeed + 1u), /*keepOutside=*/true,
              sxyz, selem, sip, srcShape);
    NSrcPoints = (dlong)selem.size();   // = nx*ny*nz

    srcGridXYZ.malloc(NSrcPoints * 3);
    srcGridElem.malloc(NSrcPoints);
    srcGridIP.malloc(NSrcPoints * mesh.Np);
    srcGridVals.malloc(NSrcPoints, 0.0f);
    for (dlong p = 0; p < NSrcPoints; ++p) {
      srcGridXYZ[p*3+0] = sxyz[p*3+0];
      srcGridXYZ[p*3+1] = sxyz[p*3+1];
      srcGridXYZ[p*3+2] = sxyz[p*3+2];
      srcGridElem[p]    = selem[p];
    }
    for (dlong n = 0; n < NSrcPoints * mesh.Np; ++n) srcGridIP[n] = sip[n];

    if (mesh.rank == 0)
      printf("Data-gen: source grid %d x %d x %d = %lld points, dx=%.4g m (ppw=%.3g, jitter=%.3g)\n",
             srcShape[0], srcShape[1], srcShape[2], (long long)NSrcPoints,
             c / (dataFmax * srcPPW), srcPPW, srcJitter);
  }
}

// ---------------------------------------------------------------------------
// SampleSourceGrid — interpolate the initial condition onto the source grid.
// Call once at t=0 (o_q holds the IC). DG reinterpolation of pressure (field 0);
// points outside the (local) domain stay 0.
// NOTE: single-rank oriented (like the rest of the data-gen write path).
// ---------------------------------------------------------------------------
void acoustics_t::SampleSourceGrid() {
  if (!dataGenEnabled || NSrcPoints == 0) return;

  o_q.copyTo(q);                       // host copy; pressure is field 0
  const int Np = mesh.Np;
  for (dlong p = 0; p < NSrcPoints; ++p) {
    const dlong e = srcGridElem[p];
    if (e < 0) { srcGridVals[p] = 0.0f; continue; }   // outside the domain
    const dlong qoff = e * Np * Nfields;
    dfloat val = 0.0;
    for (int n = 0; n < Np; ++n) val += srcGridIP[p*Np + n] * q[qoff + n];
    srcGridVals[p] = (float)val;
  }
}

// ---------------------------------------------------------------------------
// SampleDataGrid — interpolate the pressure field into frame `frameIdx`.
// ---------------------------------------------------------------------------
void acoustics_t::SampleDataGrid(int frameIdx, dfloat time) {
  if (!dataGenEnabled || NDataPoints == 0) return;
  if (frameIdx >= NDataSamples) return;
  dataKernel(NDataPoints, (dlong)frameIdx,
             o_dataVals, o_dataElements, o_dataElementsIdx,
             o_dataIP, o_q);
  dataTimes.push_back(time);
  dataSampleIdx = frameIdx + 1;
}
