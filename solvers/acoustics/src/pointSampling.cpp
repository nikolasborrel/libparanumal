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

#include "pointSampling.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace pointSampling {

namespace {

constexpr dfloat kRefTol = 1.0e-4;

bool insideReference(int dim, bool simplex, dfloat r, dfloat s, dfloat t) {
  const dfloat lo = -1.0 - kRefTol;
  if (simplex) {
    if (dim == 2)
      return (r >= lo) && (s >= lo) && (r + s <= -1.0 + kRefTol);
    return (r >= lo) && (s >= lo) && (t >= lo) &&
           (r + s + t <= -1.0 + kRefTol);
  }
  const dfloat hi = 1.0 + kRefTol;
  if (dim == 2)
    return (r >= lo && r <= hi) && (s >= lo && s <= hi);
  return (r >= lo && r <= hi) && (s >= lo && s <= hi) && (t >= lo && t <= hi);
}

// Reference vertex 0 sits at (-1,-1[,-1]) and vertex k>0 advances reference
// axis k-1 to +1, so xi_d = (ref_d + 1)/2 and X(xi) = X0 + sum_d xi_d (X_d+1 - X0).
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
    rOut = -1.0 + 2.0 * (( a11 * bx - a01 * by) / det);
    sOut = -1.0 + 2.0 * ((-a10 * bx + a00 * by) / det);
    tOut = -1.0;
    return insideReference(2, true, rOut, sOut, tOut);
  }

  const dfloat z0 = mesh.EZ[v + 0];
  const dfloat J00 = mesh.EX[v+1]-x0, J01 = mesh.EX[v+2]-x0, J02 = mesh.EX[v+3]-x0;
  const dfloat J10 = mesh.EY[v+1]-y0, J11 = mesh.EY[v+2]-y0, J12 = mesh.EY[v+3]-y0;
  const dfloat J20 = mesh.EZ[v+1]-z0, J21 = mesh.EZ[v+2]-z0, J22 = mesh.EZ[v+3]-z0;

  const dfloat det = J00 * (J11 * J22 - J12 * J21)
                   - J01 * (J10 * J22 - J12 * J20)
                   + J02 * (J10 * J21 - J11 * J20);
  if (std::fabs(det) < 1e-30) return false;

  const dfloat bx = px - x0, by = py - y0, bz = pz - z0;
  const dfloat xi0 = ( bx*(J11*J22-J12*J21) - J01*(by*J22-bz*J12) + J02*(by*J21-bz*J11) ) / det;
  const dfloat xi1 = ( J00*(by*J22-bz*J12) - bx*(J10*J22-J12*J20) + J02*(J10*bz-by*J20) ) / det;
  const dfloat xi2 = ( J00*(J11*bz-by*J21) - J01*(J10*bz-by*J20) + bx*(J10*J21-J11*J20) ) / det;

  rOut = -1.0 + 2.0 * xi0;
  sOut = -1.0 + 2.0 * xi1;
  tOut = -1.0 + 2.0 * xi2;
  return insideReference(3, true, rOut, sOut, tOut);
}

// signR/S/T hold each vertex's reference coordinates, giving the multilinear
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

  dfloat r = 0.0, s = 0.0, t = 0.0;   // start at the element centre
  for (int iter = 0; iter < 20; ++iter) {
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
      LIBP_FORCE_ABORT("pointSampling: unsupported element type");
  }
}

} // namespace

void locate(mesh_t& mesh,
            const memory<dfloat>& xyz,
            dlong npts,
            memory<dlong>& elements,
            memory<dfloat>& weights) {
  const int Np    = mesh.Np;
  const bool is3D = (mesh.dim == 3);
  const bool simplex = (mesh.elementType == Mesh::TRIANGLES ||
                        mesh.elementType == Mesh::TETRAHEDRA);

  elements.malloc(npts, (dlong)-1);
  weights.malloc(npts * Np, 0.0);

  std::vector<dfloat> signR, signS, signT;
  if (!simplex) {
    signR.resize(mesh.Nverts); signS.resize(mesh.Nverts); signT.resize(mesh.Nverts);
    for (int v = 0; v < mesh.Nverts; ++v) {
      const int n = mesh.vertexNodes[v];
      signR[v] = mesh.r[n];
      signS[v] = mesh.s[n];
      signT[v] = is3D ? mesh.t[n] : 0.0;
    }
  }

  // vertex bounding boxes, a cheap reject before the per-element inversion
  std::vector<dfloat> bxlo(mesh.Nelements), bxhi(mesh.Nelements);
  std::vector<dfloat> bylo(mesh.Nelements), byhi(mesh.Nelements);
  std::vector<dfloat> bzlo(mesh.Nelements), bzhi(mesh.Nelements);
  for (dlong e = 0; e < mesh.Nelements; ++e) {
    dfloat xlo= 1e300, ylo= 1e300, zlo= 1e300;
    dfloat xhi=-1e300, yhi=-1e300, zhi=-1e300;
    for (int v = 0; v < mesh.Nverts; ++v) {
      const dlong id = e * mesh.Nverts + v;
      xlo = std::min(xlo, mesh.EX[id]); xhi = std::max(xhi, mesh.EX[id]);
      ylo = std::min(ylo, mesh.EY[id]); yhi = std::max(yhi, mesh.EY[id]);
      if (is3D) { zlo = std::min(zlo, mesh.EZ[id]); zhi = std::max(zhi, mesh.EZ[id]); }
    }
    const dfloat pad = 1e-6 * (std::fabs(xhi-xlo) + std::fabs(yhi-ylo) + 1.0);
    bxlo[e]=xlo-pad; bxhi[e]=xhi+pad;
    bylo[e]=ylo-pad; byhi[e]=yhi+pad;
    bzlo[e]=zlo-pad; bzhi[e]=zhi+pad;
  }

  memory<dfloat> row;
  for (dlong i = 0; i < npts; ++i) {
    const dfloat px = xyz[i*3+0];
    const dfloat py = xyz[i*3+1];
    const dfloat pz = is3D ? xyz[i*3+2] : 0.0;

    dfloat r = 0, s = 0, t = -1;
    for (dlong e = 0; e < mesh.Nelements; ++e) {
      if (px < bxlo[e] || px > bxhi[e] ||
          py < bylo[e] || py > byhi[e] ||
          (is3D && (pz < bzlo[e] || pz > bzhi[e]))) continue;

      const bool in = simplex
        ? locateSimplex(mesh, e, px, py, pz, r, s, t)
        : locateTensor(mesh, e, signR, signS, signT, px, py, pz, r, s, t);
      if (!in) continue;

      interpolationRow(mesh, r, s, t, row);
      elements[i] = e;
      for (int n = 0; n < Np; ++n) weights[i*Np + n] = row[n];
      break;
    }
  }
}

} // namespace pointSampling
