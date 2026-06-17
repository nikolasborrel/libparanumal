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

#include "acoustics.hpp"
#include <cstdio>

// Brute-force point-in-tet search using half-space tests. Tet-only.
static void findReceiverElements(acoustics_t& ac) {
  mesh_t& mesh = ac.mesh;

  // face vertex ordering: [face_vertex_0, face_vertex_1, face_vertex_2, opposite_vertex]
  const int fv[4][4] = {{0,1,2,3},{0,1,3,2},{1,2,3,0},{2,0,3,1}};

  for (dlong k = 0; k < ac.NReceivers; k++) {
    const dfloat px = ac.recvXYZ[k*3+0];
    const dfloat py = ac.recvXYZ[k*3+1];
    const dfloat pz = ac.recvXYZ[k*3+2];

    for (dlong i = 0; i < mesh.Nelements; i++) {
      bool inside = true;
      for (int j = 0; j < mesh.Nfaces; j++) {
        const dfloat bx = mesh.EX[i*mesh.Nverts+fv[j][0]];
        const dfloat by = mesh.EY[i*mesh.Nverts+fv[j][0]];
        const dfloat bz = mesh.EZ[i*mesh.Nverts+fv[j][0]];

        const dfloat rx = mesh.EX[i*mesh.Nverts+fv[j][1]] - bx;
        const dfloat ry = mesh.EY[i*mesh.Nverts+fv[j][1]] - by;
        const dfloat rz = mesh.EZ[i*mesh.Nverts+fv[j][1]] - bz;

        const dfloat sx = mesh.EX[i*mesh.Nverts+fv[j][2]] - bx;
        const dfloat sy = mesh.EY[i*mesh.Nverts+fv[j][2]] - by;
        const dfloat sz = mesh.EZ[i*mesh.Nverts+fv[j][2]] - bz;

        // face normal
        const dfloat nx = ry*sz - rz*sy;
        const dfloat ny = rz*sx - rx*sz;
        const dfloat nz = rx*sy - ry*sx;

        // opposite vertex and receiver, both relative to face base
        const dfloat dx = mesh.EX[i*mesh.Nverts+fv[j][3]] - bx;
        const dfloat dy = mesh.EY[i*mesh.Nverts+fv[j][3]] - by;
        const dfloat dz = mesh.EZ[i*mesh.Nverts+fv[j][3]] - bz;

        const dfloat signedRecv  = nx*(px-bx) + ny*(py-by) + nz*(pz-bz);
        const dfloat signedOther = nx*dx + ny*dy + nz*dz;

        const dfloat absRecv = (signedRecv >= 0.0) ? signedRecv : -signedRecv;
        const int sideRecv   = (signedRecv  > 0.0) ? 1 : -1;
        const int sideOther  = (signedOther > 0.0) ? 1 : -1;

        if (sideRecv != sideOther && absRecv > 1.0e-15) {
          inside = false;
          break;
        }
      }

      if (inside) {
        ac.recvElements[k] = i;
        ac.recvElementsIdx[ac.NReceiversLocal] = k;
        ac.NReceiversLocal++;
        break;
      }
    }
  }
}

// Compute Lagrange interpolation weights at a receiver using InterpolationMatrixTet3D.
// Converts physical receiver coordinate to reference (r,s,t) via Cramer's rule on the
// linear geometric map, then calls the mesh static method.
static void buildInterpolationOperators(acoustics_t& ac) {
  mesh_t& mesh = ac.mesh;
  const int Np = mesh.Np;

  memory<dfloat> intpol(ac.NReceiversLocal * Np, 0.0);

  for (int iRecv = 0; iRecv < ac.NReceiversLocal; iRecv++) {
    const dlong rIdx = ac.recvElementsIdx[iRecv];
    const dlong ele  = ac.recvElements[rIdx];

    const dfloat x0 = mesh.EX[ele*mesh.Nverts+0], y0 = mesh.EY[ele*mesh.Nverts+0], z0 = mesh.EZ[ele*mesh.Nverts+0];
    const dfloat x1 = mesh.EX[ele*mesh.Nverts+1], y1 = mesh.EY[ele*mesh.Nverts+1], z1 = mesh.EZ[ele*mesh.Nverts+1];
    const dfloat x2 = mesh.EX[ele*mesh.Nverts+2], y2 = mesh.EY[ele*mesh.Nverts+2], z2 = mesh.EZ[ele*mesh.Nverts+2];
    const dfloat x3 = mesh.EX[ele*mesh.Nverts+3], y3 = mesh.EY[ele*mesh.Nverts+3], z3 = mesh.EZ[ele*mesh.Nverts+3];

    const dfloat px = ac.recvXYZ[rIdx*3+0];
    const dfloat py = ac.recvXYZ[rIdx*3+1];
    const dfloat pz = ac.recvXYZ[rIdx*3+2];

    // Jacobian of map x(xi) = x0 + J * xi,  xi = [(r+1)/2, (s+1)/2, (t+1)/2]
    const dfloat J00=x1-x0, J01=x2-x0, J02=x3-x0;
    const dfloat J10=y1-y0, J11=y2-y0, J12=y3-y0;
    const dfloat J20=z1-z0, J21=z2-z0, J22=z3-z0;

    const dfloat detJ = J00*(J11*J22-J12*J21)
                      - J01*(J10*J22-J12*J20)
                      + J02*(J10*J21-J11*J20);

    // rhs = p - x0
    const dfloat bx = px - x0, by = py - y0, bz = pz - z0;

    // Cramer's rule: xi = J^{-1} b
    const dfloat xi0 = ( bx*(J11*J22-J12*J21) - J01*(by*J22-bz*J12) + J02*(by*J21-bz*J11) ) / detJ;
    const dfloat xi1 = ( J00*(by*J22-bz*J12) - bx*(J10*J22-J12*J20) + J02*(J10*bz-by*J20) ) / detJ;
    const dfloat xi2 = ( J00*(J11*bz-by*J21) - J01*(J10*bz-by*J20) + bx*(J10*J21-J11*J20) ) / detJ;

    // Reference element coordinates
    memory<dfloat> rOut(1), sOut(1), tOut(1);
    rOut[0] = -1.0 + 2.0*xi0;
    sOut[0] = -1.0 + 2.0*xi1;
    tOut[0] = -1.0 + 2.0*xi2;

    memory<dfloat> row;
    mesh_t::InterpolationMatrixTet3D(mesh.N, mesh.r, mesh.s, mesh.t, rOut, sOut, tOut, row);

    for (int i = 0; i < Np; i++)
      intpol[iRecv*Np + i] = row[i];
  }

  ac.o_recvIP = ac.platform.malloc<dfloat>(intpol);
}

void acoustics_t::SetupReceivers() {
  NReceivers = 0;
  NReceiversLocal = 0;
  NRecvSamples = 0;
  recvSampleIdx = 0;

  // Receivers are optional; skip if setting absent or empty
  std::string recvFile;
  if (!settings.hasSetting("RECEIVER FILE"))           return;
  settings.getSetting("RECEIVER FILE", recvFile);
  if (recvFile.empty())                                return;

  FILE* fp = fopen(recvFile.c_str(), "r");
  if (!fp) {
    if (mesh.rank == 0)
      printf("WARNING: receiver file not found: %s — skipping receivers\n",
             recvFile.c_str());
    return;
  }

  fscanf(fp, "%d", &NReceivers);
  recvXYZ.malloc(NReceivers * 3);
  for (int i = 0; i < NReceivers; i++)
    fscanf(fp, "%lf %lf %lf",
           &recvXYZ[i*3+0], &recvXYZ[i*3+1], &recvXYZ[i*3+2]);
  fclose(fp);

  recvElements.malloc(NReceivers, (dlong)-1);
  recvElementsIdx.malloc(NReceivers, (dlong)0);

  findReceiverElements(*this);

  if (NReceiversLocal == 0) return;

  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  // NRecvSamples = one sample per output step
  dfloat outputInterval;
  settings.getSetting("OUTPUT INTERVAL", outputInterval);

  const dfloat duration = finalTime - startTime;
  NRecvSamples = (dlong)(duration / outputInterval) + 2; // +2 for t=0 and rounding

  // Sample rate in Hz (outputInterval is physical seconds).
  sampleRateOut = (int)round(1.0 / outputInterval);

  qRecv.malloc(NReceiversLocal * NRecvSamples, 0.0);
  o_qRecv = platform.malloc<dfloat>(qRecv);

  o_recvElements    = platform.malloc<dlong>(recvElements);
  o_recvElementsIdx = platform.malloc<dlong>(recvElementsIdx);

  buildInterpolationOperators(*this);

  // Build the receiver interpolation kernel (needs NRecvSamples to be known)
  properties_t recvInfo = mesh.props;
  const int blockSize = 256;
  recvInfo["defines/p_blockSize"]    = blockSize;
  recvInfo["defines/p_Np"]           = mesh.Np;
  recvInfo["defines/p_Nfields"]      = Nfields;
  recvInfo["defines/p_NRecvSamples"] = (int)NRecvSamples;

  receiverKernel = platform.buildKernel(
      DACOUSTICS "okl/acousticsReceiverKernel.okl",
      "acousticsReceiverInterpolation",
      recvInfo);
}
