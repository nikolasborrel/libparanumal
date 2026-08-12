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
#include <algorithm>
#include <cmath>
#include <cstdio>

// LR vectorfit file format (matches DTU libparanumal-dtu convention):
//   line 1: NLRNpoles  NLRNRealPoles  NLRNImagPoles
//   then NLRNRealPoles  A      values (real residues)
//   then NLRNImagPoles  B      values (imag residues, real part)
//   then NLRNImagPoles  C      values (imag residues, imag part)
//   then NLRNRealPoles  lambda values (real pole decay rates, positive = decaying)
//   then NLRNImagPoles  alpha  values (imag pole real parts)
//   then NLRNImagPoles  beta   values (imag pole imag parts)
//   then 1              Y_inf  value  (direct / instantaneous admittance)
//
// LR array storage layout (size = 1 + 2*NReal + 4*NImag):
//   [p_LRA      .. ]  NReal A
//   [p_LRB      .. ]  NImag B
//   [p_LRC      .. ]  NImag C
//   [p_LRLambda .. ]  NReal lambda
//   [p_LRAlpha  .. ]  NImag alpha
//   [p_LRBeta   .. ]  NImag beta
//   [p_LRYinf      ]  Y_inf

void acoustics_t::SetupLRBC(properties_t& kernelInfo) {

  NLRPoints = 0;

  // Count LR faces before anything can return early. Without a fit the surface
  // kernel yields vn=0 on those faces, which is a rigid wall, so a mesh asking
  // for LR walls would otherwise be simulated silently as perfectly reflective.
  // Collective, so every rank aborts together.
  dlong NLRFaces = 0;
  for (dlong e = 0; e < mesh.Nelements; ++e)
    for (int f = 0; f < mesh.Nfaces; ++f)
      if (mesh.EToB[f + mesh.Nfaces*e] == 3) ++NLRFaces;
  mesh.comm.Allreduce(NLRFaces, Comm::Sum);

  std::string lrFile;
  if (settings.hasSetting("LR VECTORFIT FILE"))
    settings.getSetting("LR VECTORFIT FILE", lrFile);

  LIBP_ABORT("Mesh has BC==3 (locally-reacting) faces but no LR VECTORFIT FILE",
             NLRFaces > 0 && lrFile.empty());

  if (lrFile.empty()) return;

  // ------------------------------------------------------------------ //
  //  Read vectorfit data
  // ------------------------------------------------------------------ //
  FILE* fp = fopen(lrFile.c_str(), "r");
  LIBP_ABORT("LR VECTORFIT FILE not found: " << lrFile
             << ", and the mesh has BC==3 faces", !fp && NLRFaces > 0);
  if (!fp) {
    if (mesh.rank == 0)
      printf("WARNING: LR vectorfit file not found: %s — skipping LR BCs\n",
             lrFile.c_str());
    return;
  }

  fscanf(fp, "%lld %lld %lld",
         (long long*)&LRNpoles, (long long*)&LRNRealPoles, (long long*)&LRNImagPoles);

  // the surface kernel walks the complex pairs as pi += 2 reading acc[base+pi+1],
  // so an inconsistent header runs off the end of the accumulator row
  LIBP_ABORT("LR vectorfit header inconsistent: Npoles (" << LRNpoles
             << ") must equal NRealPoles + 2*NImagPoles ("
             << LRNRealPoles + 2*LRNImagPoles << ") in " << lrFile,
             LRNpoles != LRNRealPoles + 2*LRNImagPoles);

  const dlong LRsize = 1 + 2*LRNRealPoles + 4*LRNImagPoles;
  LR.malloc(LRsize, 0.0);
  LRInfo.malloc(3);
  LRInfo[0] = LRNpoles;
  LRInfo[1] = LRNRealPoles;
  LRInfo[2] = LRNImagPoles;

  // p_LRA = 0
  for (dlong i = 0; i < LRNRealPoles; ++i)
    fscanf(fp, "%lf", &LR[i]);
  dlong off = LRNRealPoles; // p_LRB
  for (dlong i = 0; i < LRNImagPoles; ++i)
    fscanf(fp, "%lf", &LR[off + i]);
  off += LRNImagPoles; // p_LRC
  for (dlong i = 0; i < LRNImagPoles; ++i)
    fscanf(fp, "%lf", &LR[off + i]);
  off += LRNImagPoles; // p_LRLambda
  for (dlong i = 0; i < LRNRealPoles; ++i)
    fscanf(fp, "%lf", &LR[off + i]);
  off += LRNRealPoles; // p_LRAlpha
  for (dlong i = 0; i < LRNImagPoles; ++i)
    fscanf(fp, "%lf", &LR[off + i]);
  off += LRNImagPoles; // p_LRBeta
  for (dlong i = 0; i < LRNImagPoles; ++i)
    fscanf(fp, "%lf", &LR[off + i]);
  off += LRNImagPoles; // p_LRYinf
  fscanf(fp, "%lf", &LR[off]);
  fclose(fp);

  // ------------------------------------------------------------------ //
  //  Build mapAcc: face-node index → accumulator row (0..NLRPoints-1)
  //  -1 for non-LR face nodes.
  // ------------------------------------------------------------------ //
  const dlong mapSize = mesh.Nelements * mesh.Nfp * mesh.Nfaces;
  mapAcc.malloc(mapSize, (dlong)-1);

  dlong counter = 0;
  for (dlong e = 0; e < mesh.Nelements; ++e) {
    for (int f = 0; f < mesh.Nfaces; ++f) {
      const int bc = mesh.EToB[f + mesh.Nfaces*e];
      if (bc == 3) {
        for (int n = 0; n < mesh.Nfp; ++n) {
          const dlong id = e*mesh.Nfp*mesh.Nfaces + f*mesh.Nfp + n;
          mapAcc[id] = counter++;
        }
      }
    }
  }
  NLRPoints = counter;

  // Every rank builds the same kernels below, so the decision to set LR up has
  // to be global: a rank owning no BC==3 face still takes part.
  dlong NLRPointsGlobal = NLRPoints;
  mesh.comm.Allreduce(NLRPointsGlobal, Comm::Sum);

  if (NLRPointsGlobal == 0) {
    if (mesh.rank == 0)
      printf("WARNING: LR VECTORFIT FILE given but no BC==3 faces found — skipping LR\n");
    LRNpoles = 0;
    return;
  }

  // the accumulators are co-advanced by the fixed-step LSERK4 loop in Run()
  LIBP_ABORT("Locally-reacting boundaries require TIME INTEGRATOR LSERK4",
             !settings.compareSetting("TIME INTEGRATOR", "LSERK4"));

  // ------------------------------------------------------------------ //
  //  Allocate accumulator buffers  acc[NLRPoints * LRNpoles]
  // ------------------------------------------------------------------ //
  const dlong accSize = std::max(NLRPoints*LRNpoles, (dlong)1);
  acc   .malloc(accSize, 0.0);
  resacc.malloc(accSize, 0.0);

  o_LR     = platform.malloc<dfloat>(LR);
  o_LRInfo = platform.malloc<dlong> (LRInfo);
  o_mapAcc = platform.malloc<dlong> (mapAcc);
  o_acc    = platform.malloc<dfloat>(acc);
  o_resacc = platform.malloc<dfloat>(resacc);
  o_rhsacc = platform.malloc<dfloat>(accSize);

  // ------------------------------------------------------------------ //
  //  Kernel defines for surface kernel
  // ------------------------------------------------------------------ //
  const dlong pLRA      = 0;
  const dlong pLRB      = LRNRealPoles;
  const dlong pLRC      = LRNRealPoles   + LRNImagPoles;
  const dlong pLRLambda = LRNRealPoles   + 2*LRNImagPoles;
  const dlong pLRAlpha  = 2*LRNRealPoles + 2*LRNImagPoles;
  const dlong pLRBeta   = 2*LRNRealPoles + 3*LRNImagPoles;
  const dlong pLRYinf   = 2*LRNRealPoles + 4*LRNImagPoles;

  // fastest pole rate, used to check the accumulator ODE against the explicit
  // time step in Run()
  LRMaxPole = 0.0;
  for (dlong i = 0; i < LRNRealPoles; ++i)
    LRMaxPole = std::max(LRMaxPole, std::abs(LR[pLRLambda + i]));
  for (dlong i = 0; i < LRNImagPoles; ++i)
    LRMaxPole = std::max(LRMaxPole, std::hypot(LR[pLRAlpha + i], LR[pLRBeta + i]));

  kernelInfo["defines/p_LRLambda"] = (int)pLRLambda;
  kernelInfo["defines/p_LRAlpha"]  = (int)pLRAlpha;
  kernelInfo["defines/p_LRBeta"]   = (int)pLRBeta;
  kernelInfo["defines/p_LRA"]      = (int)pLRA;
  kernelInfo["defines/p_LRB"]      = (int)pLRB;
  kernelInfo["defines/p_LRC"]      = (int)pLRC;
  kernelInfo["defines/p_LRYinf"]   = (int)pLRYinf;

  // ------------------------------------------------------------------ //
  //  Build LR accumulator update kernel
  // ------------------------------------------------------------------ //
  properties_t lrInfo = mesh.props;
  const int blockSize = 256;
  lrInfo["defines/p_blockSize"] = blockSize;
  updateKernelLR = platform.buildKernel(
      DACOUSTICS "okl/acousticsUpdateLR.okl",
      "acousticsUpdateLRAcc",
      lrInfo);

  if (mesh.rank == 0)
    printf("LR BCs: %lld boundary points, %lld poles (%lld real, %lld imag pairs)\n",
           (long long)NLRPointsGlobal, (long long)LRNpoles,
           (long long)LRNRealPoles, (long long)LRNImagPoles);
}
