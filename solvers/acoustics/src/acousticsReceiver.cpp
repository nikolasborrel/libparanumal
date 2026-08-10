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
#include "pointSampling.hpp"
#include <cstdio>

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

  recvElementsIdx.malloc(NReceivers, (dlong)0);

  memory<dfloat> recvIP;
  pointSampling::locate(mesh, recvXYZ, NReceivers, recvElements, recvIP);
  for (dlong k = 0; k < NReceivers; ++k)
    if (recvElements[k] >= 0) recvElementsIdx[NReceiversLocal++] = k;

  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  // NRecvSamples = one sample per output step
  dfloat outputInterval;
  settings.getSetting("OUTPUT INTERVAL", outputInterval);

  const dfloat duration = finalTime - startTime;
  NRecvSamples = (dlong)(duration / outputInterval) + 2; // +2 for t=0 and rounding

  sampleRateOut = (int)round(1.0 / outputInterval);

  if (NReceiversLocal > 0) {
    qRecv.malloc(NReceiversLocal * NRecvSamples, 0.0);
    o_qRecv = platform.malloc<dfloat>(qRecv);

    o_recvElements    = platform.malloc<dlong>(recvElements);
    o_recvElementsIdx = platform.malloc<dlong>(recvElementsIdx);

    // compact the located rows into the [NReceiversLocal x Np] layout the
    // kernel indexes by local receiver
    memory<dfloat> intpol(NReceiversLocal * mesh.Np, 0.0);
    for (dlong r = 0; r < NReceiversLocal; ++r) {
      const dlong k = recvElementsIdx[r];
      for (int n = 0; n < mesh.Np; ++n)
        intpol[r*mesh.Np + n] = recvIP[k*mesh.Np + n];
    }
    o_recvIP = platform.malloc<dfloat>(intpol);
  }

  // buildKernel is collective, so every rank builds even with no local receivers
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
