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
#include <cmath>
#include <algorithm>
#include <string>
#include <cstdio>

extern "C" {
#include "tinywav.h"
}

#define WAV_BLOCK 512

void acoustics_t::WriteWavIRs() {
  if (NReceiversLocal == 0 || recvSampleIdx == 0) return;

  // Copy device buffer to host
  o_qRecv.copyTo(qRecv);

  std::string outDir = ".";
  if (settings.hasSetting("OUTPUT DIRECTORY"))
    settings.getSetting("OUTPUT DIRECTORY", outDir);

  std::string baseName = "acoustics";
  if (settings.hasSetting("OUTPUT FILE NAME"))
    settings.getSetting("OUTPUT FILE NAME", baseName);

  const dlong nSamples = recvSampleIdx; // actual samples recorded

  for (dlong iRecv = 0; iRecv < NReceiversLocal; iRecv++) {
    // Find peak for normalization
    float maxAmpl = 1.0f;
    for (dlong s = 0; s < nSamples; s++) {
      const float v = (float)qRecv[iRecv * NRecvSamples + s];
      if (std::fabs(v) > maxAmpl) maxAmpl = std::fabs(v);
    }

    char path[BUFSIZ];
    const dlong gIdx = recvElementsIdx[iRecv]; // global receiver index
    snprintf(path, sizeof(path), "%s/%s_recv%04lld.wav",
             outDir.c_str(), baseName.c_str(), (long long)gIdx);

    TinyWav tw;
    tinywav_open_write(&tw,
        1,
        sampleRateOut,
        TW_INT16,
        TW_INLINE,
        path);

    dlong processed = 0;
    while (processed < nSamples) {
      const dlong chunk = std::min((dlong)WAV_BLOCK, nSamples - processed);
      float buf[WAV_BLOCK];
      for (dlong i = 0; i < chunk; i++)
        buf[i] = (float)(qRecv[iRecv * NRecvSamples + processed + i]) / maxAmpl;
      tinywav_write_f(&tw, buf, (int)chunk);
      processed += chunk;
    }

    tinywav_close_write(&tw);

    if (mesh.rank == 0)
      printf("  wrote %s (%lld samples @ %d Hz)\n",
             path, (long long)nSamples, sampleRateOut);
  }
}
