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

// HDF5 output paths, compiled only when built with "make HDF5=1" (LIBP_HDF5).

#include "acoustics.hpp"

#ifdef LIBP_HDF5

#include <limits>
#include <string>
#include <vector>
#include <highfive/H5Easy.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

using namespace HighFive;

// Collective. Gathers the distributed receiver record, then rank 0 writes
// <OUTPUT DIRECTORY>/<SIMULATION ID>_receivers.h5 with /impulse_responses
// [NReceivers x nSamples] and /positions [NReceivers x 3]. Receivers keep their
// RECEIVER FILE order, so the layout is rank-count independent.
void acoustics_t::WriteReceiverIRs() {
  if (NReceivers == 0) return;   // every rank must reach the reductions below

  // ranks holding no receiver never advance recvSampleIdx
  dlong nSamplesG = recvSampleIdx;
  mesh.comm.Allreduce(nSamplesG, Comm::Max);
  if (nSamplesG == 0) return;

  // a receiver on a shared face is located by every adjacent rank; keep the
  // lowest so the sum below counts it once
  memory<int> owner(NReceivers, std::numeric_limits<int>::max());
  for (dlong r = 0; r < NReceiversLocal; ++r)
    owner[recvElementsIdx[r]] = mesh.rank;
  mesh.comm.Allreduce(owner, Comm::Min);

  if (NReceiversLocal > 0) o_qRecv.copyTo(qRecv);
  memory<dfloat> all(NReceivers * nSamplesG, 0.0);
  for (dlong r = 0; r < NReceiversLocal; ++r) {
    const dlong g = recvElementsIdx[r];
    if (owner[g] != mesh.rank) continue;
    for (dlong s = 0; s < recvSampleIdx; ++s)
      all[g * nSamplesG + s] = qRecv[r * NRecvSamples + s];
  }
  mesh.comm.Allreduce(all, Comm::Sum);

  if (mesh.rank != 0) return;

  const std::string path = outDir + "/" + simulationID + "_receivers.h5";
  const size_t nSamples = static_cast<size_t>(nSamplesG);
  const size_t nRecv    = static_cast<size_t>(NReceivers);

  std::vector<std::vector<double>> irs(nRecv, std::vector<double>(nSamples));
  std::vector<std::vector<double>> pos(nRecv, std::vector<double>(3));
  for (size_t r = 0; r < nRecv; ++r) {
    for (size_t s = 0; s < nSamples; ++s)
      irs[r][s] = static_cast<double>(all[r * nSamples + s]);
    pos[r][0] = static_cast<double>(recvXYZ[r * 3 + 0]);
    pos[r][1] = static_cast<double>(recvXYZ[r * 3 + 1]);
    pos[r][2] = static_cast<double>(recvXYZ[r * 3 + 2]);
  }

  File file(path, File::Overwrite);
  file.createDataSet<double>("/impulse_responses", DataSpace::From(irs)).write(irs);
  file.createDataSet<double>("/positions",         DataSpace::From(pos)).write(pos);

  H5Easy::dumpAttribute(file, "/impulse_responses", "sample_rate_hz", sampleRateOut);
  H5Easy::dumpAttribute(file, "/impulse_responses", "n_samples",   static_cast<int>(nSamples));
  H5Easy::dumpAttribute(file, "/impulse_responses", "n_receivers", static_cast<int>(nRecv));

  Logf("  wrote %s (%zu receivers x %zu samples @ %d Hz)\n",
       path.c_str(), nRecv, nSamples, sampleRateOut);
}

#else  // !LIBP_HDF5

void acoustics_t::WriteReceiverIRs() {
  if (NReceivers == 0) return;
  Logf("  receiver impulse responses not written: built without HDF5 "
       "(rebuild with \"make HDF5=1\")\n");
}

#endif // LIBP_HDF5
