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
#include "acousticsWriters.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <cstdio>
#include <unordered_map>
#include <vector>
#include <highfive/H5Easy.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5PropertyList.hpp>

using namespace HighFive;

// ---------------------------------------------------------------------------
// WriteReceiverIRs — persist sampled receiver impulse responses to HDF5
//
// Writes <OUTPUT DIRECTORY>/<SIMULATION ID>_receivers.h5 containing:
//   /impulse_responses  [NReceiversLocal x nSamples]  raw pressure time series
//   /positions          [NReceiversLocal x 3]         receiver XYZ (physical)
//   /receiver_index     [NReceiversLocal]             global receiver index
// with attributes on /impulse_responses: sample_rate_hz, n_samples, n_receivers.
// ---------------------------------------------------------------------------

void acoustics_t::WriteReceiverIRs() {
  if (NReceiversLocal == 0 || recvSampleIdx == 0) return;

  // Pull the sampled receiver buffer back to the host
  o_qRecv.copyTo(qRecv);

  const std::string dir = outDir.empty() ? std::string(".") : outDir;
  std::string simID = "acoustics";
  if (settings.hasSetting("SIMULATION ID"))
    settings.getSetting("SIMULATION ID", simID);
  const std::string path = dir + "/" + simID + "_receivers.h5";

  const size_t nSamples = static_cast<size_t>(recvSampleIdx); // samples recorded
  const size_t nRecv    = static_cast<size_t>(NReceiversLocal);

  // Raw (un-normalised) impulse responses, plus positions and global indices
  std::vector<std::vector<double>> irs(nRecv, std::vector<double>(nSamples));
  std::vector<std::vector<double>> pos(nRecv, std::vector<double>(3));
  std::vector<int> gIdx(nRecv);
  for (size_t r = 0; r < nRecv; ++r) {
    for (size_t s = 0; s < nSamples; ++s)
      irs[r][s] = static_cast<double>(qRecv[r * NRecvSamples + s]);
    const dlong g = recvElementsIdx[r];          // global receiver index
    gIdx[r]   = static_cast<int>(g);
    pos[r][0] = static_cast<double>(recvXYZ[g * 3 + 0]);
    pos[r][1] = static_cast<double>(recvXYZ[g * 3 + 1]);
    pos[r][2] = static_cast<double>(recvXYZ[g * 3 + 2]);
  }

  File file(path, File::Overwrite);
  file.createDataSet<double>("/impulse_responses", DataSpace::From(irs)).write(irs);
  file.createDataSet<double>("/positions",         DataSpace::From(pos)).write(pos);
  file.createDataSet<int>   ("/receiver_index",    DataSpace::From(gIdx)).write(gIdx);

  H5Easy::dumpAttribute(file, "/impulse_responses", "sample_rate_hz", sampleRateOut);
  H5Easy::dumpAttribute(file, "/impulse_responses", "n_samples",   static_cast<int>(nSamples));
  H5Easy::dumpAttribute(file, "/impulse_responses", "n_receivers", static_cast<int>(nRecv));

  if (mesh.rank == 0)
    printf("  wrote %s (%zu receivers x %zu samples @ %d Hz)\n",
           path.c_str(), nRecv, nSamples, sampleRateOut);
}

// ---------------------------------------------------------------------------
// WriteDataGrid — persist the jittered PPW sample grid for ML dataset use to
// <OUTPUT DIRECTORY>/<DATA SIMULATION ID>_data.h5 containing:
//   /coordinates  [NDataPoints x 3]        jittered sample-point XYZ (physical)
//   /pressures    [nFrames x NDataPoints]  pressure at each sample point per frame
//   /times        [nFrames]                actual sample times
// with attributes on /pressures: dx, fmax, spatial_ppw, jitter_frac, seed.
// Orthogonal to the visualization writers — never coupled to OUTPUT INTERVAL.
// ---------------------------------------------------------------------------
void acoustics_t::WriteDataGrid() {
  if (!dataGenEnabled) return;
  const bool haveQuery  = (NDataPoints > 0 && dataSampleIdx > 0);
  const bool haveSource = (NSrcPoints  > 0);
  if (!haveQuery && !haveSource) return;

  const std::string dir = outDir.empty() ? std::string(".") : outDir;
  const std::string path = dir + "/" + dataSimID + "_data.h5";

  File file(path, File::Overwrite);

  // --- query / output grid (DeepONet trunk): /mesh + /pressures over time ---
  // Datasets mirror the DTU layout (/mesh, /pressures with a time_steps attr).
  if (haveQuery) {
    o_dataVals.copyTo(dataVals);    // host layout [point * NDataSamples + frame]
    const size_t nFrames = static_cast<size_t>(dataSampleIdx);
    const size_t nPts    = static_cast<size_t>(NDataPoints);

    std::vector<std::vector<float>> mesh1d(nPts, std::vector<float>(3));
    for (size_t p = 0; p < nPts; ++p) {
      mesh1d[p][0] = static_cast<float>(dataXYZ[p*3+0]);
      mesh1d[p][1] = static_cast<float>(dataXYZ[p*3+1]);
      mesh1d[p][2] = static_cast<float>(dataXYZ[p*3+2]);
    }
    // frame-major [nFrames x nPts]
    std::vector<std::vector<float>> pres(nFrames, std::vector<float>(nPts));
    for (size_t f = 0; f < nFrames; ++f)
      for (size_t p = 0; p < nPts; ++p)
        pres[f][p] = static_cast<float>(dataVals[p * NDataSamples + f]);

    file.createDataSet<float>("/mesh",      DataSpace::From(mesh1d)).write(mesh1d);
    file.createDataSet<float>("/pressures", DataSpace::From(pres)).write(pres);

    H5Easy::dumpAttribute(file, "/pressures", "time_steps",  dataTimes);
    H5Easy::dumpAttribute(file, "/pressures", "dx",          static_cast<double>(dataDx));
    H5Easy::dumpAttribute(file, "/pressures", "fmax",        static_cast<double>(dataFmax));
    H5Easy::dumpAttribute(file, "/pressures", "spatial_ppw", static_cast<double>(dataPPW));
    H5Easy::dumpAttribute(file, "/pressures", "jitter_frac", static_cast<double>(dataJitter));
    H5Easy::dumpAttribute(file, "/pressures", "seed",        dataSeed);
    H5Easy::dumpAttribute(file, "/pressures", "n_frames",    static_cast<int>(nFrames));
    H5Easy::dumpAttribute(file, "/pressures", "n_points",    static_cast<int>(nPts));
  }

  // --- source grid (DeepONet branch): /umesh + /upressures (IC), DTU layout ---
  if (haveSource) {
    const size_t nu = static_cast<size_t>(NSrcPoints);
    std::vector<std::vector<float>> umesh(nu, std::vector<float>(3));
    std::vector<float> upres(nu);
    for (size_t p = 0; p < nu; ++p) {
      umesh[p][0] = static_cast<float>(srcGridXYZ[p*3+0]);
      umesh[p][1] = static_cast<float>(srcGridXYZ[p*3+1]);
      umesh[p][2] = static_cast<float>(srcGridXYZ[p*3+2]);
      upres[p]    = srcGridVals[p];
    }
    file.createDataSet<float>("/umesh",      DataSpace::From(umesh)).write(umesh);
    file.createDataSet<float>("/upressures", DataSpace::From(upres)).write(upres);

    std::vector<int> umesh_shape = {srcShape[0], srcShape[1], srcShape[2]};
    H5Easy::dumpAttribute(file, "/umesh", "umesh_shape", umesh_shape);
    H5Easy::dumpAttribute(file, "/umesh", "source_ppw",  static_cast<double>(srcPPW));
    H5Easy::dumpAttribute(file, "/umesh", "jitter_frac", static_cast<double>(srcJitter));

    std::vector<float> src_pos = {static_cast<float>(srcX),
                                  static_cast<float>(srcY),
                                  static_cast<float>(srcZ)};
    file.createDataSet<float>("/source_position", DataSpace::From(src_pos)).write(src_pos);
  }

  if (mesh.rank == 0)
    printf("  wrote %s (query %lld pts x %lld frames; source %lld pts %dx%dx%d)\n",
           path.c_str(),
           (long long)(haveQuery ? NDataPoints : 0),
           (long long)(haveQuery ? dataSampleIdx : 0),
           (long long)(haveSource ? NSrcPoints : 0),
           srcShape[0], srcShape[1], srcShape[2]);
}

// ---------------------------------------------------------------------------
// extractUniquePoints
// ---------------------------------------------------------------------------

void extractUniquePoints(acoustics_t& ac,
                         std::vector<std::vector<unsigned int>>& conn,
                         std::vector<float>& x1d,
                         std::vector<float>& y1d,
                         std::vector<float>& z1d,
                         std::vector<float>& p1d)
{
  // dfloat is defined as a macro (double) by libparanumal/include/types.h
  struct Coord3D { dfloat x, y, z; };

  struct Coord3DHasher {
    size_t operator()(const Coord3D& c) const {
      // quantise to 1e-5 resolution before hashing
      auto h = [](dfloat v) -> size_t {
        return std::hash<long long>()(static_cast<long long>(std::llround(v * 1e5)));
      };
      return h(c.x) ^ (h(c.y) << 1) ^ (h(c.z) << 2);
    }
  };

  auto eq = [](const Coord3D& a, const Coord3D& b) -> bool {
    const dfloat eps = 1e-5;
    return std::fabs(a.x - b.x) < eps &&
           std::fabs(a.y - b.y) < eps &&
           std::fabs(a.z - b.z) < eps;
  };

  std::unordered_map<Coord3D, unsigned int, Coord3DHasher, decltype(eq)> coordMap(0, Coord3DHasher{}, eq);

  conn.assign(ac.mesh.Nelements, std::vector<unsigned int>(ac.mesh.Np));

  for (dlong e = 0; e < ac.mesh.Nelements; ++e) {
    for (int n = 0; n < ac.mesh.Np; ++n) {
      dfloat x = ac.mesh.x[n + ac.mesh.Np * e];
      dfloat y = ac.mesh.y[n + ac.mesh.Np * e];
      dfloat z = ac.mesh.z[n + ac.mesh.Np * e];

      auto [it, inserted] = coordMap.emplace(Coord3D{x, y, z},
                                             static_cast<unsigned int>(x1d.size()));
      if (inserted) {
        x1d.push_back(static_cast<float>(x));
        y1d.push_back(static_cast<float>(y));
        z1d.push_back(static_cast<float>(z));
        // pressure is field 0; layout: q[e*Np*Nfields + field*Np + n]
        dlong qbase = e * ac.mesh.Np * ac.Nfields + n;
        p1d.push_back(static_cast<float>(ac.q[qbase]));
      }
      conn[e][n] = it->second;
    }
  }
}

// ---------------------------------------------------------------------------
// AcousticH5CompactWriter
// ---------------------------------------------------------------------------

void AcousticH5CompactWriter::writeMesh(const std::string& filepathH5,
                                        std::vector<float>& x1d,
                                        std::vector<float>& y1d,
                                        std::vector<float>& z1d)
{
  File file(filepathH5, File::Overwrite);
  size_t N = x1d.size();
  std::vector<std::vector<float>> meshData(N, std::vector<float>(3));
  for (size_t i = 0; i < N; ++i) {
    meshData[i][0] = x1d[i];
    meshData[i][1] = y1d[i];
    meshData[i][2] = z1d[i];
  }
  DataSet ds = file.createDataSet<float>("/mesh", DataSpace::From(meshData));
  ds.write(meshData);
}

AcousticH5CompactWriter::AcousticH5CompactWriter(acoustics_t& ac)
{
  std::string simID = "acoustics";
  if (ac.settings.hasSetting("SIMULATION ID"))
    ac.settings.getSetting("SIMULATION ID", simID);

  _filepathH5 = ac.outDir + "/" + simID + ".h5";

  auto conn  = std::vector<std::vector<unsigned int>>();
  auto x1d   = std::vector<float>();
  auto y1d   = std::vector<float>();
  auto z1d   = std::vector<float>();
  auto p1d   = std::vector<float>();
  extractUniquePoints(ac, conn, x1d, y1d, z1d, p1d);

  writeMesh(_filepathH5, x1d, y1d, z1d);

  size_t nFrames = ac.timeStepsOut.size();
  size_t nPts    = x1d.size();
  _nFrames = nFrames;

  DataSetCreateProps props;
  props.add(Chunking(std::vector<hsize_t>{1, static_cast<hsize_t>(nPts)}));

  File file(_filepathH5, File::OpenOrCreate);
  _pressureDataset = file.createDataSet<float>(
      "/pressures",
      DataSpace(std::vector<size_t>{nFrames, nPts}),
      props);

  H5Easy::dumpAttribute(file, "/pressures", "time_steps", ac.timeStepsOut);
}

void AcousticH5CompactWriter::write(acoustics_t& ac, int iter)
{
  auto conn  = std::vector<std::vector<unsigned int>>();
  auto x1d   = std::vector<float>();
  auto y1d   = std::vector<float>();
  auto z1d   = std::vector<float>();
  auto p1d   = std::vector<float>();

  ac.o_q.copyTo(ac.q);
  extractUniquePoints(ac, conn, x1d, y1d, z1d, p1d);

  // Guard against emitting more frames than were preallocated (boundary
  // floating-point rounding in the output-crossing logic).
  if (static_cast<size_t>(iter) >= _nFrames) return;

  _pressureDataset.select({static_cast<size_t>(iter), 0},
                          {1, p1d.size()}).write(p1d);
}

// ---------------------------------------------------------------------------
// AcousticXdmfWriter
// ---------------------------------------------------------------------------

void AcousticXdmfWriter::writeMesh(const std::string& filepathH5,
                                   std::vector<float>& x1d,
                                   std::vector<float>& y1d,
                                   std::vector<float>& z1d,
                                   unsigned int fileAttr)
{
  File file(filepathH5, fileAttr);
  size_t N = x1d.size();

  std::vector<std::vector<double>> meshData(N, std::vector<double>(3));
  std::vector<int> vertData(N);
  for (size_t i = 0; i < N; ++i) {
    meshData[i][0] = x1d[i];
    meshData[i][1] = y1d[i];
    meshData[i][2] = z1d[i];
    vertData[i] = static_cast<int>(i);
  }

  DataSet ds = file.createDataSet<double>("/data0", DataSpace::From(meshData));
  ds.write(meshData);
  DataSet ds2 = file.createDataSet<int>("/data1", DataSpace::From(vertData));
  ds2.write(vertData);
}

void AcousticXdmfWriter::writeXdmfHeader(acoustics_t& ac,
                                          size_t Nnodes,
                                          const std::string& filepathXdmf,
                                          const std::string& filenameH5)
{
  std::ofstream ofs(filepathXdmf);

  ofs << "<?xml version=\"1.0\" ?>\n";
  ofs << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
  ofs << "<Xdmf Version=\"3.0\">\n";
  ofs << "  <Domain>\n";
  ofs << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

  // Build one temporal grid per frame that was actually written to the HDF5 file
  // (ac.outputTimes), so the series never references a /dataN that is missing.
  for (size_t i = 0; i < ac.outputTimes.size(); ++i) {
    std::string tag = "/data" + std::to_string(i + 2);
    ofs << "      <Grid>\n";
    ofs << "        <include xpointer=\"xpointer(//Grid[@Name=&quot;mesh&quot;]/*[self::Topology or self::Geometry])\" />\n";
    ofs << "        <Time Value=\"" << ac.outputTimes[i] << "\" />\n";
    ofs << "        <Attribute Name=\"p\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    ofs << "          <DataItem DataType=\"Float\" Dimensions=\"" << Nnodes << "\" Format=\"HDF\" Precision=\"8\">\n";
    ofs << "            " << filenameH5 << ":" << tag << "\n";
    ofs << "          </DataItem>\n";
    ofs << "        </Attribute>\n";
    ofs << "      </Grid>\n";
  }

  ofs << "    </Grid>\n";
  ofs << "    <Grid Name=\"mesh\" GridType=\"Uniform\">\n";
  ofs << "      <Geometry GeometryType=\"XYZ\">\n";
  ofs << "        <DataItem DataType=\"Float\" Dimensions=\"" << Nnodes << " 3\" Format=\"HDF\" Precision=\"8\">\n";
  ofs << "          " << filenameH5 << ":/data0\n";
  ofs << "        </DataItem>\n";
  ofs << "      </Geometry>\n";
  ofs << "      <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Nnodes << "\">\n";
  ofs << "        <DataItem DataType=\"Int\" Dimensions=\"" << Nnodes << " 1\" Format=\"HDF\" Precision=\"8\">\n";
  ofs << "          " << filenameH5 << ":/data1\n";
  ofs << "        </DataItem>\n";
  ofs << "      </Topology>\n";
  ofs << "    </Grid>\n";
  ofs << "  </Domain>\n";
  ofs << "</Xdmf>\n";
}

AcousticXdmfWriter::AcousticXdmfWriter(acoustics_t& ac)
{
  std::string simID = "acoustics";
  if (ac.settings.hasSetting("SIMULATION ID"))
    ac.settings.getSetting("SIMULATION ID", simID);

  _filenameH5  = simID + ".h5";
  _filepathH5  = ac.outDir + "/" + _filenameH5;
  _filepathXdmf = ac.outDir + "/" + simID + ".xdmf";

  auto conn  = std::vector<std::vector<unsigned int>>();
  auto x1d   = std::vector<float>();
  auto y1d   = std::vector<float>();
  auto z1d   = std::vector<float>();
  auto p1d   = std::vector<float>();
  extractUniquePoints(ac, conn, x1d, y1d, z1d, p1d);

  _Nnodes = x1d.size();

  writeMesh(_filepathH5, x1d, y1d, z1d, File::Overwrite);

  // The XDMF sidecar is written in finalize(), once the number of frames
  // actually produced is known (ac.outputTimes). Writing it here would bake in
  // the planned frame count, which can exceed what the run emits.
}

void AcousticXdmfWriter::write(acoustics_t& ac, int iter)
{
  auto conn  = std::vector<std::vector<unsigned int>>();
  auto x1d   = std::vector<float>();
  auto y1d   = std::vector<float>();
  auto z1d   = std::vector<float>();
  auto p1d   = std::vector<float>();

  ac.o_q.copyTo(ac.q);
  extractUniquePoints(ac, conn, x1d, y1d, z1d, p1d);

  std::string tag = "/data" + std::to_string(iter + 2);
  H5Easy::File file(_filepathH5, H5Easy::File::OpenOrCreate);
  H5Easy::dump(file, tag, p1d);
}

void AcousticXdmfWriter::finalize(acoustics_t& ac)
{
  // Emit the temporal-collection header for exactly the frames written.
  writeXdmfHeader(ac, _Nnodes, _filepathXdmf, _filenameH5);
}
