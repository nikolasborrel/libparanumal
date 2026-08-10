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

#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>
#include <highfive/H5Easy.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5PropertyList.hpp>

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

// Deduplicate DG nodes shared between elements, quantising coordinates to 1e-5
// before hashing.
uniquePoints_t buildUniquePoints(acoustics_t& ac) {
  mesh_t& mesh = ac.mesh;
  const bool is3D = (mesh.dim == 3);

  struct key_t { long long x, y, z; };
  struct hash_t {
    size_t operator()(const key_t& k) const {
      std::hash<long long> h;
      return h(k.x) ^ (h(k.y) << 1) ^ (h(k.z) << 2);
    }
  };
  struct eq_t {
    bool operator()(const key_t& a, const key_t& b) const {
      return a.x == b.x && a.y == b.y && a.z == b.z;
    }
  };
  auto quantise = [](dfloat v) { return (long long)std::llround(v * 1e5); };

  uniquePoints_t pts;
  std::unordered_map<key_t, unsigned, hash_t, eq_t> seen;

  for (dlong e = 0; e < mesh.Nelements; ++e) {
    for (int n = 0; n < mesh.Np; ++n) {
      const dlong id = n + mesh.Np * e;
      const dfloat x = mesh.x[id], y = mesh.y[id];
      const dfloat z = is3D ? mesh.z[id] : 0.0;

      auto [it, inserted] = seen.emplace(key_t{quantise(x), quantise(y), quantise(z)},
                                         (unsigned)pts.x.size());
      if (!inserted) continue;

      pts.x.push_back((float)x);
      pts.y.push_back((float)y);
      pts.z.push_back((float)z);
      pts.qidx.push_back(e * mesh.Np * ac.Nfields + n);  // pressure is field 0
    }
  }
  return pts;
}

void gatherPressures(acoustics_t& ac, const uniquePoints_t& pts,
                     std::vector<float>& p) {
  ac.o_q.copyTo(ac.q);
  p.resize(pts.qidx.size());
  for (size_t i = 0; i < pts.qidx.size(); ++i)
    p[i] = (float)ac.q[pts.qidx[i]];
}

h5CompactWriter_t::h5CompactWriter_t(acoustics_t& ac)
{
  const std::string path = ac.outDir + "/" + ac.simulationID + ".h5";
  _pts     = buildUniquePoints(ac);
  _nFrames = ac.timeStepsOut.size();

  const size_t nPts = _pts.x.size();
  std::vector<std::vector<float>> coords(nPts, std::vector<float>(3));
  for (size_t i = 0; i < nPts; ++i)
    coords[i] = {_pts.x[i], _pts.y[i], _pts.z[i]};

  File file(path, File::Overwrite);
  file.createDataSet<float>("/mesh", DataSpace::From(coords)).write(coords);

  DataSetCreateProps props;
  props.add(Chunking(std::vector<hsize_t>{1, (hsize_t)nPts}));
  _pressures = file.createDataSet<float>(
      "/pressures", DataSpace(std::vector<size_t>{_nFrames, nPts}), props);

  H5Easy::dumpAttribute(file, "/pressures", "time_steps", ac.timeStepsOut);
}

void h5CompactWriter_t::write(acoustics_t& ac, int frame)
{
  if ((size_t)frame >= _nFrames) return;   // preallocated rows are the limit

  std::vector<float> p;
  gatherPressures(ac, _pts, p);
  _pressures.select({(size_t)frame, 0}, {1, p.size()}).write(p);
}

xdmfWriter_t::xdmfWriter_t(acoustics_t& ac)
{
  _nameH5   = ac.simulationID + ".h5";
  _pathH5   = ac.outDir + "/" + _nameH5;
  _pathXdmf = ac.outDir + "/" + ac.simulationID + ".xdmf";
  _pts      = buildUniquePoints(ac);

  const size_t nPts = _pts.x.size();
  std::vector<std::vector<double>> coords(nPts, std::vector<double>(3));
  std::vector<int> verts(nPts);
  for (size_t i = 0; i < nPts; ++i) {
    coords[i] = {_pts.x[i], _pts.y[i], _pts.z[i]};
    verts[i]  = (int)i;
  }

  File file(_pathH5, File::Overwrite);
  file.createDataSet<double>("/data0", DataSpace::From(coords)).write(coords);
  file.createDataSet<int>   ("/data1", DataSpace::From(verts)).write(verts);

  // the sidecar is written in finalize(), once the frame count is known
}

void xdmfWriter_t::write(acoustics_t& ac, int frame)
{
  std::vector<float> p;
  gatherPressures(ac, _pts, p);

  H5Easy::File file(_pathH5, H5Easy::File::OpenOrCreate);
  H5Easy::dump(file, "/data" + std::to_string(frame + 2), p);
}

void xdmfWriter_t::finalize(acoustics_t& ac)
{
  const size_t nPts = _pts.x.size();
  std::ofstream ofs(_pathXdmf);

  ofs << "<?xml version=\"1.0\" ?>\n";
  ofs << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
  ofs << "<Xdmf Version=\"3.0\">\n";
  ofs << "  <Domain>\n";
  ofs << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

  // one grid per frame actually written, so the series never references a
  // /dataN that is missing from the file
  for (size_t i = 0; i < ac.outputTimes.size(); ++i) {
    ofs << "      <Grid>\n";
    ofs << "        <include xpointer=\"xpointer(//Grid[@Name=&quot;mesh&quot;]/*[self::Topology or self::Geometry])\" />\n";
    ofs << "        <Time Value=\"" << ac.outputTimes[i] << "\" />\n";
    ofs << "        <Attribute Name=\"p\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    ofs << "          <DataItem DataType=\"Float\" Dimensions=\"" << nPts << "\" Format=\"HDF\" Precision=\"8\">\n";
    ofs << "            " << _nameH5 << ":/data" << (i + 2) << "\n";
    ofs << "          </DataItem>\n";
    ofs << "        </Attribute>\n";
    ofs << "      </Grid>\n";
  }

  ofs << "    </Grid>\n";
  ofs << "    <Grid Name=\"mesh\" GridType=\"Uniform\">\n";
  ofs << "      <Geometry GeometryType=\"XYZ\">\n";
  ofs << "        <DataItem DataType=\"Float\" Dimensions=\"" << nPts << " 3\" Format=\"HDF\" Precision=\"8\">\n";
  ofs << "          " << _nameH5 << ":/data0\n";
  ofs << "        </DataItem>\n";
  ofs << "      </Geometry>\n";
  ofs << "      <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << nPts << "\">\n";
  ofs << "        <DataItem DataType=\"Int\" Dimensions=\"" << nPts << " 1\" Format=\"HDF\" Precision=\"8\">\n";
  ofs << "          " << _nameH5 << ":/data1\n";
  ofs << "        </DataItem>\n";
  ofs << "      </Topology>\n";
  ofs << "    </Grid>\n";
  ofs << "  </Domain>\n";
  ofs << "</Xdmf>\n";
}

// One group per sample set in <SIMULATION ID>_samples.h5:
//   /<name>/points [npts x 3], /<name>/values [nframes x npts], /<name>/times
// with the grid shape and generation parameters as attributes on /<name>/points.
void acoustics_t::WriteSampleSets() {
  if (sampleSets.empty() || mesh.rank != 0) return;

  const std::string path = outDir + "/" + simulationID + "_samples.h5";
  File file(path, File::Overwrite);

  for (sampleSet_t& set : sampleSets) {
    if (set.Npoints == 0 || set.frame == 0) continue;
    set.o_values.copyTo(set.values);

    const size_t npts    = (size_t)set.Npoints;
    const size_t nframes = (size_t)set.frame;

    std::vector<std::vector<float>> pts(npts, std::vector<float>(3));
    for (size_t p = 0; p < npts; ++p)
      pts[p] = {(float)set.xyz[p*3+0], (float)set.xyz[p*3+1], (float)set.xyz[p*3+2]};

    std::vector<std::vector<float>> vals(nframes, std::vector<float>(npts));
    for (size_t f = 0; f < nframes; ++f)
      for (size_t p = 0; p < npts; ++p)
        vals[f][p] = (float)set.values[p*set.Nframes + f];

    const std::string g = "/" + set.name;
    file.createDataSet<float> (g + "/points", DataSpace::From(pts)).write(pts);
    file.createDataSet<float> (g + "/values", DataSpace::From(vals)).write(vals);
    file.createDataSet<double>(g + "/times",  DataSpace::From(set.times)).write(set.times);

    // C order, so values reshape directly: x varies fastest, z slowest
    std::vector<int> shape = {set.shape[2], set.shape[1], set.shape[0]};
    H5Easy::dumpAttribute(file, g + "/points", "grid_shape",  shape);
    H5Easy::dumpAttribute(file, g + "/points", "dx",          (double)set.dx);
    H5Easy::dumpAttribute(file, g + "/points", "jitter",      (double)set.jitter);
    H5Easy::dumpAttribute(file, g + "/points", "seed",        set.seed);
    H5Easy::dumpAttribute(file, g + "/points", "keep_outside", set.keepOutside ? 1 : 0);

    Logf("  wrote %s:%s (%zu points x %zu frames)\n",
         path.c_str(), g.c_str(), npts, nframes);
  }
}

#else  // !LIBP_HDF5

void acoustics_t::WriteReceiverIRs() {
  if (NReceivers == 0) return;
  Logf("  receiver impulse responses not written: built without HDF5 "
       "(rebuild with \"make HDF5=1\")\n");
}

void acoustics_t::WriteSampleSets() {}

#endif // LIBP_HDF5
