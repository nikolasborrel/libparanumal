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

#pragma once

#include <memory>
#include <string>
#include <vector>

#ifdef LIBP_HDF5
#include <highfive/H5DataSet.hpp>
#endif

class acoustics_t;

// Wave-field snapshot format. VTU is not a member: it is an independent output
// path (OUTPUT TO FILE, shared with the other solvers) and may run alongside.
enum class OutputFormat { NONE, H5COMPACT, XDMF };

class acousticWriter_t {
public:
  virtual ~acousticWriter_t() = default;
  virtual void write(acoustics_t& ac, int frame) = 0;
  // flush metadata that depends on the number of frames actually written
  virtual void finalize(acoustics_t& ac) {}
};

#ifdef LIBP_HDF5

// DG nodes deduplicated across elements. qidx[i] indexes acoustics.q for the
// pressure at unique point i.
struct uniquePoints_t {
  std::vector<float> x, y, z;
  std::vector<dlong> qidx;
};

uniquePoints_t buildUniquePoints(acoustics_t& ac);
void gatherPressures(acoustics_t& ac, const uniquePoints_t& pts,
                     std::vector<float>& p);

// Single chunked 2-D dataset for the whole run:
//   /mesh      [N x 3]        float32 node coordinates
//   /pressures [nFrames x N]  float32 snapshots, attribute "time_steps"
class h5CompactWriter_t : public acousticWriter_t {
public:
  explicit h5CompactWriter_t(acoustics_t& ac);
  void write(acoustics_t& ac, int frame) override;

private:
  uniquePoints_t     _pts;
  HighFive::DataSet  _pressures;
  size_t             _nFrames = 0;
};

// One dataset per snapshot plus an XDMF temporal-collection sidecar for
// ParaView / VisIt:
//   /data0    [N x 3]  node coordinates
//   /data1    [N]      point-cloud topology
//   /data2..  [N]      float32 pressure per frame
class xdmfWriter_t : public acousticWriter_t {
public:
  explicit xdmfWriter_t(acoustics_t& ac);
  void write(acoustics_t& ac, int frame) override;
  void finalize(acoustics_t& ac) override;

private:
  uniquePoints_t _pts;
  std::string    _pathH5;
  std::string    _pathXdmf;
  std::string    _nameH5;
};

#endif // LIBP_HDF5
