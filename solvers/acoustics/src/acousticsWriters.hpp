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
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

class acoustics_t;

// Deduplicates DG solution nodes across all elements.
// Fills x1d/y1d/z1d (node coords), p1d (pressure at current time), and
// conn[e][n] (element-node → unique index map).  Reads acoustics.q (host copy).
void extractUniquePoints(acoustics_t& ac,
                         std::vector<std::vector<unsigned int>>& conn,
                         std::vector<float>& x1d,
                         std::vector<float>& y1d,
                         std::vector<float>& z1d,
                         std::vector<float>& p1d);

class IAcousticWriter {
public:
    virtual ~IAcousticWriter() = default;
    virtual void write(acoustics_t& ac, int iter) = 0;
};

// Writes all timesteps as a single chunked 2-D HDF5 dataset:
//   /mesh       — [N×3] float32 node coordinates
//   /pressures  — [nFrames×N] float32 pressure snapshots
//   /pressures  attribute "time_steps" — vector of output times
class AcousticH5CompactWriter : public IAcousticWriter {
public:
    explicit AcousticH5CompactWriter(acoustics_t& ac);
    void write(acoustics_t& ac, int iter) override;

private:
    std::string _filepathH5;
    HighFive::DataSet _pressureDataset;

    void writeMesh(const std::string& filepathH5,
                   std::vector<float>& x1d,
                   std::vector<float>& y1d,
                   std::vector<float>& z1d);
};

// Writes one HDF5 dataset per timestep plus an XDMF temporal-collection sidecar
// compatible with ParaView / VisIt:
//   /data0       — [N×3] float32 node coordinates
//   /data1       — [N×1] int topology (point cloud)
//   /data2..     — [N]   float32 pressure per frame
class AcousticXdmfWriter : public IAcousticWriter {
public:
    explicit AcousticXdmfWriter(acoustics_t& ac);
    void write(acoustics_t& ac, int iter) override;

private:
    std::string _filepathH5;
    size_t      _Nnodes = 0;

    void writeXdmfHeader(acoustics_t& ac, size_t Nnodes,
                         const std::string& filepathXdmf,
                         const std::string& filenameH5);

    void writeMesh(const std::string& filepathH5,
                   std::vector<float>& x1d,
                   std::vector<float>& y1d,
                   std::vector<float>& z1d,
                   unsigned int fileAttr);
};
