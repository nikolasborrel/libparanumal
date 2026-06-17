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

#ifndef ACOUSTICS_HPP
#define ACOUSTICS_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"

#include <memory>
#include <string>
#include <vector>

#define DACOUSTICS LIBP_DIR"/solvers/acoustics/"

using namespace libp;

#include "src/acousticsWriters.hpp"

// Wave-field HDF5 snapshot format. VTU is NOT a member: VTU visualization is an
// independent output path controlled by OUTPUT TO FILE (shared with other
// solvers), and may be written alongside an HDF5 format.
enum class OutputFormat { NONE, H5COMPACT, XDMF };

class acousticsSettings_t: public settings_t {
public:
  acousticsSettings_t(comm_t _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename);
};

class acoustics_t: public solver_t {
public:
  mesh_t mesh;

  int Nfields;

  dfloat rho, c;        // medium density [kg/m^3] and speed of sound [m/s]
  dfloat fmax;          // max frequency of interest [Hz] (source/mesh property)
  dfloat sigma0;        // Gaussian source width [m] = SXYZ or 2c/(pi*fmax)
  dfloat srcX, srcY, srcZ;  // Gaussian source center [m] (SOURCE POSITION)
  dfloat ZFreqIndep;    // acoustic impedance for frequency-independent BC [Pa·s/m]

  timeStepper_t timeStepper;

  ogs::halo_t traceHalo;

  memory<dfloat> q;
  deviceMemory<dfloat> o_q;

  deviceMemory<dfloat> o_Mq;

  kernel_t volumeKernel;
  kernel_t surfaceKernel;

  kernel_t initialConditionKernel;

  // Receiver interpolation
  dlong NReceivers = 0;
  dlong NReceiversLocal = 0;
  dlong NRecvSamples = 0;
  dlong recvSampleIdx = 0;
  int   sampleRateOut = 0;

  memory<dfloat> recvXYZ;
  memory<dlong>  recvElements;
  memory<dlong>  recvElementsIdx;
  memory<dfloat> qRecv;            // host: [NReceiversLocal × NRecvSamples]

  deviceMemory<dfloat> o_recvIP;          // interpolation operators [NReceiversLocal × Np]
  deviceMemory<dlong>  o_recvElements;
  deviceMemory<dlong>  o_recvElementsIdx;
  deviceMemory<dfloat> o_qRecv;           // device: [NReceiversLocal × NRecvSamples]

  kernel_t receiverKernel;

  // Locally Reacting (LR) BC — accumulator fields for frequency-dependent impedance
  dlong LRNpoles = 0;
  dlong LRNRealPoles = 0;
  dlong LRNImagPoles = 0;
  dlong NLRPoints = 0;           // total LR boundary face nodes on this rank

  memory<dfloat>       LR;       // vectorfit coefficients
  memory<dlong>        LRInfo;   // [Npoles, NRealPoles, NImagPoles]
  memory<dlong>        mapAcc;   // [Nelements*Nfp*Nfaces] face-node → acc row (-1 if not LR)
  memory<dfloat>       acc;      // [NLRPoints * LRNpoles] accumulator state (host copy)
  memory<dfloat>       resacc;   // [NLRPoints * LRNpoles] LSERK4 residual

  deviceMemory<dfloat> o_LR;
  deviceMemory<dlong>  o_LRInfo;
  deviceMemory<dlong>  o_mapAcc;
  deviceMemory<dfloat> o_acc;
  deviceMemory<dfloat> o_resacc;
  deviceMemory<dfloat> o_rhsacc; // filled by surfaceKernel each rhsf() call

  kernel_t updateKernelLR;

  // ---- ML data-generation sampling (orthogonal to visualization) --------
  // Samples the pressure field onto a uniform points-per-wavelength grid whose
  // coordinates are randomly jittered (seeded), so a downstream ML model sees
  // varied locations and does not overfit to a fixed grid. Driven on its own
  // temporal cadence in the time loop — never coupled to Report()/OUTPUT INTERVAL.
  bool   dataGenEnabled = false;
  dlong  NDataPoints   = 0;       // jittered sample points located on this rank
  dlong  NDataSamples  = 0;       // time frames allocated
  dlong  dataSampleIdx = 0;       // frames recorded so far
  int    dataStride    = 1;       // sample every dataStride solver steps
  dfloat dataDx = 0.0, dataDt = 0.0, dataFmax = 0.0, dataPPW = 0.0, dataJitter = 0.0;
  int    dataSeed = 0;
  std::string dataSimID;

  memory<dfloat> dataXYZ;            // [NDataPoints*3] jittered sample coordinates
  memory<dlong>  dataElements;       // [NDataPoints] containing element id (local)
  memory<dlong>  dataElementsIdx;    // [NDataPoints] identity index (for shared kernel)
  memory<dfloat> dataVals;           // host [NDataPoints*NDataSamples]
  std::vector<dfloat> dataTimes;     // actual sample times recorded

  deviceMemory<dfloat> o_dataIP;           // interpolation operators [NDataPoints*Np]
  deviceMemory<dlong>  o_dataElements;
  deviceMemory<dlong>  o_dataElementsIdx;
  deviceMemory<dfloat> o_dataVals;         // device [NDataPoints*NDataSamples]

  kernel_t dataKernel;

  // Rectilinear SOURCE grid (DeepONet branch input): the initial condition
  // sampled on a uniform PPW grid (full nx*ny*nz, 0 outside the domain),
  // optionally jittered by DATA SOURCE JITTER. Written as /umesh, /upressures,
  // umesh_shape, /source_position (DTU-compatible layout).
  dfloat srcPPW = 0.0, srcJitter = 0.0;
  dlong  NSrcPoints = 0;
  int    srcShape[3] = {0, 0, 0};          // nx, ny, nz of the rectilinear grid
  memory<dfloat> srcGridXYZ;               // [NSrcPoints*3] grid coordinates
  memory<dlong>  srcGridElem;              // [NSrcPoints] containing element (-1 outside)
  memory<dfloat> srcGridIP;                // [NSrcPoints*Np] Lagrange weights
  memory<float>  srcGridVals;              // [NSrcPoints] sampled IC pressure (host)

  // HDF5/XDMF output
  OutputFormat outputFormat = OutputFormat::NONE;
  std::string  simulationID;
  std::string  outDir;
  std::string  logPath;            // <outDir>/<simID>.log mirror of key diagnostics
  bool         logOpened = false;  // truncate on first write, append thereafter
  std::vector<dfloat> timeStepsOut;  // precomputed (clamped) output times for preallocation
  std::vector<dfloat> outputTimes;   // actual wave-field output times (XDMF header finalize)
  int    outputFrame = 0;            // wave-field frames written so far
  dfloat dt = 0.0;                   // explicit solver time step (set in Setup)
  std::unique_ptr<IAcousticWriter> h5Writer;

  acoustics_t() = default;
  acoustics_t(platform_t &_platform, mesh_t &_mesh,
              acousticsSettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh,
             acousticsSettings_t& _settings);

  void SetupReceivers();
  void SetupLRBC(properties_t& kernelInfo);
  void SetupHDF5Output();

  // printf-style diagnostic that (on rank 0) prints to stdout AND appends to the
  // run log file (logPath). Truncates the log on first use. No-op off rank 0.
  void Logf(const char* fmt, ...) __attribute__((format(printf, 2, 3)));

  // ML data-generation: build the jittered PPW sample grid + interpolation
  // operators (no-op unless DATA OUTPUT is TRUE).
  void SetupDataGen();
  // Interpolate the pressure field onto the sample grid into frame `frameIdx`.
  void SampleDataGrid(int frameIdx, dfloat time);
  // Sample the initial condition onto the rectilinear source grid (call at t=0).
  void SampleSourceGrid();
  // Persist the sampled grid to <DATA SIMULATION ID>_data.h5.
  void WriteDataGrid();

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(memory<dfloat> Q, const std::string fileName);

  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);

  dfloat MaxWaveSpeed();

  // Compute the explicit CFL time step from the mesh / wave speed and store it
  // in `dt`. Called in Setup() so output cadences can be clamped to it.
  void ComputeTimeStep();

  // Write sampled receiver impulse responses to <SIMULATION ID>_receivers.h5
  // (no-op if no receivers are configured). Includes sample rate + positions.
  void WriteReceiverIRs();
};

#endif
