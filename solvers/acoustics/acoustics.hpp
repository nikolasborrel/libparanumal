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
  dfloat LRMaxPole = 0.0;        // fastest pole rate [1/s] in the vectorfit

  memory<dfloat>       LR;       // vectorfit coefficients
  memory<dlong>        LRInfo;   // [Npoles, NRealPoles, NImagPoles]
  memory<dlong>        mapAcc;   // [Nelements*Nfp*Nfaces] face-node → acc row (-1 if not LR)
  memory<dlong>        mapAccToQ;// [NLRPoints] acc row → that node's pressure entry in q
  memory<dfloat>       acc;      // [NLRPoints * LRNpoles] accumulator state (host copy)
  memory<dfloat>       resacc;   // [NLRPoints * LRNpoles] LSERK4 residual

  deviceMemory<dfloat> o_LR;
  deviceMemory<dlong>  o_LRInfo;
  deviceMemory<dlong>  o_mapAcc;
  deviceMemory<dlong>  o_mapAccToQ;
  deviceMemory<dfloat> o_acc;
  deviceMemory<dfloat> o_resacc;
  deviceMemory<dfloat> o_rhsacc; // filled by surfaceKernel each rhsf() call

  kernel_t updateKernelLR;

  // EIRK4 — ARK4(3)6L[2]SA additive Runge-Kutta. Explicit for the wave field,
  // L-stable ESDIRK for the LR accumulators, so the pole rates of the vectorfit
  // no longer constrain dt. See src/acousticsEIRK4.cpp.
  bool useEIRK4 = false;

  deviceMemory<dfloat> o_erkA, o_erkB, o_esdirkA;
  deviceMemory<dfloat> o_eirkResq;  // field stage value, halo-sized
  deviceMemory<dfloat> o_eirkXacc;  // accumulator stage value
  deviceMemory<dfloat> o_eirkKq;    // [6 x N]    field stage derivatives
  deviceMemory<dfloat> o_eirkKacc;  // [6 x Nacc] accumulator stage derivatives

  kernel_t updateKernelEIRK4;
  kernel_t updateKernelEIRK4LR;

  void SetupEIRK4();
  void StepEIRK4(dfloat time, dfloat stepdt);

  // output files are written under outDir, stemmed on simulationID
  std::string simulationID;
  std::string outDir;
  std::string logPath;
  bool        logOpened = false;  // truncate on first write, append thereafter
  int         outputFrame = 0;
  dfloat      dt = 0.0;

  OutputFormat outputFormat = OutputFormat::NONE;
  std::unique_ptr<acousticWriter_t> h5Writer;
  std::vector<dfloat> timeStepsOut;  // planned output times, for preallocation
  std::vector<dfloat> outputTimes;   // times actually written

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

  // resolve output paths; must run before any Logf()
  void SetupOutput();

  // rank 0: print to stdout and append to logPath
  void Logf(const char* fmt, ...) __attribute__((format(printf, 2, 3)));

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(memory<dfloat> Q, const std::string fileName);

  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);

  // rhsf on explicit accumulator buffers, for schemes that stage the LR
  // accumulators separately from o_acc. The two-buffer form above forwards here.
  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs,
            deviceMemory<dfloat>& o_ACC, deviceMemory<dfloat>& o_RHSACC,
            const dfloat time);

  dfloat MaxWaveSpeed();

  // set dt from the CFL bound; called in Setup so cadences can be clamped to it
  void ComputeTimeStep();

  void WriteReceiverIRs();
};

#endif
