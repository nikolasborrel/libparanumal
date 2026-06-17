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
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>

void acoustics_t::Setup(platform_t& _platform, mesh_t& _mesh,
                        acousticsSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  Nfields = (mesh.dim==3) ? 4:3;

  settings.getSetting("DENSITY", rho);
  settings.getSetting("SPEED OF SOUND", c);
  settings.getSetting("FMAX", fmax);
  settings.getSetting("FREQINDEP IMPEDANCE", ZFreqIndep);

  // Gaussian source width: explicit SXYZ if given (>0), else from FMAX using the
  // DTU convention sigma = c / (pi*fmax/2) = 2c/(pi*fmax). Exposed to FMAX-driven
  // initial-condition headers as the kernel define p_sigma0.
  dfloat sxyz = 0.0;
  settings.getSetting("SXYZ", sxyz);
  sigma0 = (sxyz > 0.0) ? sxyz : c / (M_PI * fmax / 2.0);

  // Source center: parse "x y z" (defaults to the origin if absent/malformed).
  srcX = srcY = srcZ = 0.0;
  if (settings.hasSetting("SOURCE POSITION")) {
    std::string srcStr;
    settings.getSetting("SOURCE POSITION", srcStr);
    std::istringstream iss(srcStr);
    iss >> srcX >> srcY >> srcZ;
  }
  if (mesh.rank == 0)
    printf("Source: Gaussian at (%g, %g, %g) m, width sigma0 = %g m "
           "(%s; fmax=%g Hz, c=%g m/s)\n",
           srcX, srcY, srcZ, sigma0,
           (sxyz > 0.0) ? "from SXYZ" : "from FMAX", fmax, c);

  dlong Nlocal = mesh.Nelements*mesh.Np*Nfields;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np*Nfields;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd"});

  /*setup trace halo exchange */
  traceHalo = mesh.HaloTraceSetup(Nfields);

  //setup timeStepper
  if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    timeStepper.Setup<TimeStepper::ab3>(mesh.Nelements,
                                        mesh.totalHaloPairs,
                                        mesh.Np, Nfields, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    timeStepper.Setup<TimeStepper::lserk4>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    timeStepper.Setup<TimeStepper::dopri5>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, platform, comm);
  }

  // set penalty parameter
  dfloat Lambda2 = 0.5;

  // compute samples of q at interpolation nodes
  q.malloc(Nlocal+Nhalo);
  o_q = platform.malloc<dfloat>(q);

  //storage for M*q during reporting
  o_Mq = platform.malloc<dfloat>(q);
  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= Nfields;

  const dfloat p_half = 1./2.;
  kernelInfo["defines/" "p_half"]= p_half;

  kernelInfo["defines/" "p_rho"]= rho;
  kernelInfo["defines/" "p_c"]= c;
  kernelInfo["defines/" "p_AcConstant"]= rho*c*c;
  kernelInfo["defines/" "p_Z_IND"]= ZFreqIndep;
  kernelInfo["defines/" "p_sigma0"]= sigma0;
  kernelInfo["defines/" "p_srcX"]= srcX;
  kernelInfo["defines/" "p_srcY"]= srcY;
  kernelInfo["defines/" "p_srcZ"]= srcZ;

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  kernelInfo["defines/" "p_Lambda2"]= Lambda2;

  // set kernel name suffix
  std::string suffix;
  if(mesh.elementType==Mesh::TRIANGLES)
    suffix = "Tri2D";
  if(mesh.elementType==Mesh::QUADRILATERALS)
    suffix = "Quad2D";
  if(mesh.elementType==Mesh::TETRAHEDRA)
    suffix = "Tet3D";
  if(mesh.elementType==Mesh::HEXAHEDRA)
    suffix = "Hex3D";

  std::string oklFilePrefix = DACOUSTICS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  // kernels from volume file
  fileName   = oklFilePrefix + "acousticsVolume" + suffix + oklFileSuffix;
  kernelName = "acousticsVolume" + suffix;

  volumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  // Safe defaults so the bc==3 branch in surface kernels compiles even when LR
  // is not active. SetupLRBC overwrites these with the correct offsets when an
  // LR vectorfit file is provided.
  kernelInfo["defines/p_LRLambda"] = 0;
  kernelInfo["defines/p_LRAlpha"]  = 0;
  kernelInfo["defines/p_LRBeta"]   = 0;
  kernelInfo["defines/p_LRA"]      = 0;
  kernelInfo["defines/p_LRB"]      = 0;
  kernelInfo["defines/p_LRC"]      = 0;
  kernelInfo["defines/p_LRYinf"]   = 0;

  // LR BC setup — may override defines above; must happen before surface kernel JIT
  SetupLRBC(kernelInfo);

  // Allocate minimal dummy LR buffers so the surface kernel never receives null
  // pointers on non-LR runs (OCCA crashes on empty device pointers).
  if (NLRPoints == 0) {
    memory<dfloat> _lr(1, 0.0);   o_LR     = platform.malloc<dfloat>(_lr);
    memory<dlong>  _li(3, 0LL);   o_LRInfo = platform.malloc<dlong> (_li);
    memory<dlong>  _ma(1, -1LL);  o_mapAcc = platform.malloc<dlong> (_ma);
    memory<dfloat> _ac(1, 0.0);
    o_acc    = platform.malloc<dfloat>(_ac);
    o_resacc = platform.malloc<dfloat>(_ac);
    o_rhsacc = platform.malloc<dfloat>(_ac);
  }

  // kernels from surface file
  fileName   = oklFilePrefix + "acousticsSurface" + suffix + oklFileSuffix;
  kernelName = "acousticsSurface" + suffix;

  surfaceKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "acousticsInitialCondition2D" + oklFileSuffix;
    kernelName = "acousticsInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "acousticsInitialCondition3D" + oklFileSuffix;
    kernelName = "acousticsInitialCondition3D";
  }

  initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                                  kernelInfo);

  // Setup receiver interpolation (builds kernel + allocates buffers)
  SetupReceivers();

  // Compute the explicit time step now so the visualization (OUTPUT INTERVAL)
  // and data-gen (DATA TEMPORAL PPW) cadences can be clamped to what the solver
  // can actually deliver. Reused verbatim in Run().
  ComputeTimeStep();

  // Setup HDF5/XDMF output writers
  SetupHDF5Output();

  // Setup ML data-generation sampler (no-op unless DATA OUTPUT is TRUE)
  SetupDataGen();
}

void acoustics_t::Logf(const char* fmt, ...)
{
  if (mesh.rank != 0) return;

  char buf[2048];
  va_list args;
  va_start(args, fmt);
  vsnprintf(buf, sizeof(buf), fmt, args);
  va_end(args);

  fputs(buf, stdout);
  fflush(stdout);

  if (logPath.empty()) return;
  std::ofstream ofs(logPath, logOpened ? std::ios::app : std::ios::trunc);
  logOpened = true;
  ofs << buf;
}

void acoustics_t::SetupHDF5Output()
{
  outDir = ".";
  if (settings.hasSetting("OUTPUT DIRECTORY"))
    settings.getSetting("OUTPUT DIRECTORY", outDir);

  // Create the output directory if it doesn't exist — every output path writes
  // here: VTU (Report), HDF5/XDMF wave field, receiver IRs, and data-gen grids.
  // error_code overload so a benign race between ranks doesn't throw.
  if (!outDir.empty()) {
    std::error_code ec;
    std::filesystem::create_directories(outDir, ec);
  }

  simulationID = "acoustics";
  if (settings.hasSetting("SIMULATION ID"))
    settings.getSetting("SIMULATION ID", simulationID);

  // Run log mirrors key diagnostics (cadence, warnings/errors) to a file so they
  // survive beyond the console. Set before any Logf() call below.
  logPath = (outDir.empty() ? std::string(".") : outDir) + "/" + simulationID + ".log";

  // Wave-field HDF5/XDMF output is exclusively a data-generation feature: it is
  // selected by OUTPUT FORMAT, which lives in the data config (DATA CONFIG FILE)
  // and is MANDATORY there. Standalone runs (no data config) can only ever write
  // the VTU visualization common to all solvers (OUTPUT TO FILE), which is
  // orthogonal — so OUTPUT FORMAT must not appear in a parent config.
  std::string dataConfig;
  if (settings.hasSetting("DATA CONFIG FILE"))
    settings.getSetting("DATA CONFIG FILE", dataConfig);
  const bool dataGen = !dataConfig.empty();

  std::string fmt = "NONE";
  if (settings.hasSetting("OUTPUT FORMAT"))
    settings.getSetting("OUTPUT FORMAT", fmt);

  if (!dataGen) {
    // No data config → no wave-field output, ever. Reject a stray OUTPUT FORMAT
    // so it can't silently look like it would write snapshots.
    LIBP_ABORT("OUTPUT FORMAT belongs in the data config (wave-field HDF5/XDMF "
               "is data-gen only); standalone runs use VTU via OUTPUT TO FILE",
               fmt != "NONE");
    outputFormat = OutputFormat::NONE;
    return;
  }

  LIBP_ABORT("OUTPUT FORMAT is mandatory in the data config — set H5COMPACT or XDMF",
             fmt != "H5COMPACT" && fmt != "XDMF");
  outputFormat = (fmt == "H5COMPACT") ? OutputFormat::H5COMPACT : OutputFormat::XDMF;

  // Precompute the vector of output times so writers can pre-allocate.
  dfloat startTime, finalTime, outputInterval;
  settings.getSetting("START TIME",      startTime);
  settings.getSetting("FINAL TIME",      finalTime);
  settings.getSetting("OUTPUT INTERVAL", outputInterval);

  // The solver cannot emit snapshots more often than one per time step. If a
  // finer OUTPUT INTERVAL is requested, clamp it to dt — otherwise the output
  // loop drifts and writes far fewer frames than the header declares, leaving
  // the XDMF/VTU series referencing snapshots that were never written.
  const dfloat effectiveInterval = std::max(outputInterval, dt);
  Logf("Visualization cadence: solver dt=%.4g s, requested OUTPUT INTERVAL="
       "%.4g s -> effective %.4g s\n", dt, outputInterval, effectiveInterval);
  if (outputInterval < dt)
    Logf("WARNING: OUTPUT INTERVAL=%.4g s could not be met (solver dt=%.4g s), "
         "clamped to %.4g s\n", outputInterval, dt, effectiveInterval);

  timeStepsOut.clear();
  dfloat t = startTime;
  while (t <= finalTime + effectiveInterval * 1e-10) {
    timeStepsOut.push_back(t);
    t += effectiveInterval;
  }

  // Construct writer only on rank 0 — all HDF5 I/O is single-rank.
  if (mesh.rank == 0) {
    if (outputFormat == OutputFormat::H5COMPACT)
      h5Writer = std::make_unique<AcousticH5CompactWriter>(*this);
    else
      h5Writer = std::make_unique<AcousticXdmfWriter>(*this);
  }
}
