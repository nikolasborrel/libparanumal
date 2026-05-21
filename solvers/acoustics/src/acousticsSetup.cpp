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

void acoustics_t::Setup(platform_t& _platform, mesh_t& _mesh,
                        acousticsSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  Nfields = (mesh.dim==3) ? 4:3;

  settings.getSetting("DENSITY", rho);
  settings.getSetting("SPEED OF SOUND", c);
  settings.getSetting("FREQINDEP IMPEDANCE", ZFreqIndep);

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

  // Setup HDF5/XDMF output writers
  SetupHDF5Output();
}

void acoustics_t::SetupHDF5Output()
{
  outDir = ".";
  if (settings.hasSetting("OUTPUT DIRECTORY"))
    settings.getSetting("OUTPUT DIRECTORY", outDir);

  simulationID = "acoustics";
  if (settings.hasSetting("SIMULATION ID"))
    settings.getSetting("SIMULATION ID", simulationID);

  std::string fmt = "VTU";
  if (settings.hasSetting("OUTPUT FORMAT"))
    settings.getSetting("OUTPUT FORMAT", fmt);

  if (fmt == "H5COMPACT")
    outputFormat = OutputFormat::H5COMPACT;
  else if (fmt == "XDMF")
    outputFormat = OutputFormat::XDMF;
  else
    outputFormat = OutputFormat::VTU;

  if (outputFormat == OutputFormat::VTU)
    return;

  // Precompute the vector of output times so writers can pre-allocate.
  dfloat startTime, finalTime, outputInterval;
  settings.getSetting("START TIME",      startTime);
  settings.getSetting("FINAL TIME",      finalTime);
  settings.getSetting("OUTPUT INTERVAL", outputInterval);

  timeStepsOut.clear();
  dfloat t = startTime;
  while (t <= finalTime + outputInterval * 1e-10) {
    timeStepsOut.push_back(t);
    t += outputInterval;
  }

  // Construct writer only on rank 0 — all HDF5 I/O is single-rank.
  if (mesh.rank == 0) {
    if (outputFormat == OutputFormat::H5COMPACT)
      h5Writer = std::make_unique<AcousticH5CompactWriter>(*this);
    else
      h5Writer = std::make_unique<AcousticXdmfWriter>(*this);
  }
}
