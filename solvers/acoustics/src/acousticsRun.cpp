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
#include <algorithm>
#include <cmath>

// LSERK4 coefficients (Carpenter & Kennedy, 1994 — 5-stage 4th-order)
static const dfloat _lserk4_rka[5] = {
   0.0,
  -567301805773.0/1357537059087.0,
  -2404267990393.0/2016746695238.0,
  -3550918686646.0/2091501179385.0,
  -1275806237668.0/842570457699.0
};
static const dfloat _lserk4_rkb[5] = {
   1432997174477.0/9575080441755.0,
   5161836677717.0/13612068292357.0,
   1720146321549.0/2090206949498.0,
   3134564353537.0/4481467310338.0,
   2277821191437.0/14882151754819.0
};
static const dfloat _lserk4_rkc[6] = {
   0.0,
   1432997174477.0/9575080441755.0,
   2526269341429.0/6820363962896.0,
   2006345519317.0/3224310063776.0,
   2802321613138.0/2924317926251.0,
   1.0
};

void acoustics_t::Run(){

  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  initialConditionKernel(mesh.Nelements,
                         startTime,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         o_q);

  // Sample the initial condition onto the data-gen source grid (no-op if off).
  SampleSourceGrid();

  // dt (member) was computed in Setup (ComputeTimeStep) so the output cadences
  // could be clamped to it; reuse the identical value here.
  timeStepper.SetTimeStep(dt);

  if(mesh.rank==0)
    printf("Time step dt = %17.15lg\n", dt);

  // Data-gen sampling cadence (dataStride) and its diagnostics are set in
  // SetupDataGen(), where dt is already known and NDataSamples is sized to it.

  if (LRNpoles > 0 || (dataGenEnabled && NDataPoints > 0)) {
    // Custom LSERK4 loop — steps both o_q and o_acc with the same coefficients,
    // and (when enabled) samples the data-gen grid on its own cadence.
    // DOPRI5 is not supported here: LR accumulators and exact-step data
    // sampling both need fixed-step control with synchronised stages.
    if (!settings.compareSetting("TIME INTEGRATOR", "LSERK4")) {
      if (mesh.rank == 0)
        printf("ERROR: LR BCs / DATA OUTPUT require TIME INTEGRATOR LSERK4\n");
      exit(1);
    }

    // the accumulators are advanced explicitly, so the fastest pole has to fit
    // inside the LSERK4 stability region or the solution diverges
    const dfloat lserk4Stability = 2.78;
    if (mesh.rank == 0 && LRMaxPole*dt > lserk4Stability) {
      dfloat cfl = 1.0;
      settings.getSetting("CFL NUMBER", cfl);
      printf("WARNING: LR pole rate %.4lg 1/s is unstable at dt %.4lg, "
             "refine the mesh or use CFL NUMBER <= %.4lg\n",
             LRMaxPole, dt, cfl*lserk4Stability/(LRMaxPole*dt));
    }

    const dlong N    = mesh.Nelements * mesh.Np * Nfields;
    const dlong Nacc = NLRPoints * LRNpoles;

    deviceMemory<dfloat> o_resq = platform.malloc<dfloat>(N);
    deviceMemory<dfloat> o_rhsq = platform.malloc<dfloat>(N);

    // Build the same LSERK4 update kernel the built-in timeStepper uses
    properties_t kInfo = platform.props();
    kInfo["defines/p_blockSize"] = 256;
    kernel_t updateKernelQ = platform.buildKernel(
        LIBP_DIR "/libs/timeStepper/okl/timeStepperLSERK4.okl", "lserk4Update", kInfo);

    dfloat outputInterval;
    settings.getSetting("OUTPUT INTERVAL", outputInterval);
    // Cannot output more often than one snapshot per step: clamp to dt so the
    // output-crossing logic below stays ahead of `time` and emits every frame
    // (matches the clamped timeStepsOut / XDMF header). See SetupHDF5Output().
    outputInterval = std::max(outputInterval, dt);
    dfloat outputTime = startTime + outputInterval;
    dfloat time = startTime;
    int tstep = 0;

    Report(time, 0);

    // Record the initial field as data-gen frame 0 (separate cadence from Report)
    if (dataGenEnabled && NDataPoints > 0)
      SampleDataGrid(0, time);

    // Allocate save buffers for the output-at-exact-time trick
    deviceMemory<dfloat> o_saveq   = platform.malloc<dfloat>(N + mesh.totalHaloPairs*mesh.Np*Nfields);
    deviceMemory<dfloat> o_saveacc = platform.malloc<dfloat>(Nacc > 0 ? Nacc : 1);
    deviceMemory<dfloat> o_saveres = platform.malloc<dfloat>(Nacc > 0 ? Nacc : 1);

    auto stepLSERK4 = [&](dfloat stepdt) {
      for (int rk = 0; rk < 5; ++rk) {
        dfloat currentTime = time + _lserk4_rkc[rk]*stepdt;
        rhsf(o_q, o_rhsq, currentTime);  // fills o_rhsq and o_rhsacc
        updateKernelQ(N, stepdt, _lserk4_rka[rk], _lserk4_rkb[rk],
                      o_rhsq, o_resq, o_q);
        if (Nacc > 0)
          updateKernelLR(NLRPoints, LRNpoles, stepdt,
                         _lserk4_rka[rk], _lserk4_rkb[rk],
                         o_rhsacc, o_resacc, o_acc);
      }
    };

    while (time < finalTime) {
      dfloat stepdt = (time + dt > finalTime) ? finalTime - time : dt;

      if (time < outputTime && time + stepdt >= outputTime) {
        o_saveq  .copyFrom(o_q,    N + mesh.totalHaloPairs*mesh.Np*Nfields);
        if (Nacc > 0) {
          o_saveacc.copyFrom(o_acc,    Nacc);
          o_saveres.copyFrom(o_resacc, Nacc);
        }
        dfloat smalldt = outputTime - time;
        stepLSERK4(smalldt);
        Report(outputTime, tstep);
        o_q.copyFrom(o_saveq, N + mesh.totalHaloPairs*mesh.Np*Nfields);
        if (Nacc > 0) {
          o_acc   .copyFrom(o_saveacc, Nacc);
          o_resacc.copyFrom(o_saveres, Nacc);
        }
        outputTime += outputInterval;
      }

      stepLSERK4(stepdt);
      time += stepdt;
      tstep++;

      // Data-gen sampling on its own step-stride cadence.
      if (dataGenEnabled && NDataPoints > 0 &&
          (tstep % dataStride == 0) && dataSampleIdx < NDataSamples)
        SampleDataGrid((int)dataSampleIdx, time);
    }
  } else {
    timeStepper.Run(*this, o_q, startTime, finalTime);
  }

  // Finalize wave-field output: write the XDMF header from the frames actually
  // written (ac.outputTimes), so the time series never references missing data.
  if (outputFormat != OutputFormat::NONE && mesh.rank == 0 && h5Writer)
    h5Writer->finalize(*this);

  // Write receiver impulse responses to HDF5 (no-op if no receivers configured)
  WriteReceiverIRs();

  // Write the ML data-generation grid to HDF5 (no-op unless DATA OUTPUT)
  WriteDataGrid();

  // output norm of the sampled receiver record
  if (NReceivers > 0) {
    dfloat recvNorm2 = 0.0;

    if (NReceiversLocal > 0) {
      o_qRecv.copyTo(qRecv);
      for (dlong n=0;n<NReceiversLocal*NRecvSamples;++n)
        recvNorm2 += qRecv[n]*qRecv[n];
    }
    mesh.comm.Allreduce(recvNorm2, Comm::Sum);

    if(mesh.rank==0)
      printf("Receiver norm = %17.15lg\n", sqrt(recvNorm2));
  }

  // output norm of final solution
  {
    //compute q.M*q
    mesh.MassMatrixApply(o_q, o_Mq);

    dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
    dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);
  }
}
