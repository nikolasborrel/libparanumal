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

  dfloat cfl=1.0;
  settings.getSetting("CFL NUMBER", cfl);

  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat vmax = MaxWaveSpeed();
  dfloat dt = cfl*hmin/(vmax*(mesh.N+1.)*(mesh.N+1.));

  if(mesh.rank==0)
    printf("Time step dt = %17.15lg\n", dt);

  if (useEIRK4 || LRNpoles > 0) {
    // Custom fixed-step loop — steps both o_q and o_acc. SetupLRBC has already
    // rejected any integrator that cannot co-advance the accumulators.

    // LSERK4 advances the accumulators explicitly, so the fastest pole has to
    // fit inside its stability region or the solution diverges. This is a bound
    // on the material fit, not on the mesh, so refining does not always help —
    // EIRK4 treats the accumulators implicitly and removes it entirely.
    const dfloat lserk4Stability = 2.78;
    LIBP_ABORT("LR pole rate " << LRMaxPole << " 1/s exceeds the LSERK4 stability "
               "bound at dt " << dt << ". Use TIME INTEGRATOR EIRK4, or CFL NUMBER <= "
               << cfl*lserk4Stability/(LRMaxPole*dt),
               !useEIRK4 && LRMaxPole*dt > lserk4Stability);

    const dlong N    = mesh.Nelements * mesh.Np * Nfields;
    const dlong Nacc = NLRPoints * LRNpoles;

    deviceMemory<dfloat> o_resq, o_rhsq;
    kernel_t updateKernelQ;

    if (useEIRK4) {
      SetupEIRK4();
    } else {
      o_resq = platform.malloc<dfloat>(N);
      o_rhsq = platform.malloc<dfloat>(N);

      // Build the same LSERK4 update kernel the built-in timeStepper uses
      properties_t kInfo = platform.props();
      kInfo["defines/p_blockSize"] = 256;
      updateKernelQ = platform.buildKernel(
          LIBP_DIR "/libs/timeStepper/okl/timeStepperLSERK4.okl", "lserk4Update", kInfo);
    }

    dfloat outputInterval;
    settings.getSetting("OUTPUT INTERVAL", outputInterval);
    dfloat outputTime = startTime + outputInterval;
    dfloat time = startTime;
    int tstep = 0;

    Report(time, 0);

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

    auto step = [&](dfloat stepdt) {
      if (useEIRK4) StepEIRK4(time, stepdt);
      else          stepLSERK4(stepdt);
    };

    while (time < finalTime) {
      dfloat stepdt = (time + dt > finalTime) ? finalTime - time : dt;

      if (time < outputTime && time + stepdt >= outputTime) {
        // EIRK4 needs no residual register saved: its stage buffers are written
        // in full at stage 1 from o_q and o_acc, and carry nothing across steps
        o_saveq  .copyFrom(o_q,    N + mesh.totalHaloPairs*mesh.Np*Nfields);
        if (Nacc > 0) {
          o_saveacc.copyFrom(o_acc,    Nacc);
          if (!useEIRK4) o_saveres.copyFrom(o_resacc, Nacc);
        }
        dfloat smalldt = outputTime - time;
        step(smalldt);
        Report(outputTime, tstep);
        o_q.copyFrom(o_saveq, N + mesh.totalHaloPairs*mesh.Np*Nfields);
        if (Nacc > 0) {
          o_acc   .copyFrom(o_saveacc, Nacc);
          if (!useEIRK4) o_resacc.copyFrom(o_saveres, Nacc);
        }
        outputTime += outputInterval;
      }

      step(stepdt);
      time += stepdt;
      tstep++;
    }
  } else {
    timeStepper.SetTimeStep(dt);
    timeStepper.Run(*this, o_q, startTime, finalTime);
  }

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
