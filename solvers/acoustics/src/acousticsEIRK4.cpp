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

// Kennedy & Carpenter ARK4(3)6L[2]SA — a 6-stage explicit tableau additively
// coupled to a 6-stage L-stable ESDIRK with gamma = 1/4, sharing b and c.
//
// Both halves satisfy the order conditions through order 4; both row sums equal
// c; sum(b) = 1. The ESDIRK is stiffly accurate (its last row equals b), which
// StepEIRK4 relies on. The explicit half is not, so there is no FSAL and a step
// costs 6 rhs evaluations.
//
// erkc == esdirkc is load-bearing: one surface-kernel call evaluates the
// explicit and implicit right-hand sides at a single stage time.
//
// The embedded error weights of the pair are deliberately absent — this is a
// fixed-step scheme, and the adaptive variant in the DTU fork this was ported
// from is a non-functional stub.
static const dfloat _ark4_erkA[36] = {
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  1.0/2.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  13861.0/62500.0, 6889.0/62500.0, 0.0, 0.0, 0.0, 0.0,
  -116923316275.0/2393684061468.0, -2731218467317.0/15368042101831.0,
   9408046702089.0/11113171139209.0, 0.0, 0.0, 0.0,
  -451086348788.0/2902428689909.0, -2682348792572.0/7519795681897.0,
   12662868775082.0/11960479115383.0, 3355817975965.0/11060851509271.0, 0.0, 0.0,
   647845179188.0/3216320057751.0, 73281519250.0/8382639484533.0,
   552539513391.0/3454668386233.0, 3354512671639.0/8306763924573.0, 4040.0/17871.0, 0.0
};

static const dfloat _ark4_esdirkA[36] = {
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  1.0/4.0, 1.0/4.0, 0.0, 0.0, 0.0, 0.0,
  8611.0/62500.0, -1743.0/31250.0, 1.0/4.0, 0.0, 0.0, 0.0,
  5012029.0/34652500.0, -654441.0/2922500.0, 174375.0/388108.0, 1.0/4.0, 0.0, 0.0,
  15267082809.0/155376265600.0, -71443401.0/120774400.0, 730878875.0/902184768.0,
   2285395.0/8070912.0, 1.0/4.0, 0.0,
  82889.0/524892.0, 0.0, 15625.0/83664.0, 69875.0/102672.0, -2260.0/8211.0, 1.0/4.0
};

static const dfloat _ark4_b[6] = {
  82889.0/524892.0, 0.0, 15625.0/83664.0, 69875.0/102672.0, -2260.0/8211.0, 1.0/4.0
};

static const dfloat _ark4_c[6] = {
  0.0, 1.0/2.0, 83.0/250.0, 31.0/50.0, 17.0/20.0, 1.0
};

static const int _ark4_Nstages = 6;

void acoustics_t::SetupEIRK4() {

  const dlong N     = mesh.Nelements*mesh.Np*Nfields;
  const dlong Nhalo = mesh.totalHaloPairs*mesh.Np*Nfields;
  const dlong Nacc  = NLRPoints*LRNpoles;

  // the stage field is handed to rhsf, which exchanges the trace halo into it
  o_eirkResq = platform.malloc<dfloat>(N + Nhalo);
  o_eirkKq   = platform.malloc<dfloat>(_ark4_Nstages*N);
  o_eirkKacc = platform.malloc<dfloat>(std::max((dlong)_ark4_Nstages*Nacc, (dlong)1));
  o_eirkXacc = platform.malloc<dfloat>(std::max(Nacc, (dlong)1));

  memory<dfloat> erkA(36); for (int i=0;i<36;++i) erkA[i] = _ark4_erkA[i];
  memory<dfloat> esdA(36); for (int i=0;i<36;++i) esdA[i] = _ark4_esdirkA[i];
  memory<dfloat> b(6);     for (int i=0;i<6;++i)  b[i]    = _ark4_b[i];

  o_erkA    = platform.malloc<dfloat>(erkA);
  o_esdirkA = platform.malloc<dfloat>(esdA);
  o_erkB    = platform.malloc<dfloat>(b);

  // collective; reached by every rank, and needs no LR defines
  updateKernelEIRK4 = platform.buildKernel(DACOUSTICS "okl/acousticsEIRK4Update.okl",
                                           "acousticsEIRK4Update", mesh.props);
}

// One EIRK4 step. The explicit field update runs before the implicit
// accumulator solve at every stage: the ERK stage value depends only on
// derivatives already computed, so the wall pressure the implicit stage needs is
// known by then and the coupled system decouples into independent per-pole
// blocks. That ordering is a correctness requirement, not an optimisation.
void acoustics_t::StepEIRK4(dfloat time, dfloat stepdt) {

  const dlong N    = mesh.Nelements*mesh.Np*Nfields;
  const dlong Nacc = NLRPoints*LRNpoles;

  for (int s = 0; s < _ark4_Nstages; ++s) {

    // stage 1 reads the step-start state; both tableaux have a zero first row,
    // so Q_1 = q_n and X_1 = acc_n with nothing to compute
    deviceMemory<dfloat>& o_Qs   = (s==0) ? o_q   : o_eirkResq;
    deviceMemory<dfloat>& o_Xs   = (s==0) ? o_acc : o_eirkXacc;

    deviceMemory<dfloat>  o_kq   = o_eirkKq   + (ptrdiff_t)s*N;
    deviceMemory<dfloat>  o_kacc = o_eirkKacc + (ptrdiff_t)s*Nacc;

    // one evaluation fills both stage derivatives: the surface kernel writes the
    // accumulator rhs while computing the impedance flux
    rhsf(o_Qs, o_kq, o_Xs, o_kacc, time + _ark4_c[s]*stepdt);

    updateKernelEIRK4(N, s, stepdt, o_erkA, o_erkB, o_eirkKq, o_q, o_eirkResq);

    if (Nacc > 0) {
      if (s < _ark4_Nstages-1)
        updateKernelEIRK4LR(NLRPoints, LRNpoles, LRNRealPoles, s, stepdt,
                            o_esdirkA, o_mapAccToQ, o_LR, o_eirkKacc,
                            o_eirkResq, o_acc, o_eirkXacc);
      else
        // stiffly accurate: the final stage value is already the new solution
        o_acc.copyFrom(o_eirkXacc, Nacc);
    }
  }
}
