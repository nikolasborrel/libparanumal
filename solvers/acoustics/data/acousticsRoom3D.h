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

// Room acoustics boundary conditions (3D)
//
// BC types:
//   1 — Perfect reflection (specular)
//   2 — Frequency-independent impedance (Z = p_Z_IND)
//   3 — Locally reacting (LR) — handled in kernel body (PR 5)
//   4 — Extended reaction (ER) — handled in kernel body (PR 6)
//
// Activates upwindBC flux path in surface kernels via p_ROOM_ACOUSTICS.
// p_Z_IND is set as a kernel define from the FREQINDEP IMPEDANCE setting.

#define p_ROOM_ACOUSTICS 1

// All BCs: set ghost state to specular reflection.
// Impedance velocity (vn) is computed separately in the kernel body
// using the BC type, so the macro only needs to set the plus state.
#define acousticsDirichletConditions3D(bc, t, x, y, z, nx, ny, nz, rM, uM, vM, wM, rB, uB, vB, wB) \
{                                                                                                     \
  if(bc > 0) {                                                                                       \
    *(rB) = rM;                                                                                      \
    *(uB) = -uM;                                                                                     \
    *(vB) = -vM;                                                                                     \
    *(wB) = -wM;                                                                                     \
  }                                                                                                  \
}

// Initial condition: Gaussian pressure pulse at the origin
#define acousticsInitialConditions3D(t, x, y, z, r, u, v, w) \
{                                                              \
  *(r) = exp(-3.0*(x*x + y*y + z*z));                         \
  *(u) = 0.0;                                                  \
  *(v) = 0.0;                                                  \
  *(w) = 0.0;                                                  \
}
