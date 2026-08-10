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

// Room acoustics with an FMAX-parameterized Gaussian source (3D).
//
// The pulse width is the kernel define p_sigma0, set in Setup() from [SXYZ] if
// given, else from [FMAX] as 2c/(pi*FMAX). BC types as in acousticsRoom3D.h.

#define p_ROOM_ACOUSTICS 1

#define acousticsDirichletConditions3D(bc, t, x, y, z, nx, ny, nz, rM, uM, vM, wM, rB, uB, vB, wB) \
{                                                                                                     \
  if(bc > 0) {                                                                                       \
    *(rB) = rM;                                                                                      \
    *(uB) = -uM;                                                                                     \
    *(vB) = -vM;                                                                                     \
    *(wB) = -wM;                                                                                     \
  }                                                                                                  \
}

#define acousticsInitialConditions3D(t, x, y, z, r, u, v, w)              \
{                                                                        \
  const dfloat _dx = x - (p_srcX);                                       \
  const dfloat _dy = y - (p_srcY);                                       \
  const dfloat _dz = z - (p_srcZ);                                       \
  *(r) = exp(-(_dx*_dx + _dy*_dy + _dz*_dz) / (p_sigma0*p_sigma0));       \
  *(u) = 0.0;                                                            \
  *(v) = 0.0;                                                            \
  *(w) = 0.0;                                                            \
}
