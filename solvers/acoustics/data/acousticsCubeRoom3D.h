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

// Room acoustics BCs for a unit cube mesh [0,1]^3 (e.g. cube_500hz_p4_5ppw.msh).
// Matches the DTU test setup: source centered at (0.5, 0.5, 0.5), width sigma=0.4.
// Use for cross-repo comparison against libparanumal-dtu reference WAV files.

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

// Gaussian pulse centered at (0.5, 0.5, 0.5), sigma=0.4 — matches DTU SXYZ=0.4
#define acousticsInitialConditions3D(t, x, y, z, r, u, v, w)  \
{                                                               \
  const dfloat _sx = x - 0.5;                                  \
  const dfloat _sy = y - 0.5;                                  \
  const dfloat _sz = z - 0.5;                                  \
  const dfloat _sig = 0.4;                                      \
  *(r) = exp(-(_sx*_sx + _sy*_sy + _sz*_sz) / (_sig*_sig));   \
  *(u) = 0.0;                                                   \
  *(v) = 0.0;                                                   \
  *(w) = 0.0;                                                   \
}
