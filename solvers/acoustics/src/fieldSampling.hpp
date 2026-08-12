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

#pragma once

#include <string>
#include <vector>

// A sample set records the pressure field on a point set with its own temporal
// cadence, into its own group of <SIMULATION ID>_samples.h5. Sets are declared
// in the file named by SAMPLE SETS FILE; see ParseSampleSets for the syntax.
struct sampleSet_t {
  std::string name;                 // HDF5 group name

  // point set: a uniform Cartesian grid over the mesh bounding box
  dfloat dx     = 0.0;              // spacing [m], or derived from ppw
  dfloat ppw    = 0.0;              // points per wavelength at FMAX, if used
  dfloat jitter = 0.0;              // coordinate offset ~ U(-jitter*dx, +jitter*dx)
  int    seed   = 0;
  bool   keepOutside = false;       // keep points outside the domain as zeros,
                                    // preserving the rectilinear shape

  // temporal cadence
  enum class Cadence { INITIAL, INTERVAL } cadence = Cadence::INITIAL;
  dfloat interval = 0.0;            // [s], for Cadence::INTERVAL

  // resolved at setup
  int    shape[3] = {0, 0, 0};      // nx, ny, nz of the generated grid
  dlong  Npoints  = 0;              // points held by this rank
  int    stride   = 1;              // solver steps between samples
  dlong  Nframes  = 0;              // frames allocated
  dlong  frame    = 0;              // frames recorded so far

  memory<dfloat> xyz;               // [Npoints*3]
  memory<dlong>  elements;          // [Npoints]
  memory<dlong>  identity;          // [Npoints], indirection for the shared kernel
  memory<dfloat> values;            // host [Npoints*Nframes]
  std::vector<dfloat> times;

  deviceMemory<dfloat> o_ip;
  deviceMemory<dlong>  o_elements;
  deviceMemory<dlong>  o_identity;
  deviceMemory<dfloat> o_values;

  kernel_t kernel;
};
