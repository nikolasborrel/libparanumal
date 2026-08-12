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

#include "mesh.hpp"

// Locating arbitrary physical points in a DG mesh and building the Lagrange
// operator that evaluates a field there. Depends only on mesh_t, so it is not
// specific to the acoustics solver.
namespace pointSampling {

using namespace libp;

// Locate npts points given as xyz[3*i .. 3*i+2] (z ignored for 2D meshes).
//
// elements[i] is the local element containing point i, or -1 when the point
// lies outside this rank's partition. weights[i*mesh.Np + n] is the Lagrange
// weight of node n, zero for points that were not located.
//
// Supports Tri2D, Quad2D, Tet3D and Hex3D. Simplices invert the affine vertex
// map directly; tensor elements invert the (bi/tri)linear map by Newton.
// Points on a shared face may be located by more than one rank.
void locate(mesh_t& mesh,
            const memory<dfloat>& xyz,
            dlong npts,
            memory<dlong>& elements,
            memory<dfloat>& weights);

} // namespace pointSampling
