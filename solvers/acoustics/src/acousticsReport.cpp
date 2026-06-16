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

void acoustics_t::Report(dfloat time, int tstep){

  static int frame=0;

  // Sample receivers at this output instant
  if (NReceiversLocal > 0 && recvSampleIdx < NRecvSamples) {
    receiverKernel(NReceiversLocal, recvSampleIdx,
                   o_qRecv, o_recvElements, o_recvElementsIdx,
                   o_recvIP, o_q);
    recvSampleIdx++;
  }

  //compute q.M*q
  mesh.MassMatrixApply(o_q, o_Mq);

  dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
  dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm));

  if(mesh.rank==0)
    printf("%5.2f (%d), %5.2f (time, timestep, norm)\n", time, tstep, norm2);

  if (outputFormat != OutputFormat::VTU) {
    // HDF5/XDMF: rank 0 writes; writer handles its own o_q.copyTo(q)
    if (mesh.rank == 0 && h5Writer)
      h5Writer->write(*this, frame);
    frame++;
  } else if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {
    // copy data back to host
    o_q.copyTo(q);

    // output field files (honor OUTPUT DIRECTORY, like the HDF5/XDMF writers)
    std::string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    const std::string& dir = outDir.empty() ? std::string(".") : outDir;
    char fname[BUFSIZ];
    snprintf(fname, sizeof(fname), "%s/%s_%04d_%04d.vtu",
             dir.c_str(), name.c_str(), mesh.rank, frame++);

    PlotFields(q, std::string(fname));
  }
}
