#!/usr/bin/env python3

#####################################################################################
#
#The MIT License (MIT)
#
#Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#
#####################################################################################

from test import *

data2D = acousticsDir + "/data/acousticsGaussian2D.h"
data3D = acousticsDir + "/data/acousticsGaussian3D.h"

def acousticsSettings(rcformat="2.0", data_file=data2D,
                     mesh="BOX", dim=2, element=4, nx=10, ny=10, nz=10, boundary_flag=-1,
                     degree=4, thread_model=device, platform_number=0, device_number=0,
                      time_integrator="DOPRI5", cfl=1.0, start_time=0.0, final_time=1.0,
                      output_to_file="FALSE", density=1.0, sound_speed=1.0):
  return [setting_t("FORMAT", rcformat),
          setting_t("DATA FILE", data_file),
          setting_t("DENSITY", density),
          setting_t("SPEED OF SOUND", sound_speed),
          setting_t("MESH FILE", mesh),
          setting_t("MESH DIMENSION", dim),
          setting_t("ELEMENT TYPE", element),
          setting_t("BOX NX", nx),
          setting_t("BOX NY", ny),
          setting_t("BOX NZ", nz),
          setting_t("BOX BOUNDARY FLAG", boundary_flag),
          setting_t("POLYNOMIAL DEGREE", degree),
          setting_t("THREAD MODEL", thread_model),
          setting_t("PLATFORM NUMBER", platform_number),
          setting_t("DEVICE NUMBER", device_number),
          setting_t("TIME INTEGRATOR", time_integrator),
          setting_t("CFL NUMBER", cfl),
          setting_t("START TIME", start_time),
          setting_t("FINAL TIME", final_time),
          setting_t("OUTPUT TO FILE", output_to_file)]

def main():
  failCount=0;

  failCount += test(name="testAcousticsTri",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2D,dim=2),
                    referenceNorm=10.1302322430996)

  failCount += test(name="testAcousticsQuad",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=4,data_file=data2D,dim=2),
                    referenceNorm=10.1299609797959)

  failCount += test(name="testAcousticsTet",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=6,data_file=data3D,dim=3,
                                               degree=2),
                    referenceNorm=31.6577046152384)

  failCount += test(name="testAcousticsHex",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=12,data_file=data3D,dim=3,
                                               degree=2),
                    referenceNorm=31.6576028812776)

  failCount += test(name="testAcousticsTri_MPI", ranks=4,
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2D,dim=2,output_to_file="TRUE"),
                    referenceNorm=10.1300558638317)

  #fixed step LSERK4 so dt is constant; dt = cfl*hmin/(c*(N+1)^2) scales as 1/c,
  #so referenceDt is the c=1 Tet dt / 343
  failCount += test(name="testAcousticsMedium",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=6,data_file=data3D,dim=3,
                                               degree=2,sound_speed=343.0,density=1.2,
                                               time_integrator="LSERK4",final_time=0.05),
                    referenceNorm=31.6579637397127,
                    referenceDt=0.000229059533912066)

  #clean up
  for file_name in os.listdir(testDir):
    if file_name.endswith('.vtu'):
      os.remove(testDir + "/" + file_name)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)
