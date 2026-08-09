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
data2DRoom = acousticsDir + "/data/acousticsRoom2D.h"
data3DRoom = acousticsDir + "/data/acousticsRoom3D.h"

def acousticsSettings(rcformat="2.0", data_file=data2D,
                     mesh="BOX", dim=2, element=4, nx=10, ny=10, nz=10, boundary_flag=-1,
                     box_dim=10,
                     degree=4, thread_model=device, platform_number=0, device_number=0,
                      time_integrator="DOPRI5", cfl=1.0, start_time=0.0, final_time=1.0,
                      output_to_file="FALSE", density=1.0, sound_speed=1.0,
                      impedance=415.0, surface_flux="UPWIND"):
  return [setting_t("FORMAT", rcformat),
          setting_t("DATA FILE", data_file),
          setting_t("DENSITY", density),
          setting_t("SPEED OF SOUND", sound_speed),
          setting_t("FREQINDEP IMPEDANCE", impedance),
          setting_t("SURFACE FLUX", surface_flux),
          setting_t("MESH FILE", mesh),
          setting_t("MESH DIMENSION", dim),
          setting_t("ELEMENT TYPE", element),
          setting_t("BOX DIMX", box_dim),
          setting_t("BOX DIMY", box_dim),
          setting_t("BOX DIMZ", box_dim),
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
                    referenceNorm=31.6551427875736,
                    referenceDt=0.000229059533912066)

  failCount += test(name="testAcousticsCentralTri",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2D,dim=2,
                                               surface_flux="CENTRAL"),
                    referenceNorm=10.1302597729239)

  failCount += test(name="testAcousticsCentralHex",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=12,data_file=data3D,dim=3,
                                               degree=2,surface_flux="CENTRAL"),
                    referenceNorm=31.6581825895929)

  #boundary flag 1 is a rigid wall, flag 2 the frequency-independent impedance
  #BC. The pulse sits at the centre of a 2m box, so it hits the walls at t=1
  failCount += test(name="testAcousticsRoomRigidTri",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=1,final_time=2.0),
                    referenceNorm=0.723214962153181)

  failCount += test(name="testAcousticsRoomRigidTet",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=6,data_file=data3DRoom,dim=3,
                                               degree=2,box_dim=2,boundary_flag=1,
                                               final_time=2.0),
                    referenceNorm=0.614784622759662)

  failCount += test(name="testAcousticsRoomImpedanceTri",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=2,final_time=2.0,
                                               impedance=5.0),
                    referenceNorm=0.387485800176955)

  failCount += test(name="testAcousticsRoomImpedanceQuad",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=4,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=2,final_time=2.0,
                                               impedance=5.0),
                    referenceNorm=0.387479993745552)

  failCount += test(name="testAcousticsRoomImpedanceTet",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=6,data_file=data3DRoom,dim=3,
                                               degree=2,box_dim=2,boundary_flag=2,
                                               final_time=2.0,impedance=5.0),
                    referenceNorm=0.298868490129398)

  failCount += test(name="testAcousticsRoomImpedanceHex",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=12,data_file=data3DRoom,dim=3,
                                               degree=2,box_dim=2,boundary_flag=2,
                                               final_time=2.0,impedance=5.0),
                    referenceNorm=0.298590700093030)

  #Z=rho*c is anechoic
  failCount += test(name="testAcousticsRoomAnechoicTri",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=2,final_time=2.0,
                                               impedance=1.0),
                    referenceNorm=0.133067270222819)

  #Z>>rho*c must recover the rigid wall norm above
  failCount += test(name="testAcousticsRoomRigidLimitTri",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=2,final_time=2.0,
                                               impedance=1.0e8),
                    referenceNorm=0.723214962153181)

  failCount += test(name="testAcousticsRoomCentralTri",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=2,final_time=2.0,
                                               impedance=5.0,surface_flux="CENTRAL"),
                    referenceNorm=0.387483415651966)

  failCount += test(name="testAcousticsRoomCentralHex",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=12,data_file=data3DRoom,dim=3,
                                               degree=2,box_dim=2,boundary_flag=2,
                                               final_time=2.0,impedance=5.0,
                                               surface_flux="CENTRAL"),
                    referenceNorm=0.298738225769498)

  #rho*c=1 leaves the velocity unscaled, so c=2 at t=1 must reproduce the
  #c=1 norm at t=2
  failCount += test(name="testAcousticsRoomScaledTri",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=2,final_time=1.0,
                                               impedance=5.0,density=0.5,sound_speed=2.0),
                    referenceNorm=0.387485800176955)

  failCount += test(name="testAcousticsRoomImpedanceTri_MPI", ranks=4,
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=2,final_time=2.0,
                                               impedance=5.0),
                    referenceNorm=0.387475178751856)

  #clean up
  for file_name in os.listdir(testDir):
    if file_name.endswith('.vtu'):
      os.remove(testDir + "/" + file_name)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)
