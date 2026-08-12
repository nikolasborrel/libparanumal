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

recvFile = testDir + "/receivers.dat"

#surface admittance of a 5cm porous layer, fitted with 14 poles over 50-2000Hz,
#see solvers/acoustics/vectorfit/README.md
lrMaterialFile = acousticsDir + "/data/LRDATA14.dat"
lrFile         = testDir + "/lrVectorfit.dat"

def acousticsSettings(rcformat="2.0", data_file=data2D,
                     mesh="BOX", dim=2, element=4, nx=10, ny=10, nz=10, boundary_flag=-1,
                     box_dim=10,
                     degree=4, thread_model=device, platform_number=0, device_number=0,
                      time_integrator="DOPRI5", cfl=1.0, start_time=0.0, final_time=1.0,
                      output_to_file="FALSE", density=1.0, sound_speed=1.0,
                      impedance=415.0, surface_flux="UPWIND",
                      receiver_file=None, lr_file=None):
  settings = [setting_t("FORMAT", rcformat),
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

  if receiver_file is not None:
    settings.append(setting_t("RECEIVER FILE", receiver_file))

  if lr_file is not None:
    settings.append(setting_t("LR VECTORFIT FILE", lr_file))

  return settings

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

  #a mesh with LR walls and no fit must be rejected: the surface kernel would
  #otherwise yield vn=0 there, silently simulating a perfectly rigid wall
  failCount += test(name="testAcousticsRoomLRMissingFit",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=3,final_time=0.1,
                                               time_integrator="LSERK4"),
                    referenceNorm=None,
                    expectAbort="Mesh has BC==3 (locally-reacting) faces but no LR VECTORFIT FILE")

  #the surface kernel walks the complex pairs two at a time, so a header whose
  #pole counts do not add up would read past the end of each accumulator row
  writeLRData(lrFile, [4,1,1], [0.0, 0.0, 0.0, 100.0, 100.0, 200.0, 0.2])
  failCount += test(name="testAcousticsRoomLRBadPoleCount",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=3,final_time=0.1,
                                               time_integrator="LSERK4",lr_file=lrFile),
                    referenceNorm=None,
                    expectAbort="LR vectorfit header inconsistent")

  #boundary flag 3 is the locally-reacting BC. A constant surface admittance Y
  #must reproduce the frequency-independent BC with Z = 1/Y
  writeLRData(lrFile, *constantAdmittanceLRData(1.0/5.0))
  failCount += test(name="testAcousticsRoomLRConstantTri",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=3,final_time=2.0,
                                               time_integrator="LSERK4",lr_file=lrFile),
                    referenceNorm=0.387485800176955)

  #the accumulator ODE is co-advanced by a fixed-step loop, so an adaptive
  #integrator must be rejected rather than silently ignored
  failCount += test(name="testAcousticsRoomLRDopri5",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=3,final_time=0.1,
                                               lr_file=lrFile),
                    referenceNorm=None,
                    expectAbort="require TIME INTEGRATOR LSERK4")

  failCount += test(name="testAcousticsRoomLRConstantQuad",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=4,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=3,final_time=2.0,
                                               time_integrator="LSERK4",lr_file=lrFile),
                    referenceNorm=0.387479993745552)

  failCount += test(name="testAcousticsRoomLRConstantTet",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=6,data_file=data3DRoom,dim=3,
                                               degree=2,box_dim=2,boundary_flag=3,
                                               final_time=2.0,time_integrator="LSERK4",
                                               lr_file=lrFile),
                    referenceNorm=0.298868490129398)

  failCount += test(name="testAcousticsRoomLRConstantHex",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=12,data_file=data3DRoom,dim=3,
                                               degree=2,box_dim=2,boundary_flag=3,
                                               final_time=2.0,time_integrator="LSERK4",
                                               lr_file=lrFile),
                    referenceNorm=0.298590700093030)

  failCount += test(name="testAcousticsRoomLRCentralTri",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=3,final_time=2.0,
                                               time_integrator="LSERK4",lr_file=lrFile,
                                               surface_flux="CENTRAL"),
                    referenceNorm=0.387483415651966)

  #Y=0 is a rigid wall
  writeLRData(lrFile, *constantAdmittanceLRData(0.0))
  failCount += test(name="testAcousticsRoomLRRigidTri",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=3,final_time=2.0,
                                               time_integrator="LSERK4",lr_file=lrFile),
                    referenceNorm=0.723214962153181)

  #a fitted material: the poles sit at 300-21000 rad/s, so the medium has to be
  #air for the pulse spectrum to reach the frequency range that was fitted
  lrHeader, lrCoeffs = readLRData(lrMaterialFile)
  writeLRData(lrFile, lrHeader, lrCoeffs)
  failCount += test(name="testAcousticsRoomLRMaterialTri",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=3,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=3,final_time=0.02,
                                               density=1.2,sound_speed=343.0,
                                               time_integrator="LSERK4",lr_file=lrFile),
                    referenceNorm=0.484349091289788)

  failCount += test(name="testAcousticsRoomLRMaterialQuad",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=4,data_file=data2DRoom,dim=2,
                                               box_dim=2,boundary_flag=3,final_time=0.02,
                                               density=1.2,sound_speed=343.0,
                                               time_integrator="LSERK4",lr_file=lrFile),
                    referenceNorm=0.484348255295611)

  failCount += test(name="testAcousticsRoomLRMaterialTet",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=6,data_file=data3DRoom,dim=3,
                                               degree=2,nx=6,ny=6,nz=6,box_dim=2,
                                               boundary_flag=3,final_time=0.02,
                                               density=1.2,sound_speed=343.0,
                                               time_integrator="LSERK4",lr_file=lrFile),
                    referenceNorm=0.350474716365334)

  failCount += test(name="testAcousticsRoomLRMaterialHex",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=12,data_file=data3DRoom,dim=3,
                                               degree=2,nx=6,ny=6,nz=6,box_dim=2,
                                               boundary_flag=3,final_time=0.02,
                                               density=1.2,sound_speed=343.0,
                                               time_integrator="LSERK4",lr_file=lrFile),
                    referenceNorm=0.352709954891316)

  #halving the time scale and moving the poles up by the same factor is an exact
  #invariance of the LR BC, so this reproduces the material norm above
  writeLRData(lrFile, lrHeader, scaleLRData(lrCoeffs, 2.0))
  failCount += test(name="testAcousticsRoomLRScaledTet",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=6,data_file=data3DRoom,dim=3,
                                               degree=2,nx=6,ny=6,nz=6,box_dim=2,
                                               boundary_flag=3,final_time=0.01,
                                               density=0.6,sound_speed=686.0,
                                               time_integrator="LSERK4",lr_file=lrFile),
                    referenceNorm=0.350474716365334)

  #the same halved time scale with the poles left in place: a constant
  #admittance would still give the material norm, a frequency-dependent one not
  writeLRData(lrFile, lrHeader, lrCoeffs)
  failCount += test(name="testAcousticsRoomLRUnscaledTet",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=6,data_file=data3DRoom,dim=3,
                                               degree=2,nx=6,ny=6,nz=6,box_dim=2,
                                               boundary_flag=3,final_time=0.01,
                                               density=0.6,sound_speed=686.0,
                                               time_integrator="LSERK4",lr_file=lrFile),
                    referenceNorm=0.264633120701544)

  failCount += test(name="testAcousticsRoomLRMaterialTet_MPI", ranks=4,
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=6,data_file=data3DRoom,dim=3,
                                               degree=2,nx=4,ny=4,nz=4,box_dim=2,
                                               boundary_flag=3,final_time=0.02,
                                               density=1.2,sound_speed=343.0,
                                               time_integrator="LSERK4",lr_file=lrFile),
                    referenceNorm=0.350247744152265)

  #the box centre is a mesh node, so the t=0 sample is exp(0)=1 exactly, and
  #final_time=0 leaves that sample as the whole record
  writeReceivers(recvFile, [(0.0, 0.0, 0.0)])
  failCount += test(name="testAcousticsReceiverNode",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=6,data_file=data3DRoom,dim=3,
                                               degree=2,box_dim=2,boundary_flag=1,
                                               time_integrator="LSERK4",final_time=0.0,
                                               receiver_file=recvFile),
                    referenceNorm=0.614876573596189,
                    referenceRecvNorm=1.0)

  #off-node receivers exercise the interpolation weights; sampling must not
  #perturb the solution, so the norm matches testAcousticsRoomRigidTet
  writeReceivers(recvFile, [(0.3, 0.15, -0.22), (-0.55, 0.42, 0.61)])
  failCount += test(name="testAcousticsReceiverTet",
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=6,data_file=data3DRoom,dim=3,
                                               degree=2,box_dim=2,boundary_flag=1,
                                               final_time=2.0,receiver_file=recvFile),
                    referenceNorm=0.614784622759662,
                    referenceRecvNorm=1.45569073560319)

  #ranks owning no receiver still take part in the collective kernel build and
  #the norm reduction
  failCount += test(name="testAcousticsReceiverTet_MPI", ranks=4,
                    cmd=acousticsBin,
                    settings=acousticsSettings(element=6,data_file=data3DRoom,dim=3,
                                               degree=2,box_dim=2,boundary_flag=1,
                                               final_time=2.0,receiver_file=recvFile),
                    referenceNorm=0.614966266247452,
                    referenceRecvNorm=1.45552027978335)

  #clean up
  os.remove(recvFile)
  os.remove(lrFile)
  for file_name in os.listdir(testDir):
    if file_name.endswith('.vtu'):
      os.remove(testDir + "/" + file_name)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)
