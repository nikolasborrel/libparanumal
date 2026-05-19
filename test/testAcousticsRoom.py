#!/usr/bin/env python3

#####################################################################################
#
#The MIT License (MIT)
#
#Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus
#
#Permission is hereby granted, free of charge, to any person without restriction, including
#the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software...
#
#####################################################################################

# Room-acoustics regression tests for PR-series acoustics features:
#   testAcousticsRoomPerfRefl  — BOX, BC 1 (perfect reflector), LSERK4
#   testAcousticsRoomFreqIndep — BOX, BC 2 (freq-independent impedance), LSERK4
#   testAcousticsRoomLR        — BOX, BC 3 (locally-reacting LR), LSERK4
#
# All three use:
#   - Tet BOX mesh [-1,1]^3 (NX=NY=NZ=2, p=4)
#   - acousticsRoom3D.h (Gaussian source at origin, room BCs)
#   - FINAL TIME = 0.05 s  (> first reflection ≈ 0.006 s)
#   - DENSITY=1.2 kg/m³, SPEED OF SOUND=343 m/s (air)
#
# The ordering norm(perf_refl) > norm(freq_indep) confirms energy absorption at walls.
# norm(LR) differs from norm(freq_indep) due to frequency-dependent admittance.

from test import *

roomData3D    = acousticsDir + "/data/acousticsRoom3D.h"
lrVectfitFile = acousticsDir + "/data/freq_dep_lr.dat"

def roomAcousticsSettings(boundary_flag=1,
                          lr_vectorfit_file=None,
                          freqindep_impedance=415.0,
                          final_time=0.05,
                          cfl=1.0):
  settings = [
    setting_t("FORMAT",            "2.0"),
    setting_t("DATA FILE",         roomData3D),
    setting_t("MESH FILE",         "BOX"),
    setting_t("MESH DIMENSION",    3),
    setting_t("ELEMENT TYPE",      6),       # tet
    setting_t("BOX DIMX",          2),
    setting_t("BOX DIMY",          2),
    setting_t("BOX DIMZ",          2),
    setting_t("BOX NX",            2),
    setting_t("BOX NY",            2),
    setting_t("BOX NZ",            2),
    setting_t("BOX BOUNDARY FLAG", boundary_flag),
    setting_t("POLYNOMIAL DEGREE", 4),
    setting_t("THREAD MODEL",      device),
    setting_t("PLATFORM NUMBER",   0),
    setting_t("DEVICE NUMBER",     0),
    setting_t("TIME INTEGRATOR",   "LSERK4"),
    setting_t("CFL NUMBER",        cfl),
    setting_t("START TIME",        0.0),
    setting_t("FINAL TIME",        final_time),
    setting_t("OUTPUT INTERVAL",   0.01),
    setting_t("OUTPUT TO FILE",    "FALSE"),
    setting_t("DENSITY",           1.2),
    setting_t("SPEED OF SOUND",    343.0),
    setting_t("FREQINDEP IMPEDANCE", freqindep_impedance),
  ]
  if lr_vectorfit_file:
    settings.append(setting_t("LR VECTORFIT FILE", lr_vectorfit_file))
  return settings


def main():
  failCount = 0

  failCount += test(name="testAcousticsRoomPerfRefl",
                    cmd=acousticsBin,
                    settings=roomAcousticsSettings(boundary_flag=1),
                    referenceNorm=0.597430370725016)

  failCount += test(name="testAcousticsRoomFreqIndep",
                    cmd=acousticsBin,
                    settings=roomAcousticsSettings(boundary_flag=2,
                                                   freqindep_impedance=415.0),
                    referenceNorm=0.596889525574821)

  # LR BCs use CFL=0.1: the LSERK4 explicit integrator for the accumulator ODE
  # requires |lambda_max * dt| < ~2.8 (stability region). The vectorfit data has
  # complex poles with |lambda| ~ 180,000 s^-1; the coarse 2-element BOX mesh
  # gives dt ~ 5.8e-5 s at CFL=1, so CFL <= 0.1 is needed for stability.
  # Production meshes (cm-scale elements) satisfy this automatically.
  failCount += test(name="testAcousticsRoomLR",
                    cmd=acousticsBin,
                    settings=roomAcousticsSettings(boundary_flag=3,
                                                   lr_vectorfit_file=lrVectfitFile,
                                                   cfl=0.1),
                    referenceNorm=0.597434850050192)

  return failCount

if __name__ == "__main__":
  failCount = 0
  failCount += main()
  sys.exit(failCount)
