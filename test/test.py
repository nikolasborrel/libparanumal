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

import os
import sys
import subprocess
import re
from pathlib import Path

#can only run the tests from LIBP_DIR/test
pwd = Path(os.getcwd())
testDir = str(pwd)
libPDir = str(pwd.parent)
solverDir = libPDir + "/solvers"

gradientDir      = solverDir + "/gradient"
advectionDir     = solverDir + "/advection"
acousticsDir     = solverDir + "/acoustics"
ellipticDir      = solverDir + "/elliptic"
fokkerPlanckDir  = solverDir + "/fokkerPlanck"
cnsDir           = solverDir + "/cns"
bnsDir           = solverDir + "/bns"
lbsDir           = solverDir + "/lbs"
insDir           = solverDir + "/ins"

gradientBin  = gradientDir      + "/gradientMain"
advectionBin = advectionDir     + "/advectionMain"
acousticsBin = acousticsDir     + "/acousticsMain"
ellipticBin  = ellipticDir      + "/ellipticMain"
fpeBin       = fokkerPlanckDir  + "/fpeMain"
cnsBin       = cnsDir           + "/cnsMain"
bnsBin       = bnsDir           + "/bnsMain"
lbsBin       = lbsDir           + "/lbsMain"
insBin       = insDir           + "/insMain"

inputRC = testDir + "/setup.rc"

TOL = 1.0e-5
alignWidth = 40

#integrators that adapt the time step; the printed dt is not constant for these
adaptiveIntegrators = ["DOPRI5", "SARK4", "SARK5"]

numeric_const_pattern = r"[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?"

if len(sys.argv)>1:
  device = sys.argv[1];
  if device!="Serial" and \
     device!="OpenMP" and \
     device!="CUDA"   and \
     device!="HIP"    and \
     device!="OpenCL":
    exit("Invalid mode requested.")
else:
  device="Serial"

class bcolors:
  TEST    = '\033[35m'
  PASS    = '\033[92m'
  WARNING = '\033[93m'
  FAIL    = '\033[91m'
  ENDC    = '\033[0m'

class setting_t:
  def __init__(self, name, value):
    self.name = name
    self.value = value

def writeSetup(filename, settings):
  str_settings=""
  for setting in settings:
    str_settings += "[" + setting.name + "]\n"
    str_settings += str(setting.value) + "\n\n"

  file = open(filename+".rc", "w")
  file.write(str_settings)
  file.close()

#HDF5 output is an opt-in build (make HDF5=1), so .h5 checks skip when the
#reader or the build is unavailable
try:
  import h5py
  haveH5py = True
except ImportError:
  haveH5py = False

def checkH5(path, datasets={}, attrs={}):
  """Verify an HDF5 output file. Returns None on success, else a message.

  datasets — {name: expected shape}, an entry of None matches any extent
  attrs    — {(dataset name, attribute name): expected value}
  """
  if not os.path.isfile(path):
    return "missing output file " + path

  with h5py.File(path, "r") as f:
    for name, shape in datasets.items():
      if name not in f:
        return "missing dataset " + name + " in " + os.path.basename(path)
      if shape is None:
        continue
      got = f[name].shape
      if len(got) != len(shape) or \
         any(e is not None and g != e for g, e in zip(got, shape)):
        return "dataset " + name + " has shape " + str(got) + \
               ", expected " + str(tuple("*" if e is None else e for e in shape))

    for (name, attr), expected in attrs.items():
      if name not in f:
        return "missing dataset " + name + " in " + os.path.basename(path)
      if attr not in f[name].attrs:
        return "missing attribute " + attr + " on " + name
      got = f[name].attrs[attr]
      if hasattr(got, "tolist"):
        got = got.tolist()
      if got != expected:
        return "attribute " + attr + " on " + name + " is " + str(got) + \
               ", expected " + str(expected)

  return None

def writeReceivers(filename, points):
  file = open(filename, "w")
  file.write(str(len(points)) + "\n")
  for point in points:
    file.write("%.16g %.16g %.16g\n" % point)
  file.close()

#a vectorfit file is a "Npoles NRealPoles NImagPoles" header followed by one
#coefficient per line in the order A, B, C, lambda, alpha, beta, Yinf, then a
#"-----" footer recording the material that was fitted
def readLRData(filename):
  lines = []
  for line in open(filename):
    if line.startswith("-----"): break
    if line.strip(): lines.append(line.strip())
  return [int(v) for v in lines[0].split()], [float(v) for v in lines[1:]]

def writeLRData(filename, header, coeffs):
  file = open(filename, "w")
  file.write("%d %d %d\n" % tuple(header))
  for coeff in coeffs:
    file.write("%.16g\n" % coeff)
  file.close()

#the residues and poles are rates in 1/s and scale with the time scale, the
#instantaneous admittance Yinf (last entry) is a ratio and does not
def scaleLRData(coeffs, kappa):
  return [kappa*coeff for coeff in coeffs[:-1]] + [coeffs[-1]]

#one real pole and one complex pair, all residues zero, so the admittance is the
#constant Yinf while both pole loops still run
def constantAdmittanceLRData(Yinf):
  return [3,1,1], [0.0, 0.0, 0.0, 100.0, 100.0, 200.0, Yinf]

def test(name, cmd, settings, referenceNorm, ranks=1, referenceDt=None,
         referenceRecvNorm=None, expectAbort=None, h5Checks=None):

  #referenceDt is None for integrators with a variable sized dt
  #expectAbort is a message the run must fail with, for configurations the
  #solver is required to reject rather than silently reinterpret
  #h5Checks is a list of (path, datasets, attrs) tuples passed to checkH5 once
  #the run has passed its norm checks; see checkH5 for the dict formats

  #create input file
  writeSetup("setup",settings)

  #print test name
  print(bcolors.TEST + f"{name:.<{alignWidth}}" + bcolors.ENDC, end="", flush=True)

  #run test
  run = subprocess.run(["mpirun", "--oversubscribe", "-np", str(ranks), cmd, inputRC],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  if expectAbort is not None:
    output = run.stdout.decode() + run.stderr.decode()
    if run.returncode == 0:
      print(bcolors.FAIL + "FAIL" + bcolors.ENDC)
      print(bcolors.WARNING + "Expected the run to be rejected, but it succeeded" + bcolors.ENDC)
      writeSetup(name,settings)
      return 1
    if expectAbort not in output:
      print(bcolors.FAIL + "FAIL" + bcolors.ENDC)
      print(bcolors.WARNING + "Expected message: " + expectAbort + bcolors.ENDC)
      writeSetup(name,settings)
      return 1
    print(bcolors.PASS + "PASS" + bcolors.ENDC)
    return 0

  if len(run.stdout.decode().splitlines())==0:
    #this failure is bad, dump the whole output for debug
    print(bcolors.FAIL + "FAIL" + bcolors.ENDC)
    print(bcolors.WARNING + name + " stdout:" + bcolors.ENDC)
    print(run.stdout.decode())
    print(bcolors.WARNING + name + " stderr:" + bcolors.ENDC)
    print(run.stderr.decode())
    #save the setup for reproducibility
    writeSetup(name,settings)
    failed = 1
  else:
    #collect last line of output
    output = run.stdout.decode().splitlines()[-1]

    #check last line's syntax
    failed=0;
    if "Solution norm = " in output:
      norm = float(output.split()[3])

      #the printed dt is only constant, and thus a meaningful reference, for
      #fixed step integrators
      integrator = next((s.value for s in settings if s.name=="TIME INTEGRATOR"), None)
      if integrator in adaptiveIntegrators:
        referenceDt = None

      #collect the time step
      dt=None
      if referenceDt is not None:
        for line in run.stdout.decode().splitlines():
          if "Time step dt = " in line:
            dt = float(line.split()[4])

      recvNorm=None
      if referenceRecvNorm is not None:
        for line in run.stdout.decode().splitlines():
          if "Receiver norm = " in line:
            recvNorm = float(line.split()[3])

      if abs(norm - referenceNorm) >= TOL:
        #failed residual check
        print(bcolors.FAIL + "FAIL" + bcolors.ENDC)
        print(bcolors.WARNING + "Expected Result: " + str(referenceNorm) + bcolors.ENDC)
        print(bcolors.WARNING + "Observed Result: " + str(norm) + bcolors.ENDC)
        #save the setup for reproducibility
        writeSetup(name,settings)
        failed = 1
      elif referenceDt is not None and abs(dt - referenceDt) >= TOL:
        #failed time step check
        print(bcolors.FAIL + "FAIL" + bcolors.ENDC)
        print(bcolors.WARNING + "Expected dt: " + str(referenceDt) + bcolors.ENDC)
        print(bcolors.WARNING + "Observed dt: " + str(dt) + bcolors.ENDC)
        #save the setup for reproducibility
        writeSetup(name,settings)
        failed = 1
      elif referenceRecvNorm is not None and abs(recvNorm - referenceRecvNorm) >= TOL:
        #failed receiver check
        print(bcolors.FAIL + "FAIL" + bcolors.ENDC)
        print(bcolors.WARNING + "Expected receiver norm: " + str(referenceRecvNorm) + bcolors.ENDC)
        print(bcolors.WARNING + "Observed receiver norm: " + str(recvNorm) + bcolors.ENDC)
        #save the setup for reproducibility
        writeSetup(name,settings)
        failed = 1
      else:
        #the norm checks passed; now inspect any requested HDF5 output
        noHDF5 = "built without HDF5" in run.stdout.decode()
        h5Error = None
        if h5Checks and haveH5py and not noHDF5:
          for path, datasets, attrs in h5Checks:
            h5Error = checkH5(path, datasets, attrs)
            if h5Error is not None:
              break

        if h5Error is not None:
          print(bcolors.FAIL + "FAIL" + bcolors.ENDC)
          print(bcolors.WARNING + "HDF5 output check: " + h5Error + bcolors.ENDC)
          #save the setup for reproducibility
          writeSetup(name,settings)
          failed = 1
        elif h5Checks and (not haveH5py or noHDF5):
          reason = "no h5py" if not haveH5py else "built without HDF5"
          print(bcolors.PASS + "PASS" + bcolors.ENDC +
                bcolors.WARNING + " (HDF5 checks skipped: " + reason + ")" + bcolors.ENDC)
        else:
          print(bcolors.PASS + "PASS" + bcolors.ENDC)
    else:
      #this failure is worse, so dump the whole output for debug
      print(bcolors.FAIL + "FAIL" + bcolors.ENDC)
      print(bcolors.WARNING + name + " stdout:" + bcolors.ENDC)
      print(run.stdout.decode())
      print(bcolors.WARNING + name + " stderr:" + bcolors.ENDC)
      print(run.stderr.decode())
      #save the setup for reproducibility
      writeSetup(name,settings)
      failed = 1

  # writeSetup(name,settings)
  # print(bcolors.WARNING + name + " stdout:" + bcolors.ENDC)
  # print(run.stdout.decode())
  # print(bcolors.WARNING + name + " stderr:" + bcolors.ENDC)
  # print(run.stderr.decode())

  #clean up
  os.remove(inputRC)

  return failed

if __name__ == "__main__":
  import testMesh
  import testGradient
  import testAdvection
  import testAcoustics
  import testElliptic
  import testFokkerPlanck
  import testCns
  import testBns
  import testLbs
  import testIns
  import testTimeStepper
  import testLinearSolver
  import testParAlmond
  import testParAdogs
  import testInitialGuess

  failCount=0;
  failCount+=testMesh.main()
  failCount+=testParAdogs.main()
  failCount+=testGradient.main()
  failCount+=testAdvection.main()
  failCount+=testAcoustics.main()
  failCount+=testElliptic.main()
  failCount+=testFokkerPlanck.main()
  failCount+=testCns.main()
  failCount+=testBns.main()
  failCount+=testLbs.main()
  failCount+=testIns.main()
  failCount+=testInitialGuess.main()
  failCount+=testTimeStepper.main()
  failCount+=testLinearSolver.main()
  failCount+=testParAlmond.main()

  sys.exit(failCount)
