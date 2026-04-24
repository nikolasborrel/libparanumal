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

#include "catch_amalgamated.hpp"
#include "acoustics.hpp"      // pulls in libp namespace, Comm, and DACOUSTICS macro
#include <cstdio>             // popen / pclose / fgets / snprintf
#include <cstdlib>            // mkstemp
#include <cstring>            // strerror
#include <string>
#include <stdexcept>
#include <unistd.h>           // close

// Write a temporary .rc settings file using mkstemp (TOCTOU-safe).
// Returns the path to the created file.
static std::string writeTempRcFile(int elementType, int dim, int degree,
                                    int nx, int ny, int nz,
                                    const std::string& dataFile,
                                    int boundaryFlag) {
  char tmpPath[] = "/tmp/acousticsTestXXXXXX";
  int fd = mkstemp(tmpPath);
  if (fd == -1) {
    throw std::runtime_error(std::string("mkstemp failed: ") + strerror(errno));
  }

  std::string content;
  content += "[FORMAT]\n2.0\n\n";
  content += "[DATA FILE]\n" + dataFile + "\n\n";
  content += "[MESH FILE]\nBOX\n\n";
  content += "[MESH DIMENSION]\n" + std::to_string(dim) + "\n\n";
  content += "[ELEMENT TYPE]\n" + std::to_string(elementType) + "\n\n";
  content += "[BOX NX]\n" + std::to_string(nx) + "\n\n";
  content += "[BOX NY]\n" + std::to_string(ny) + "\n\n";
  content += "[BOX NZ]\n" + std::to_string(nz) + "\n\n";
  content += "[BOX BOUNDARY FLAG]\n" + std::to_string(boundaryFlag) + "\n\n";
  content += "[POLYNOMIAL DEGREE]\n" + std::to_string(degree) + "\n\n";
  content += "[THREAD MODEL]\nSerial\n\n";
  content += "[PLATFORM NUMBER]\n0\n\n";
  content += "[DEVICE NUMBER]\n0\n\n";
  content += "[TIME INTEGRATOR]\nDOPRI5\n\n";
  content += "[CFL NUMBER]\n1.0\n\n";
  content += "[START TIME]\n0\n\n";
  content += "[FINAL TIME]\n1.0\n\n";
  content += "[OUTPUT TO FILE]\nFALSE\n\n";

  if (write(fd, content.c_str(), content.size()) == -1) {
    close(fd);
    throw std::runtime_error(std::string("write to temp file failed: ") + strerror(errno));
  }
  close(fd);

  // Rename to add .rc extension so acousticsMain can parse the file
  std::string rcPath = std::string(tmpPath) + ".rc";
  if (rename(tmpPath, rcPath.c_str()) != 0) {
    throw std::runtime_error(std::string("rename temp file failed: ") + strerror(errno));
  }
  return rcPath;
}

// Invoke acousticsMain as a subprocess (mirrors test/test.py subprocess pattern)
// and parse the "Solution norm = " value from its stdout.
// Critical: acousticsRun.cpp uses printf() which writes to C FILE* fd 1, not std::cout.
// std::cout.rdbuf() redirect cannot capture it. Subprocess via popen is the correct approach.
static double runAndGetNorm(int elementType, int dim, int degree,
                             int nx, int ny, int nz,
                             const std::string& dataFile,
                             int boundaryFlag) {
  std::string rcPath = writeTempRcFile(elementType, dim, degree,
                                        nx, ny, nz, dataFile, boundaryFlag);

  std::string cmd = std::string("mpirun --oversubscribe -np 1 ")
                  + DACOUSTICS "acousticsMain "
                  + rcPath
                  + " 2>/dev/null";

  FILE* pipe = popen(cmd.c_str(), "r");
  if (!pipe) {
    std::remove(rcPath.c_str());
    throw std::runtime_error(std::string("popen failed: ") + strerror(errno));
  }

  char buf[4096];
  std::string matchedLine;
  while (fgets(buf, sizeof(buf), pipe)) {
    std::string line(buf);
    if (line.find("Solution norm = ") != std::string::npos) {
      matchedLine = line;
    }
  }

  int rc = pclose(pipe);
  std::remove(rcPath.c_str());

  if (rc != 0) {
    throw std::runtime_error("acousticsMain exited non-zero: rc=" + std::to_string(rc));
  }

  if (matchedLine.empty()) {
    throw std::runtime_error("Solution norm not found in output");
  }

  auto pos = matchedLine.find("Solution norm = ");
  std::string normStr = matchedLine.substr(pos + std::string("Solution norm = ").size());
  return std::stod(normStr);
}

// Unit tests: assert Gaussian pulse L2 norm within 1e-5 of reference values
// Reference norms from test/testAcoustics.py (verbatim)

TEST_CASE("Gaussian pulse L2 norm - Tri2D", "[unit][tri2d]") {
  double norm = runAndGetNorm(3, 2, 4, 10, 10, 10,
                               DACOUSTICS "data/acousticsGaussian2D.h", -1);
  REQUIRE(norm == Catch::Approx(10.1302322430996).margin(1e-5));
}

TEST_CASE("Gaussian pulse L2 norm - Quad2D", "[unit][quad2d]") {
  double norm = runAndGetNorm(4, 2, 4, 10, 10, 10,
                               DACOUSTICS "data/acousticsGaussian2D.h", -1);
  REQUIRE(norm == Catch::Approx(10.1299609797959).margin(1e-5));
}

TEST_CASE("Gaussian pulse L2 norm - Tet3D", "[unit][tet3d]") {
  double norm = runAndGetNorm(6, 3, 2, 10, 10, 10,
                               DACOUSTICS "data/acousticsGaussian3D.h", -1);
  REQUIRE(norm == Catch::Approx(31.6577046152384).margin(1e-5));
}

TEST_CASE("Gaussian pulse L2 norm - Hex3D", "[unit][hex3d]") {
  double norm = runAndGetNorm(12, 3, 2, 10, 10, 10,
                               DACOUSTICS "data/acousticsGaussian3D.h", -1);
  REQUIRE(norm == Catch::Approx(31.6576028812776).margin(1e-5));
}

// BC smoke tests: assert solver runs without throwing (crash-only, no norm assertion)
// FreqIndep and PerfRefl both use boundary_flag=1 (wall) in Phase 1.
// Norm assertions for these BC types are added in Phase 3 when the flux variants land.

TEST_CASE("BC smoke - PerfRefl Tri2D", "[integration][perf_refl][tri2d]") {
  REQUIRE_NOTHROW(runAndGetNorm(3, 2, 4, 10, 10, 10,
                                 DACOUSTICS "data/acousticsGaussian2D.h", 1));
}

TEST_CASE("BC smoke - FreqIndep Tri2D", "[integration][freq_indep][tri2d]") {
  // Phase 1 placeholder: wall BC; norm assertion added in Phase 3 when upwindBC lands
  REQUIRE_NOTHROW(runAndGetNorm(3, 2, 4, 10, 10, 10,
                                 DACOUSTICS "data/acousticsGaussian2D.h", 1));
}

TEST_CASE("BC smoke - PerfRefl Tet3D", "[integration][perf_refl][tet3d]") {
  REQUIRE_NOTHROW(runAndGetNorm(6, 3, 2, 10, 10, 10,
                                 DACOUSTICS "data/acousticsGaussian3D.h", 1));
}

TEST_CASE("BC smoke - FreqIndep Tet3D", "[integration][freq_indep][tet3d]") {
  // Phase 1 placeholder: wall BC; norm assertion added in Phase 3
  REQUIRE_NOTHROW(runAndGetNorm(6, 3, 2, 10, 10, 10,
                                 DACOUSTICS "data/acousticsGaussian3D.h", 1));
}

// Custom main with MPI bracketing.
// Comm::Init / Comm::Finalize wrap MPI_Init / MPI_Finalize (from include/comm.hpp).
// The entire Catch2 session runs inside MPI scope because the solver requires it.
int main(int argc, char* argv[]) {
  Comm::Init(argc, argv);
  int result = Catch::Session().run(argc, argv);
  Comm::Finalize();
  return result;
}
