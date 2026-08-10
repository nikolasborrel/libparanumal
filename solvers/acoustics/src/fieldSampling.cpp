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
#include "pointSampling.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <random>
#include <sstream>

namespace {

// key=value token, e.g. "ppw=4" or "cadence=interval:0.05"
bool splitToken(const std::string& tok, std::string& key, std::string& val) {
  const size_t eq = tok.find('=');
  if (eq == std::string::npos) return false;
  key = tok.substr(0, eq);
  val = tok.substr(eq + 1);
  return true;
}

} // namespace

// One set per line, "#" starts a comment:
//   name=volume ppw=4 jitter=0.5 seed=0 cadence=interval:0.05
//   name=source ppw=2 outside=keep cadence=initial
// Spacing is ppw=<points per wavelength at FMAX> or dx=<metres>.
void acoustics_t::ParseSampleSets(const std::string& path) {
  std::ifstream file(path);
  LIBP_ABORT("Cannot open SAMPLE SETS FILE: " << path, !file.is_open());

  std::string line;
  int lineNo = 0;
  while (std::getline(file, line)) {
    ++lineNo;
    const size_t hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);

    std::istringstream iss(line);
    std::vector<std::string> tokens{std::istream_iterator<std::string>(iss),
                                    std::istream_iterator<std::string>()};
    if (tokens.empty()) continue;

    sampleSet_t set;
    bool haveCadence = false;
    for (const std::string& tok : tokens) {
      std::string key, val;
      LIBP_ABORT("Sample set " << path << ":" << lineNo
                 << ": expected key=value, got '" << tok << "'",
                 !splitToken(tok, key, val));

      if      (key == "name")    set.name   = val;
      else if (key == "ppw")     set.ppw    = std::stod(val);
      else if (key == "dx")      set.dx     = std::stod(val);
      else if (key == "jitter")  set.jitter = std::stod(val);
      else if (key == "seed")    set.seed   = std::stoi(val);
      else if (key == "outside") {
        LIBP_ABORT("Sample set " << path << ":" << lineNo
                   << ": outside must be drop or keep", val != "drop" && val != "keep");
        set.keepOutside = (val == "keep");
      }
      else if (key == "cadence") {
        haveCadence = true;
        if (val == "initial") {
          set.cadence = sampleSet_t::Cadence::INITIAL;
        } else {
          const size_t colon = val.find(':');
          LIBP_ABORT("Sample set " << path << ":" << lineNo
                     << ": cadence must be initial, interval:<s> or ppp:<n>",
                     colon == std::string::npos);
          const std::string mode = val.substr(0, colon);
          const dfloat      num  = std::stod(val.substr(colon + 1));
          set.cadence = sampleSet_t::Cadence::INTERVAL;
          if (mode == "interval") {
            set.interval = num;
          } else if (mode == "ppp") {
            LIBP_ABORT("Sample set " << path << ":" << lineNo
                       << ": ppp must be >= 2 to resolve FMAX", num < 2.0);
            set.interval = 1.0/(fmax*num);
          } else {
            LIBP_FORCE_ABORT("Sample set " << path << ":" << lineNo
                             << ": unknown cadence '" << mode << "'");
          }
        }
      }
      else LIBP_FORCE_ABORT("Sample set " << path << ":" << lineNo
                            << ": unknown key '" << key << "'");
    }

    LIBP_ABORT("Sample set " << path << ":" << lineNo << ": missing name", set.name.empty());
    LIBP_ABORT("Sample set '" << set.name << "': give exactly one of ppw= or dx=",
               (set.ppw > 0.0) == (set.dx > 0.0));
    LIBP_ABORT("Sample set '" << set.name << "': missing cadence", !haveCadence);

    if (set.ppw > 0.0) {
      LIBP_ABORT("Sample set '" << set.name << "': ppw needs FMAX > 0", fmax <= 0.0);
      set.dx = c/(fmax*set.ppw);
    }
    sampleSets.push_back(std::move(set));
  }
}

void acoustics_t::SetupFieldSampling() {
  std::string path;
  if (settings.hasSetting("SAMPLE SETS FILE"))
    settings.getSetting("SAMPLE SETS FILE", path);
  if (path.empty()) return;

  ParseSampleSets(path);
  if (sampleSets.empty()) return;

  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  // global bounding box, so every rank generates the identical grid
  dfloat lo[3] = { 1e300,  1e300,  1e300};
  dfloat hi[3] = {-1e300, -1e300, -1e300};
  const dlong Nnode = mesh.Nelements*mesh.Np;
  for (dlong n = 0; n < Nnode; ++n) {
    lo[0] = std::min(lo[0], mesh.x[n]); hi[0] = std::max(hi[0], mesh.x[n]);
    lo[1] = std::min(lo[1], mesh.y[n]); hi[1] = std::max(hi[1], mesh.y[n]);
    if (mesh.dim == 3) { lo[2] = std::min(lo[2], mesh.z[n]); hi[2] = std::max(hi[2], mesh.z[n]); }
  }
  if (mesh.dim != 3) { lo[2] = 0.0; hi[2] = 0.0; }
  for (int d = 0; d < 3; ++d) {
    mesh.comm.Allreduce(lo[d], Comm::Min);
    mesh.comm.Allreduce(hi[d], Comm::Max);
  }

  for (sampleSet_t& set : sampleSets) {
    const int nx = std::max(1, (int)std::floor((hi[0]-lo[0])/set.dx) + 1);
    const int ny = std::max(1, (int)std::floor((hi[1]-lo[1])/set.dx) + 1);
    const int nz = (mesh.dim == 3)
                 ? std::max(1, (int)std::floor((hi[2]-lo[2])/set.dx) + 1) : 1;
    set.shape[0] = nx; set.shape[1] = ny; set.shape[2] = nz;

    std::mt19937 rng((unsigned)set.seed);
    std::uniform_real_distribution<dfloat> jit(-set.jitter*set.dx, set.jitter*set.dx);

    memory<dfloat> grid((dlong)nx*ny*nz*3);
    dlong g = 0;
    for (int k = 0; k < nz; ++k)
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
          dfloat jx = 0.0, jy = 0.0, jz = 0.0;
          if (set.jitter > 0.0) {
            jx = jit(rng); jy = jit(rng);
            if (mesh.dim == 3) jz = jit(rng);
          }
          grid[g*3+0] = lo[0] + i*set.dx + jx;
          grid[g*3+1] = lo[1] + j*set.dx + jy;
          grid[g*3+2] = (mesh.dim == 3) ? (lo[2] + k*set.dx + jz) : 0.0;
          ++g;
        }

    memory<dlong>  elem;
    memory<dfloat> weights;
    pointSampling::locate(mesh, grid, g, elem, weights);

    // keep every point (zeros outside) or only those located on this rank
    std::vector<dlong> keep;
    keep.reserve(g);
    for (dlong p = 0; p < g; ++p)
      if (set.keepOutside || elem[p] >= 0) keep.push_back(p);

    set.Npoints = (dlong)keep.size();
    Logf("Sample set '%s': %d x %d x %d grid, dx=%.4g m, %lld point(s) kept\n",
         set.name.c_str(), nx, ny, nz, set.dx, (long long)set.Npoints);
    if (set.Npoints == 0) continue;

    if (set.cadence == sampleSet_t::Cadence::INITIAL) {
      set.stride  = 0;
      set.Nframes = 1;
    } else {
      set.stride = std::max(1, (int)std::lround(set.interval/dt));
      const dfloat achieved = set.stride*dt;
      set.Nframes = (dlong)((finalTime - startTime)/achieved) + 2;
      if (std::fabs(achieved - set.interval) > 1e-9*set.interval)
        Logf("  requested %.4g s between samples, solver dt=%.4g s -> using %.4g s\n",
             set.interval, dt, achieved);
      // rounding the stride can push the rate below two samples per period
      if (fmax > 0.0 && achieved > 1.0/(2.0*fmax))
        Logf("  WARNING: %.4g s between samples aliases FMAX=%.4g Hz "
             "(needs <= %.4g s)\n", achieved, fmax, 1.0/(2.0*fmax));
    }

    const int Np = mesh.Np;
    set.xyz.malloc(set.Npoints*3);
    set.elements.malloc(set.Npoints);
    set.identity.malloc(set.Npoints);
    memory<dfloat> ip(set.Npoints*Np, 0.0);
    for (dlong r = 0; r < set.Npoints; ++r) {
      const dlong p = keep[r];
      set.xyz[r*3+0] = grid[p*3+0];
      set.xyz[r*3+1] = grid[p*3+1];
      set.xyz[r*3+2] = grid[p*3+2];
      // points outside the domain keep element 0 with zero weights, so the
      // kernel reads a valid offset and contributes nothing
      set.elements[r] = (elem[p] >= 0) ? elem[p] : 0;
      set.identity[r] = r;
      if (elem[p] >= 0)
        for (int n = 0; n < Np; ++n) ip[r*Np+n] = weights[p*Np+n];
    }

    set.values.malloc(set.Npoints*set.Nframes, 0.0);
    set.o_ip        = platform.malloc<dfloat>(ip);
    set.o_elements  = platform.malloc<dlong>(set.elements);
    set.o_identity  = platform.malloc<dlong>(set.identity);
    set.o_values    = platform.malloc<dfloat>(set.values);
    set.times.reserve(set.Nframes);

    // the receiver kernel already interpolates field 0 at arbitrary points; it
    // pins the frame count as a define, so each set gets its own build
    properties_t info = mesh.props;
    info["defines/p_blockSize"]    = 256;
    info["defines/p_Np"]           = Np;
    info["defines/p_Nfields"]      = Nfields;
    info["defines/p_NRecvSamples"] = (int)set.Nframes;
    set.kernel = platform.buildKernel(DACOUSTICS "okl/acousticsReceiverKernel.okl",
                                      "acousticsReceiverInterpolation", info);
  }

  // sampling on a step cadence needs the fixed-step loop in Run()
  LIBP_ABORT("Sample sets with a step cadence require TIME INTEGRATOR LSERK4 or EIRK4",
             SampleSetsNeedStepControl() &&
             !settings.compareSetting("TIME INTEGRATOR", "LSERK4") && !useEIRK4);
}

void acoustics_t::SampleFields(int tstep, dfloat time) {
  for (sampleSet_t& set : sampleSets) {
    if (set.Npoints == 0 || set.frame >= set.Nframes) continue;

    const bool due = (set.cadence == sampleSet_t::Cadence::INITIAL)
                   ? (tstep == 0)
                   : (tstep % set.stride == 0);
    if (!due) continue;

    set.kernel(set.Npoints, (dlong)set.frame,
               set.o_values, set.o_elements, set.o_identity, set.o_ip, o_q);
    set.times.push_back(time);
    set.frame++;
  }
}

bool acoustics_t::SampleSetsNeedStepControl() const {
  for (const sampleSet_t& set : sampleSets)
    if (set.Npoints > 0 && set.cadence != sampleSet_t::Cadence::INITIAL) return true;
  return false;
}
