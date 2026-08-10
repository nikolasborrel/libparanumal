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

//settings for acoustics solver
acousticsSettings_t::acousticsSettings_t(comm_t _comm):
  settings_t(_comm) {

  newSetting("DATA FILE",
             "data/acousticsGaussian2D.h",
             "Boundary and Initial conditions header");

  newSetting("DENSITY",
             "1.0",
             "Medium density [kg/m^3]");

  newSetting("SPEED OF SOUND",
             "1.0",
             "Speed of sound in the medium [m/s]");

  newSetting("FMAX",
             "1000.0",
             "Maximum frequency of interest [Hz]; the mesh is expected to resolve it");

  newSetting("SXYZ",
             "0",
             "Gaussian source width sigma [m]. If <= 0, derived from FMAX as "
             "sigma = 2c/(pi*FMAX). Exposed to initial-condition headers as p_sigma0");

  newSetting("SOURCE X",
             "0.0",
             "Gaussian source center x [m], exposed as p_srcX");

  newSetting("SOURCE Y",
             "0.0",
             "Gaussian source center y [m], exposed as p_srcY");

  newSetting("SOURCE Z",
             "0.0",
             "Gaussian source center z [m], exposed as p_srcZ");

  newSetting("FREQINDEP IMPEDANCE",
             "415.0",
             "Acoustic impedance for frequency-independent BC [Pa·s/m] (default: air at 20°C, rho*c)");

  newSetting("SURFACE FLUX",
             "UPWIND",
             "Numerical flux used in the surface integral",
             {"UPWIND", "CENTRAL"});

  newSetting("TIME INTEGRATOR",
             "DOPRI5",
             "Time integration method. EIRK4 is an additive Runge-Kutta scheme "
             "that treats the locally-reacting accumulators implicitly, so a "
             "stiff surface admittance does not restrict the time step",
             {"AB3", "DOPRI5", "LSERK4", "EIRK4"});

  newSetting("CFL NUMBER",
             "1.0",
             "Multiplier for timestep stability bound");

  // dt = cfl*hmin/(c*(N+1)^2) uses the physical speed of sound, so these are
  // physical seconds
  newSetting("START TIME",
             "0",
             "Start time for time integration [s]");

  newSetting("FINAL TIME",
             "10",
             "End time for time integration [s]");

  newSetting("OUTPUT INTERVAL",
             ".1",
             "Time between output snapshots [s]");

  newSetting("OUTPUT TO FILE",
             "FALSE",
             "Flag for writing fields to VTU files. Independent of OUTPUT FORMAT: "
             "a run may write VTU, an HDF5 format, both, or neither",
             {"TRUE", "FALSE"});

  newSetting("OUTPUT FORMAT",
             "NONE",
             "Wave-field snapshot format, written on the OUTPUT INTERVAL cadence: "
             "H5COMPACT (single /pressures dataset) or XDMF (per-step datasets "
             "plus a ParaView sidecar). Requires a build with HDF5=1",
             {"NONE", "H5COMPACT", "XDMF"});

  newSetting("OUTPUT FILE NAME",
             "acoustics");

  newSetting("LR VECTORFIT FILE",
             "",
             "Path to vectorfit data file for LR BCs (optional; required when BOX BOUNDARY FLAG or mesh has BC type 3)");

  newSetting("RECEIVER FILE",
             "",
             "Path to receiver locations file (optional)");

  newSetting("OUTPUT DIRECTORY",
             ".",
             "Directory for field/receiver output files (created if absent)");

  newSetting("SIMULATION ID",
             "acoustics",
             "Stem name for run-scoped output files (<id>.log, <id>_receivers.h5)");
}

void acousticsSettings_t::report() {

  if (comm.rank()==0) {
    std::cout << "Acoustics Settings:\n\n";
    reportSetting("DATA FILE");
    reportSetting("DENSITY");
    reportSetting("SPEED OF SOUND");
    reportSetting("FMAX");
    reportSetting("SURFACE FLUX");
    reportSetting("TIME INTEGRATOR");
    reportSetting("START TIME");
    reportSetting("FINAL TIME");
    reportSetting("OUTPUT INTERVAL");
    reportSetting("OUTPUT TO FILE");
    reportSetting("OUTPUT FORMAT");
    reportSetting("OUTPUT FILE NAME");
  }
}

void acousticsSettings_t::parseFromFile(platformSettings_t& platformSettings,
                                  meshSettings_t& meshSettings,
                                  const std::string filename) {
  //read all settings from file
  settings_t s(comm);
  s.readSettingsFromFile(filename);

  for(auto it = s.settings.begin(); it != s.settings.end(); ++it) {
    setting_t& set = it->second;
    const std::string name = set.getName();
    const std::string val = set.getVal<std::string>();
    if (platformSettings.hasSetting(name))
      platformSettings.changeSetting(name, val);
    else if (meshSettings.hasSetting(name))
      meshSettings.changeSetting(name, val);
    else if (hasSetting(name)) //self
      changeSetting(name, val);
    else  {
      LIBP_FORCE_ABORT("Unknown setting: [" << name << "] requested");
    }
  }
}
