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
             "Maximum frequency of interest [Hz]. A source/mesh property: the "
             "mesh is resolved to FMAX and the data sampler derives its grid "
             "spacing (dx = c/(FMAX*ppw)) and temporal cadence from it.");

  newSetting("SXYZ",
             "0",
             "Gaussian source width sigma [m]. If <= 0, computed from FMAX as "
             "sigma = 2c/(pi*FMAX) (DTU convention). Used by FMAX-parameterized "
             "initial-condition headers via the kernel define p_sigma0.");

  newSetting("SOURCE POSITION",
             "0.0 0.0 0.0",
             "Gaussian source center 'x y z' [m]. Exposed to FMAX-parameterized "
             "initial-condition headers via the defines p_srcX/p_srcY/p_srcZ.");

  newSetting("FREQINDEP IMPEDANCE",
             "415.0",
             "Acoustic impedance for frequency-independent BC [Pa·s/m] (default: air at 20°C, rho*c)");

  newSetting("SURFACE FLUX",
             "UPWIND",
             "Numerical flux used in the surface integral",
             {"UPWIND", "CENTRAL"});

  newSetting("TIME INTEGRATOR",
             "DOPRI5",
             "Time integration method",
             {"AB3", "DOPRI5", "LSERK4"});

  newSetting("CFL NUMBER",
             "1.0",
             "Multiplier for timestep stability bound");

  // Time settings are used directly by the time loop, exactly as in the other
  // solvers. Because dt = cfl*hmin/(c*(N+1)^2) uses the physical speed of sound,
  // these values are in physical seconds.
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
             "Flag for writing wave-field snapshots to VTU files (visualization), "
             "using the same VTU path as the other solvers. Independent of "
             "OUTPUT FORMAT (HDF5/XDMF): a run may write VTU, an HDF5 format, "
             "both, or neither. (Name kept for cross-solver consistency.)",
             {"TRUE", "FALSE"});

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
             "Directory for field/receiver output files");

  newSetting("SIMULATION ID",
             "acoustics",
             "Stem name for .h5 / .xdmf output files");

  // ---- ML data-generation output (orthogonal to visualization) ----------
  // These typically live in a separate file referenced by DATA CONFIG FILE,
  // so train/val/test splits share the main mesh/physics/viz setup and swap
  // only the data config (e.g. different DATA SEED / SPATIAL PPW DATA).
  newSetting("DATA CONFIG FILE",
             "",
             "Path to a separate .rc with data-generation settings. Setting this "
             "(in the main config) enables ML data-generation; leaving it empty "
             "disables it. VTU visualization is controlled independently by "
             "OUTPUT TO FILE in the main config.");

  newSetting("OUTPUT FORMAT",
             "NONE",
             "Wave-field snapshot format: H5COMPACT (single-file /pressures) or "
             "XDMF (per-step datasets + ParaView sidecar). Lives in the data config "
             "and is MANDATORY for data-generation runs. Independent of VTU output "
             "(OUTPUT TO FILE). NONE = no wave-field output.",
             {"NONE", "H5COMPACT", "XDMF"});

  newSetting("SPATIAL PPW DATA",
             "4",
             "Points-per-wavelength for the data sample grid at FMAX (sets dx = c/(FMAX*ppw))");

  newSetting("DATA TEMPORAL PPW",
             "20.0",
             "Temporal samples per period at FMAX (>=2 for Nyquist); sets data sampling dt");

  newSetting("DATA JITTER",
             "0.5",
             "Sample-point coordinate jitter as a fraction of dx; offset ~ Uniform(-f*dx, +f*dx) per axis");

  newSetting("DATA SOURCE PPW",
             "2",
             "Points-per-wavelength for the rectilinear SOURCE grid at FMAX (the "
             "initial condition sampled for ML input; DTU MESH_RECTILINEAR_PPW). "
             "<=0 disables source-grid output.");

  newSetting("DATA SOURCE JITTER",
             "0",
             "Source-grid coordinate jitter as a fraction of dx (default 0 = uniform rectilinear grid)");

  newSetting("DATA SEED",
             "0",
             "RNG seed for the sample-grid coordinate jitter (vary across dataset splits)");

  newSetting("DATA SIMULATION ID",
             "",
             "Stem name for the _data.h5 file (defaults to SIMULATION ID when empty)");
}

void acousticsSettings_t::report() {

  if (comm.rank()==0) {
    std::cout << "Acoustics Settings:\n\n";
    reportSetting("DATA FILE");
    reportSetting("DENSITY");
    reportSetting("SPEED OF SOUND");
    reportSetting("SURFACE FLUX");
    reportSetting("TIME INTEGRATOR");
    reportSetting("START TIME");
    reportSetting("FINAL TIME");
    reportSetting("OUTPUT INTERVAL");
    reportSetting("OUTPUT TO FILE");
    reportSetting("OUTPUT FILE NAME");
  }
}

void acousticsSettings_t::parseFromFile(platformSettings_t& platformSettings,
                                  meshSettings_t& meshSettings,
                                  const std::string filename) {

  // Apply every setting from one .rc file, routing each key to the settings
  // object that owns it (platform / mesh / self).
  auto applyFile = [&](const std::string& fname) {
    settings_t s(comm);
    s.readSettingsFromFile(fname);

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
  };

  // Main setup file first ...
  applyFile(filename);

  // ... then merge an optional modular data-generation config on top, so a
  // dataset split overrides only data-gen keys (DATA SEED, SPATIAL PPW DATA, …)
  // while sharing the main mesh/physics/visualization setup.
  std::string dataConfig;
  if (hasSetting("DATA CONFIG FILE")) {
    getSetting("DATA CONFIG FILE", dataConfig);
    if (!dataConfig.empty())
      applyFile(dataConfig);
  }
}
