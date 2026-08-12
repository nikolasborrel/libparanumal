# Acoustics solver — output and data generation

This covers what the acoustics solver can write and how to configure it: VTU
snapshots, HDF5/XDMF wave fields, receiver impulse responses, and field sampling
onto configurable point sets for machine-learning datasets.

## Building

HDF5 output is opt-in, so the default build keeps the original dependency set:

```
make                 # no HDF5; VTU output only
make HDF5=1          # enables the .h5 / .xdmf output paths
```

`HDF5=1` needs `libhdf5-dev` and the header-only
[HighFive](https://github.com/BlueBrain/HighFive) library. Fetch HighFive once:

```
make fetch-highfive          # clones a pinned tag into third_party/
```

or point `HIGHFIVE_INC` at an existing checkout. Fetching is never a build
prerequisite, so an ordinary build does not touch the network. Toggling `HDF5`
recompiles rather than silently relinking stale objects.

Selecting an HDF5 output in a build without it aborts with that instruction,
except receiver impulse responses, which degrade to a notice.

## Output files

Everything is written under `OUTPUT DIRECTORY` (created if absent) and stemmed
on `SIMULATION ID`:

| File | Written when |
|---|---|
| `<name>_<rank>_<frame>.vtu` | `OUTPUT TO FILE = TRUE` |
| `<id>.h5`, `<id>.xdmf` | `OUTPUT FORMAT = H5COMPACT` or `XDMF` |
| `<id>_receivers.h5` | `RECEIVER FILE` is set |
| `<id>_samples.h5` | `SAMPLE SETS FILE` is set |
| `<id>.log` | always; mirrors the key diagnostics printed to stdout |

These are independent. A run may write all of them, or none.

## Source configuration

Initial-condition headers that use the `p_sigma0` and `p_srcX/Y/Z` defines —
`data/acousticsGaussianFmax3D.h` is the reference one — take their source from
settings instead of hard-coded constants:

| Setting | Meaning |
|---|---|
| `FMAX` | highest frequency of interest [Hz] |
| `SXYZ` | Gaussian width sigma [m]; if `<= 0`, derived as `2c/(pi*FMAX)` |
| `SOURCE X`, `SOURCE Y`, `SOURCE Z` | pulse centre [m] |

The centre is three scalars, not one `"x y z"` string, because the settings
reader strips whitespace from values.

`FMAX` also sets the sample grid spacings and the `ppp` cadence below, so it is
usually the one number to change when moving to a different frequency range.

## Wave-field snapshots

`OUTPUT FORMAT` writes the pressure field on the `OUTPUT INTERVAL` cadence:

- `H5COMPACT` — one chunked `/pressures [nFrames x N]` dataset plus `/mesh`,
  with the output times as an attribute.
- `XDMF` — one dataset per snapshot plus a temporal-collection sidecar. Open the
  `.xdmf` in ParaView or VisIt.

The sidecar is written after the time loop from the frames that actually
exist, so it never references a missing snapshot. `OUTPUT INTERVAL` is clamped
to the solver time step, since the solver cannot emit more than one snapshot per
step.

Example: `setups/setupWavefieldTet3D.rc`.

## Receivers

`RECEIVER FILE` points at a text file holding a count followed by one `x y z`
per line. Receivers are sampled on the `OUTPUT INTERVAL` cadence and written to
`<id>_receivers.h5`:

```
/impulse_responses  [NReceivers x nSamples]   attrs: sample_rate_hz, n_samples, n_receivers
/positions          [NReceivers x 3]
```

The record is gathered before rank 0 writes, and receivers keep their file
order, so the layout does not depend on how many ranks the run used.

## Field sampling

`SAMPLE SETS FILE` names a file declaring any number of sample sets. Each set
records the pressure on its own point grid, with its own temporal cadence, into
its own group of `<id>_samples.h5`. One set per line, `#` starts a comment:

```
name=trunk ppw=4 jitter=0.5 seed=0 cadence=ppp:10
name=branch ppw=2 outside=keep cadence=initial
```

| Key | Meaning |
|---|---|
| `name` | group name in the output file |
| `ppw=<n>` | spacing as points per wavelength at `FMAX`, i.e. `dx = c/(FMAX*n)` |
| `dx=<m>` | spacing in metres, instead of `ppw` |
| `jitter=<f>` | offset each point by `U(-f*dx, +f*dx)` per axis |
| `seed=<i>` | RNG seed for the jitter |
| `outside=drop` | drop points outside the domain (default) — a scattered set |
| `outside=keep` | keep them as zeros, preserving the rectilinear shape |
| `cadence=initial` | sample the initial condition only |
| `cadence=interval:<s>` | sample every `<s>` seconds |
| `cadence=ppp:<n>` | sample `n` times per period at `FMAX` (`n >= 2`) |

Give exactly one of `ppw` or `dx`. A set with a step cadence requires
`TIME INTEGRATOR LSERK4`, since sampling at exact steps needs fixed step
control. Sets are independent of `OUTPUT INTERVAL` and of the wave-field
snapshots.

The grid always spans the mesh bounding box.

### Output layout

Per set `<name>`:

```
/<name>/points   [npts x 3]        attrs: grid_shape, dx, jitter, seed, keep_outside
/<name>/values   [nframes x npts]
/<name>/times    [nframes]
```

`grid_shape` is in C order `(nz, ny, nx)`, so for an `outside=keep` set
`values[f].reshape(grid_shape)` gives the field as a 3-D tensor directly. Points
are stored with x varying fastest.

## Generating a dataset

`setups/setupRoomDataGen.rc` with `setups/sampleSetsDeepONet.txt` is a worked
example: a 3 m impedance-walled room at `FMAX = 500 Hz`, a jittered evaluation
grid over time, and a coarser rectilinear grid of the initial condition.

```
$ make HDF5=1
$ ./acousticsMain setups/setupRoomDataGen.rc
Sample set 'trunk': 18 x 18 x 18 grid, dx=0.1715 m, 5374 point(s) kept
  requested 0.0002 s between samples, solver dt=4.123e-05 s -> using 0.0002062 s
Sample set 'branch': 9 x 9 x 9 grid, dx=0.343 m, 729 point(s) kept
  wrote out/room_train_000_samples.h5:/trunk (5374 points x 98 frames)
  wrote out/room_train_000_samples.h5:/branch (729 points x 1 frames)
```

`trunk` kept 5374 of 5832 points because jitter pushes boundary points outside
the room and `outside=drop` removes them. `branch` kept all 729.

Reading it back, for an operator-learning model that maps the initial condition
to the field at arbitrary points in space and time:

```python
import h5py, numpy as np

with h5py.File("out/room_train_000_samples.h5") as f:
    u = f["branch/values"][0].reshape(f["branch/points"].attrs["grid_shape"])
    y = f["trunk/points"][:]        # (5374, 3)
    t = f["trunk/times"][:]         # (98,)
    s = f["trunk/values"][:]        # (98, 5374)

# trunk coordinates as space-time pairs (x, y, z, t)
yt = np.concatenate([np.repeat(y, len(t), 0), np.tile(t, len(y))[:, None]], 1)
target = s.T.reshape(-1)
```

The group names above are only HDF5 group names; nothing in the solver is tied
to a model architecture. Rename them, add a third set, or sample a single plane
by choosing `dx` and a thin domain.

### Sweeping over runs

Vary the source position and the jitter seed per run, and give each its own
`SIMULATION ID`:

```bash
for i in $(seq 0 199); do
  id=room_train_$(printf %03d $i)
  sed "s/seed=0/seed=$i/" setups/sampleSetsDeepONet.txt > /tmp/ss_$i.txt
  x=$(python3 -c "import random;random.seed($i);print(round(random.uniform(-1.2,1.2),3))")
  sed -e "s|setups/sampleSetsDeepONet.txt|/tmp/ss_$i.txt|" \
      -e "s/^room_train_000$/$id/" -e "s/^0.6$/$x/" \
      setups/setupRoomDataGen.rc > /tmp/run_$i.rc
  ./acousticsMain /tmp/run_$i.rc
done
```

Different source positions are what give the model different input functions;
different seeds move the evaluation points between runs.

### Sizing

Cost and file size scale as `points x frames`. The example above is 2.2 MB for
0.02 s. Doubling `FMAX` shrinks `dx` by two in each axis and doubles the `ppp`
rate, so roughly 16x; a full half-second impulse response is another 25x. The
grid size is printed before the time loop — check it before launching a sweep.

### Cadence rounding

Sampling lands on whole solver steps, so a requested interval is rounded to the
nearest multiple of `dt` and the achieved value is reported. If rounding ever
pushes the rate below two samples per period at `FMAX`, the run warns that the
result aliases.
