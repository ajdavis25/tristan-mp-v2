# Phase 8A Notes

Phase 8A is the first PBH-propagation redesign after both 7C and 7D produced
near-global synchronous response.

This run is a source-plus-buffer control:

- inner-only seed
- inner-only drive
- uniform-temperature plasma
- passive outer wings
- explicit x-edge field/particle damping
- absorbing-x build (`-Dabsorb=ON`)

## Files

- input: `input_phase8a_control`
- analysis: `validate_phase8a.py`
- launcher: `/work/vmo703/tristan-mp-v2/slurm/phase8a.slurm`
- userfile: `/work/vmo703/tristan-mp-v2/user/user_mri_8a.F90`

## Build

```bash
source /apps/hpcx/2.18.1/hpcx-init.sh && hpcx_load
export HDF5_ROOT=/work/vmo703/hdf5-parallel
cmake -S /work/vmo703/tristan-mp-v2 -B /work/vmo703/tristan-mp-v2/build_3d_8a \
  -Duser=user_mri_8a -Ddim=3 -Dhdf5=ON -Dmpi08=ON -Dabsorb=ON
cmake --build /work/vmo703/tristan-mp-v2/build_3d_8a -j
```

## Analysis

```bash
cd /work/vmo703/tristan-mp-v2/runs/phase8
python3 validate_phase8a.py out_phase8a_control \
  --input input_phase8a_control \
  --save phase8a_control.png \
  --csv phase8a_control_summary.csv
```

The validator excludes the x-edge buffers from the passive-wing averages and lag
metrics. It keeps the source-interface Poynting-flux measurement at the source
edges.

## Current status

The completed `phase8a_control` run is an operational success but a physics
failure for the PBH-propagation isolate.

Observed summary from `phase8a_control_summary.csv`:

- `peak_delay_left = 40`, `peak_delay_right = 0`
- `best_lag_left = 0`, `best_lag_right = 0`
- `zero_lag_corr_left = 0.9999974`
- `zero_lag_corr_right = 0.9999963`
- `mean_outward_flux_avg = 7.41e-10`

Interpretation:

- the left passive wing peaked one dump later than the source, but the lag
  metric still selects zero lag
- both passive wings remain nearly perfectly synchronous with the source
- the source-interface flux still averages to nearly zero, so the run does not
  isolate directional outward propagation

This means 8A did not yet resolve the PBH-centered question. The immediate
follow-up is a longer x-domain with the same control physics.

## Phase 8A-long

`phase8a_long` keeps the same source-plus-buffer control design and the same
validator, but extends the x-domain so the passive wings are physically farther
from the source before the absorbing buffers begin.

- input: `input_phase8a_long_control`
- launcher: `/work/vmo703/tristan-mp-v2/slurm/phase8a_long.slurm`
- analysis:

```bash
cd /work/vmo703/tristan-mp-v2/runs/phase8
python3 validate_phase8a.py out_phase8a_long_control \
  --input input_phase8a_long_control \
  --save phase8a_long_control.png \
  --csv phase8a_long_control_summary.csv
```

Long-box control parameters:

- `sizex = 24`, `mx0 = 384`
- `sizey = 4`, `my0 = 64`
- `sizez = 2`, `mz0 = 32`
- `last = 15000`
- same source fractions, buffer fractions, damping, and uniform-temperature
  control physics as 8A

## Success criterion

8A is a success only if at least one passive wing shows:

- `best_lag >= 40`
- `peak_delay >= 40`
- `corr_peak >= 0.8`
- `zero_lag_corr <= 0.995`

and the plots agree with that interpretation:

- `E_perp(x,t)` shows source-to-wing progression
- the passive-wing signal appears before the edge buffers dominate
- source-interface outward flux is more than equal-and-opposite oscillation

If this still collapses to `best_lag = 0` and near-unity zero-lag correlation,
the immediate follow-up is to revisit the experiment design itself rather than
continuing to increase box size without changing the source/outer-medium
coupling.
