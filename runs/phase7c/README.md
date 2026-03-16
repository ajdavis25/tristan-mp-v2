# Phase 7C Notes

This directory contains the first Phase 7C demonstrator and its first
same-geometry control:

- input: `input_phase7c_demo`
- control input: `input_phase7c_control_uniformT`
- analysis: `validate_phase7c.py`

## Geometry and intent

Phase 7C is a two-zone shearing box with:

- inner core: `x/Lx in [0.30, 0.70]`
- transition bands: `[0.20, 0.30]` and `[0.70, 0.80]`
- outer wings: `[0.00, 0.20]` and `[0.80, 1.00]`

The inner region uses `backgr_T_inner = 1.0`, the outer region uses
`backgr_T_outer = 4.0`, and the MRI seed is localized to the inner region
with a cosine taper across the transition bands.

This is a propagation analogue test, not a strict "MRI off outside" setup.

Scientific question:

- does an inner-seeded favorable region produce a delayed or outward-spreading
  response in less favorable outer wings?
- if the response is not delayed, does the setup simply behave like a nearly
  global mode despite the localized seed?

## Build / run

Fresh build target:

```bash
source /apps/hpcx/2.18.1/hpcx-init.sh && hpcx_load
export HDF5_ROOT=/work/vmo703/hdf5-parallel
cmake -S /work/vmo703/tristan-mp-v2 -B /work/vmo703/tristan-mp-v2/build_3d_7c \
  -Duser=user_mri_7c -Ddim=3 -Dhdf5=ON -Dmpi08=ON
cmake --build /work/vmo703/tristan-mp-v2/build_3d_7c -j
```

SLURM launcher:

- `/work/vmo703/tristan-mp-v2/slurm/phase7c.slurm`
- `/work/vmo703/tristan-mp-v2/slurm/phase7c_control.slurm`

## Analysis

From this directory:

```bash
python3 validate_phase7c.py out_phase7c_demo \
  --input input_phase7c_demo \
  --save phase7c_demo.png \
  --csv phase7c_demo_summary.csv
```

Uniform-temperature control:

```bash
python3 validate_phase7c.py out_phase7c_control_uniformT \
  --input input_phase7c_control_uniformT \
  --save phase7c_control_uniformT.png \
  --csv phase7c_control_uniformT_summary.csv
```

The analysis uses only standard `flds.tot.*` dumps and reports:

- `E_perp(x,t)` heatmap
- `-BxBy(x,t)` heatmap
- region-averaged inner/outer `E_perp(t)`
- outward `Sx(t)` at the inner/outer interfaces
- peak delays and non-negative-lag inner-to-outer correlations

Correlation details:

- `zero_lag_corr_*` is the Pearson correlation over the full post-onset window
- `best_lag_*` / `corr_peak_*` are found from overlap-local Pearson correlations
- the lag search requires at least 50% overlap, so late-window artifacts do not dominate the result

## Current run set

Completed outputs in this directory:

- demonstrator: `out_phase7c_demo`
- uniform-temperature control: `out_phase7c_control_uniformT`

Saved analysis outputs:

- demonstrator:
  - `phase7c_demo.png`
  - `phase7c_demo_summary.csv`
- control:
  - `phase7c_control_uniformT.png`
  - `phase7c_control_uniformT_summary.csv`

## Pass / fail criteria used here

Implementation / verification pass:

- the run completes through `last = 7000`
- field dumps and params dumps are written through the final step
- the validator runs without errors and produces PNG + CSV outputs

Physics / validation pass for the intended propagation claim:

- the inner region peaks measurably earlier than the outer wings
- `best_lag_left` and/or `best_lag_right` is clearly positive
- zero-lag correlation is not effectively unity across the full post-onset window
- the `E_perp(x,t)` heatmap shows a visible inner-to-outer front or delayed outer brightening
- interface `Sx` has a directional interpretation beyond equal-and-opposite oscillations

## Current status

Implementation status: PASS

- both runs complete and are analyzable
- the analysis script now produces bounded, interpretable lag/correlation values

Physics status for the propagation claim: FAIL / NOT YET VALIDATED

- the localized seed does not produce a usable propagation experiment in the current configuration

### Demonstrator result

From `phase7c_demo_summary.csv`:

- `t_onset_inner = 20`
- `t_peak_inner = 540`
- `t_peak_outer_left = 540`
- `t_peak_outer_right = 540`
- `peak_delay_left = 0`
- `peak_delay_right = 0`
- `zero_lag_corr_left = 0.9994499501448587`
- `zero_lag_corr_right = 0.9998326261187758`
- `best_lag_left = 0`
- `best_lag_right = 0`

Interpretation:

- the hot-wing demonstrator responds essentially synchronously across inner and outer regions
- the best lag is zero on both sides
- the post-onset inner/outer correlation is effectively unity
- the mean outward interface flux is close to zero:
  - `mean_outward_flux_avg = -8.331900624374385e-09`

### Uniform-temperature control result

From `phase7c_control_uniformT_summary.csv`:

- `t_onset_inner = 20`
- `t_peak_inner = 580`
- `t_peak_outer_left = 600`
- `t_peak_outer_right = 600`
- `peak_delay_left = 20`
- `peak_delay_right = 20`
- `zero_lag_corr_left = 0.9996139329325208`
- `zero_lag_corr_right = 0.9994522918877216`
- `best_lag_left = 0`
- `best_lag_right = 0`

Interpretation:

- the control behaves the same way in the metric that matters here
- the one-output-step peak delay is not supported by the lag metric
- the response is still best described as near-synchronous and near-global
- the mean outward interface flux is again small:
  - `mean_outward_flux_avg = 1.7693098205902767e-08`

### Comparison conclusion

The key comparison is:

- demonstrator: `best_lag = 0`, zero-lag correlation ~ `0.9995-0.9998`
- control: `best_lag = 0`, zero-lag correlation ~ `0.9995-0.9996`

That means the current two-zone localized-seed setup does not isolate inward-to-outward propagation. The temperature contrast changes details only weakly; it does not create a clean delayed outer response.

## What this means for project flow

What is already established:

- Phase 7C is implemented
- the analysis pipeline is usable
- the current configuration gives a scientifically interpretable null result:
  the box responds nearly globally even with a localized seed

What is not established:

- a clean propagation front
- delayed triggering in the outer wings
- a configuration that can be presented as a positive Phase 7C propagation result

## PBH-centered interpretation

The notebook does not frame 7C as an optional side quest. It frames 7C as the
black-hole analogue test:

- phase 7C outcome: a surrogate for "MRI propagation away from a black hole"
- follow-up question: "does MRI propagate outward or die outside the black hole?"

Given that project framing, the current 7C result is not a final success. It is
the first null design. The current setup is useful because it rules out one
obvious interpretation:

- localizing the seed alone is not enough
- adding a hot unfavorable outer region is not enough
- in this periodic shearing-box configuration, the response is still nearly global

So yes: if the PBH-centered claim is the real target, propagation still needs to
be isolated before the project can move into a paper-production phase.

## What happened after 7C

Phase 7D was created as a stronger isolate:

- it kept the same geometry and localized seed
- it also localized the shearing-box drive to the inner zone
- the outer wings were therefore passive by construction

That still failed:

- see `/work/vmo703/tristan-mp-v2/runs/phase7d/phase7d_demo_summary.csv`
- `best_lag_left = 0`
- `best_lag_right = 0`
- `zero_lag_corr_left = 0.9999400970713795`
- `zero_lag_corr_right = 0.9998983034875565`

That means the next phase now has to be more structural than either 7C or 7D.

## Concrete next phase: proposed Phase 8

Phase 8 should be a redesigned PBH-propagation isolate, not a production sweep.

### 8A. source-plus-buffer control

Start from `/work/vmo703/tristan-mp-v2/user/user_mri_7d.F90` and build a new
`user_mri_8a.F90` with:

- inner-only seed
- inner-only drive
- wider `x` domain than the current `sizex = 8`
- uniform temperature everywhere for the first control
- explicit damping / absorbing buffer zones near the `x` edges

Implementation path in this codebase:

- add field damping in `userFieldBoundaryConditions()`
- add particle removal / damping in `userParticleBoundaryConditions()`
- reuse the absorbing-layer patterns already present in
  `/work/vmo703/tristan-mp-v2/user/user_psr_simple.F90`

### 8B. PBH-analogue hot-wing run

Only after 8A shows a real delayed outer response, restore the unfavorable outer
plasma:

- `backgr_T_inner = 1.0`
- `backgr_T_outer = 4.0`
- same source width, same damping layers, same analysis

This is the run that should answer the PBH-analogue question: can MRI-driven
activity launched near the source region survive and move outward into a less
favorable plasma?

### 8C. production-style parameter scan

Only after 8A and 8B pass should the project move into a paper-style scan:

- outer temperature
- outer-zone width
- source-zone width
- total box length in `x`

## Success criterion for the next phase

The next propagation design should count as successful only if it shows:

- positive propagation lag larger than a one-dump ambiguity
- delayed outer peaks relative to the inner source region
- zero-lag correlation below the current near-unity behavior
- a visible inner-to-outer front in `E_perp(x,t)`
- outward interface flux with a directional interpretation

Until then, Phase 7C should be read as a useful null result, not as proof that
MRI propagation away from the PBH analogue has been captured.
