# Phase 7D Notes

Phase 7D preserves the Phase 7C null result and introduces a new design
whose purpose is narrower and cleaner:

- Phase 7C tested whether a localized seed plus two-zone temperature split
  was enough to isolate inward-to-outward propagation.
- It was not. Both the hot-wing demonstrator and the uniform-temperature
  control responded nearly synchronously.
- Phase 7D therefore changes the experiment definition rather than
  retroactively changing Phase 7C.

## Design change

Phase 7D keeps:

- the same `8 x 4 x 2` box
- the same `128 x 64 x 32` mesh
- the same localized seed geometry
- the same uniform vertical field `B0_ext = 0.06`

Phase 7D changes one essential thing:

- the shearing-box drive is localized to the inner zone with the same cosine
  taper used for the seed

That means:

- inner band: actively driven MRI source region
- outer wings: passive magnetized medium

This is a cleaner propagation isolate than Phase 7C because the outer region
is no longer locally MRI-active by construction.

## Files

- input: `input_phase7d_demo`
- analysis: `validate_phase7d.py`
- launcher: `/work/vmo703/tristan-mp-v2/slurm/phase7d.slurm`
- userfile: `/work/vmo703/tristan-mp-v2/user/user_mri_7d.F90`

## Build

```bash
source /apps/hpcx/2.18.1/hpcx-init.sh && hpcx_load
export HDF5_ROOT=/work/vmo703/hdf5-parallel
cmake -S /work/vmo703/tristan-mp-v2 -B /work/vmo703/tristan-mp-v2/build_3d_7d \
  -Duser=user_mri_7d -Ddim=3 -Dhdf5=ON -Dmpi08=ON
cmake --build /work/vmo703/tristan-mp-v2/build_3d_7d -j
```

## Analysis

The Phase 7D validator is the same metric set as Phase 7C, but the result
is interpreted differently:

- a positive result should now show delayed outer response because the outer
  zone is passive
- if Phase 7D still gives `best_lag = 0` and near-unity zero-lag correlation,
  the conclusion is stronger: this setup is exciting a nearly global response
  even without outer-region drive

Run analysis:

```bash
python3 validate_phase7d.py out_phase7d_demo \
  --input input_phase7d_demo \
  --save phase7d_demo.png \
  --csv phase7d_demo_summary.csv
```

## Current result

The first full 7D demonstrator completed and analyzed cleanly, but it did not
produce a propagation isolate.

From `phase7d_demo_summary.csv`:

- `t_onset_inner = 20`
- `t_peak_inner = 580`
- `t_peak_outer_left = 580`
- `t_peak_outer_right = 580`
- `peak_delay_left = 0`
- `peak_delay_right = 0`
- `zero_lag_corr_left = 0.9999400970713795`
- `zero_lag_corr_right = 0.9998983034875565`
- `best_lag_left = 0`
- `best_lag_right = 0`
- `mean_outward_flux_avg = 1.5236480937027477e-08`

Interpretation:

- localizing the drive did not create a delayed outer response
- inner and outer regions still evolve almost synchronously
- the average outward-flux signal is tiny compared with the oscillatory left/right
  interface behavior

So 7D is a stronger null result than 7C, not a positive propagation result.

## PBH-centered interpretation

This matters because the notebook frames the propagation branch in black-hole
language:

- Phase 7C is the surrogate for MRI propagation away from a black hole
- the stated question is whether MRI propagates outward or dies outside the
  black-hole source region

Under that success criterion, 7D does not close the project. It shows that even
after making the outer wings passive, this shearing-box surrogate still behaves
too globally to support the PBH-propagation claim.

## Concrete next phase: proposed Phase 8

The next phase should be a propagation redesign, not a production run.

### 8A. redesign the box to suppress global recirculation

Start from `user_mri_7d.F90` and move to a new `user_mri_8a.F90` with:

- the same inner-only source concept
- a longer `x` domain
- passive outer wings
- damping / absorbing buffer layers near the `x` edges

Implementation should stay inside existing user hooks:

- `userFieldBoundaryConditions()`
- `userParticleBoundaryConditions()`

There is already absorbing-boundary logic elsewhere in the repository,
especially `user/user_psr_simple.F90`, so this is a codebase-native redesign.

### 8B. rerun the source-only control

Before restoring any hot or unfavorable outer plasma, rerun the redesigned box
with uniform temperature everywhere. The purpose is to answer:

- can the new geometry produce a real, delayed outer response at all?

### 8C. reintroduce the unfavorable outer plasma

Only after 8B passes should the PBH-analogue temperature contrast be restored:

- `backgr_T_inner = 1.0`
- `backgr_T_outer = 4.0`

Then compare the propagation lag, front morphology, and outward flux against the
8B control.

## Success criterion for moving beyond 8C

Do not treat the next phase as production-ready unless the redesigned setup shows:

- positive best-lag larger than a one-dump ambiguity
- delayed outer peaks relative to the inner source region
- non-unity zero-lag correlation
- visible inward-to-outward or source-to-wing progression in `E_perp(x,t)`
- outward interface flux with a stable directional interpretation

If those conditions are met, the next step can become a real paper-level scan.
If they are not met, the project still has not resolved the PBH propagation
question.
