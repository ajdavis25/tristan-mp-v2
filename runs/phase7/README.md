# Phase 7 Notes

This directory contains the late Phase 7 restart point:

- `validate_phase7a.py` / `validate_phase7a.png`: field-only temperature sweep
- `validate_phase7b.py` / `validate_phase7b.png`: particle-enabled boundary diagnostics

## Sweep layout

- Phase 7A field-only runs:
  - `out_ppc08_B006_T05`, `out_ppc08_B006_T1`, `out_ppc08_B006_T2`, `out_ppc08_B006_T4`
  - `out_ppc16_B006_T05`, `out_ppc16_B006_T1`, `out_ppc16_B006_T2`, `out_ppc16_B006_T4`
- Phase 7B particle-enabled runs:
  - `out_ppc08_T05`, `out_ppc08_T1`, `out_ppc08_T2`, `out_ppc08_T4`
  - `out_ppc16_T05`, `out_ppc16_T1`, `out_ppc16_T2`, `out_ppc16_T4`

All current Phase 7 inputs use:

- `B0_ext = 0.06`
- `Omega = 0.004`
- `q_shear = 1.5`
- `last = 7000`

## Build provenance

- Active 3D build cache: `/work/vmo703/tristan-mp-v2/build_3d/CMakeCache.txt`
- Cached user target: `user_mri_6`
- `/work/vmo703/tristan-mp-v2/user/user_mri_6.F90` and `/work/vmo703/tristan-mp-v2/user/user_mri_7.F90` are identical in the current tree, so the 6 vs 7 distinction is naming/provenance rather than physics.

## Log provenance

- Field-only Phase 7A logs are clean and map one-to-one to current outputs:
  - `/work/vmo703/tristan-mp-v2/logs/phase7_181754` through `/work/vmo703/tristan-mp-v2/logs/phase7_181761`
- Particle-enabled Phase 7B provenance is split:
  - clean/current naming:
    - `/work/vmo703/tristan-mp-v2/logs/phase7b_181804` -> `out_ppc16_T05`
    - `/work/vmo703/tristan-mp-v2/logs/phase7b_181805` -> `out_ppc08_T1`
    - `/work/vmo703/tristan-mp-v2/logs/phase7b_181806` -> `out_ppc08_T05`
    - `/work/vmo703/tristan-mp-v2/logs/phase7b_181807` -> `out_ppc16_T2`
  - older/ambiguous naming:
    - `/work/vmo703/tristan-mp-v2/logs/phase7b_181799` -> current `out_ppc08_T2`
    - `/work/vmo703/tristan-mp-v2/logs/phase7b_181800` -> current `out_ppc08_T4`
    - `/work/vmo703/tristan-mp-v2/logs/phase7b_181801` -> current `out_ppc16_T1`
    - `/work/vmo703/tristan-mp-v2/logs/phase7b_181802` -> current `out_ppc16_T4`

Those older logs refer to `_prtl` output names and echo mismatched log-directory names inside the files, so treat them as evidence of successful runs but not as clean provenance.

## Analysis conventions

- `validate_phase7a.py` keeps the original field-only sweep metrics.
- `validate_phase7b.py` uses the full-run `corr(E_perp, |<B_x B_y>|)` as the canonical transport metric because that matches the saved legacy figure and notebook text.
- The same script still identifies the linear MRI window from `E_perp`, but only for field-based `k_z` summaries.
- Current particle dumps are written every 1000 steps, so no particle dumps overlap the linear MRI window (typically `t ~ 120-220`); anisotropy summaries therefore use the available sparse particle history rather than a true linear-window particle diagnostic.

## Phase 7 pass / fail status

Current interpretation of Phase 7:

- Phase 7A: PASS
- Phase 7B: PROVISIONAL / INCOMPLETE
- Phase 7C: IMPLEMENTED, BUT FAILS AS A PBH-PROPAGATION ISOLATE
- Phase 7D: IMPLEMENTED, BUT FAILS AS A STRONGER PBH-PROPAGATION ISOLATE

Overall project status for the PBH-centered question:

- unresolved / incomplete
- the code now supports the MRI-exists and MRI-becomes-kinetic claims
- it does not yet support the PBH-analogue claim that MRI-driven activity is
  launched near the source region and then propagates outward into the
  surrounding plasma in a controlled, delayed way

Phase 7A:

- the field-only temperature sweep supports the claim that MRI-like growth persists across the tested temperature range

Phase 7B:

- the particle-enabled sweep supports the weaker claim that higher-temperature / lower-collectivity cases look less MHD-like
- it does not yet define a clean breakdown boundary

Phase 7C:

- detailed notes live in `/work/vmo703/tristan-mp-v2/runs/phase7c/README.md`
- the current localized-seed two-zone design completed and analyzed successfully
- both the hot-wing demonstrator and the uniform-temperature control give:
  - `best_lag = 0`
  - near-unity zero-lag inner/outer correlation
- therefore the current 7C configuration does not isolate inward-to-outward
  propagation and should not be treated as a validated positive result

Phase 7D:

- detailed notes live in `/work/vmo703/tristan-mp-v2/runs/phase7d/README.md`
- the shearing-box drive was localized to the inner zone to make the outer wings
  passive by construction
- the first 7D demonstrator still gives:
  - `t_peak_inner = 580`
  - `t_peak_outer_left = 580`
  - `t_peak_outer_right = 580`
  - `best_lag_left = 0`
  - `best_lag_right = 0`
  - `zero_lag_corr_left = 0.9999400970713795`
  - `zero_lag_corr_right = 0.9998983034875565`
- therefore 7D strengthens the negative conclusion from 7C: the current
  shearing-box analogue is still exciting a nearly global response rather than
  a clean delayed outward signal

## PBH-centered success criterion

The notebook makes the intended interpretation explicit:

- `plan.ipynb` phase 7C: "provide a controlled surrogate for MRI propagation
  away from a black hole"
- `plan.ipynb` follow-up questions: "does MRI propagate outward or die outside
  the black hole?"
- later notebook note: "mri propagates all the way out to r_0 but it starts by
  the black hole"

Given that framing, Phase 7 propagation work is not complete just because the
plots exist. For the PBH-analogue branch to count as a success, the surrogate
experiment must show all of the following:

- an inner source region that becomes active first
- an outer region that responds later rather than synchronously
- a positive best-lag signal that survives comparison to controls
- a visible inner-to-outer front or delayed outer brightening in `E_perp(x,t)`
- an outward-flux interpretation that is stronger than equal-and-opposite
  interface oscillations

The existing 7C and 7D runs fail that criterion.

## Concrete next phase: proposed Phase 8

Because 7C and 7D both fail for the same reason, the next phase should not be a
parameter tweak. It should be a new propagation-isolate design whose purpose is
to break the current near-global mode behavior.

### Phase 8A. source-plus-buffer isolate

Create a new user problem, likely `user/user_mri_8a.F90`, starting from
`user/user_mri_7d.F90` and using the existing user hooks already available in
this codebase:

- keep the inner localized seed and inner-only drive
- enlarge the radial direction beyond the current `sizex = 8`
- keep the first test uniform in temperature so the propagation mechanism is
  isolated before reintroducing unfavorable outer plasma
- add explicit edge buffer layers near the `x` boundaries
- implement damping in:
  - `userFieldBoundaryConditions()`
  - `userParticleBoundaryConditions()`

Grounding for feasibility:

- those hooks already exist in the MRI userfiles
- the repository already contains absorbing-layer logic examples in
  `user/user_psr_simple.F90`

The goal of 8A is simple: suppress recirculation / standing-mode behavior so
that a source-launched signal can separate from the boundaries.

### Phase 8B. source-only control

Under `runs/phase8`, first run a uniform-temperature source-only control:

- inner-only drive
- inner-only seed
- passive outer wings
- no temperature penalty yet

This run should answer one question only: can the redesigned box produce a
nonzero propagation lag when the outer plasma is not also made unfavorable?

### Phase 8C. hot-wing PBH analogue

Only if 8B shows real propagation, reintroduce the unfavorable outer wings:

- `backgr_T_inner = 1.0`
- `backgr_T_outer = 4.0`
- same geometry and same damping/buffer setup as 8B

This is the first configuration that should be interpreted as the PBH analogue:
MRI launched near the source region, then tested for survival as it moves into
less favorable plasma.

### Phase 8D. production scan

Only after 8B and 8C pass should the project move into a production-style
research phase:

- vary outer temperature or outer collectivity
- vary source width and outer-wing width
- vary box length in `x`
- keep the same propagation metrics and controls

At that point, the result is closer to paper-grade because the experiment
definition has already passed its control test.

## Phase 8 pass criteria

Treat the next propagation design as a success only if all of the following are
true in the redesigned setup:

- `best_lag_left` or `best_lag_right` is positive in a way that is larger than a
  one-dump ambiguity
- the outer peak time is delayed relative to the inner peak
- zero-lag correlation is no longer effectively unity across the post-onset
  window
- `E_perp(x,t)` shows a visible inner-to-outer progression
- the same analysis on the control run supports the same interpretation rather
  than collapsing back to the current near-global response

## Regeneration commands

From `/work/vmo703/tristan-mp-v2/runs/phase7`:

```bash
python3 validate_phase7a.py \
  out_ppc08_B006_T05 out_ppc08_B006_T1 out_ppc08_B006_T2 out_ppc08_B006_T4 \
  out_ppc16_B006_T05 out_ppc16_B006_T1 out_ppc16_B006_T2 out_ppc16_B006_T4 \
  --save validate_phase7a.png

python3 validate_phase7b.py \
  out_ppc08_T05 out_ppc08_T1 out_ppc08_T2 out_ppc08_T4 \
  out_ppc16_T05 out_ppc16_T1 out_ppc16_T2 out_ppc16_T4 \
  --save validate_phase7b.png \
  --csv validate_phase7b_summary.csv
```
