# MRI Persistence Recovery From Phase 8A-long

This document records the current project truth from the repo, evaluates the
original JPP submission and referee report against that truth, and defines the
shortest defensible path from the current state to publishable evidence about
whether MRI persists in the early universe.

## Bottom line

- `phase8a_long` is another propagation-isolation null result, not progress
  toward the PBH-propagation claim.
- The decisive evidence is in
  [runs/phase8/phase8a_long_control_summary.csv](/work/vmo703/tristan-mp-v2/runs/phase8/phase8a_long_control_summary.csv):
  - `peak_delay_left = 0`, `peak_delay_right = 0`
  - `best_lag_left = 0`, `best_lag_right = 0`
  - `zero_lag_corr_left = 0.9999965209199224`
  - `zero_lag_corr_right = 0.9999956346050219`
  - `mean_outward_flux_avg = -8.1659727194631e-10`
- The project currently supports a **local MRI-persistence** story much better
  than a **PBH outward-propagation** story.
- The original paper should not be treated as fully rescuable in its current
  claim scope. The referee's PIC objection is only partially addressable by the
  repo as it stands today.

## Current project state

The original roadmap in
[notebooks/assets/map.md](/work/vmo703/tristan-mp-v2/notebooks/assets/map.md)
is useful historical context only. It is no longer the authoritative roadmap.

### Phase status

- Phase 6: completed robustness sweep.
  - Evidence:
    [runs/phase6/phase6.png](/work/vmo703/tristan-mp-v2/runs/phase6/phase6.png)
  - Interpretation: local MRI-like response is not a single-point artifact in
    `ppc`/`B0` space.
- Phase 7A: `PASS`.
  - Evidence:
    [runs/phase7/validate_phase7a.png](/work/vmo703/tristan-mp-v2/runs/phase7/validate_phase7a.png)
  - Interpretation: MRI-like growth persists across `Theta = 0.5-4.0` in the
    tested local pair-plasma setup.
- Phase 7B: `PROVISIONAL`.
  - Evidence:
    [runs/phase7/validate_phase7b.png](/work/vmo703/tristan-mp-v2/runs/phase7/validate_phase7b.png),
    [runs/phase7/validate_phase7b_summary.csv](/work/vmo703/tristan-mp-v2/runs/phase7/validate_phase7b_summary.csv)
  - Interpretation: morphology, coherence, and transport shift with
    temperature, but the regime-boundary claim is not fully closed because the
    particle diagnostics were undersampled in the linear window.
- Phase 7C: implemented, but failed as a PBH-propagation isolate.
  - Evidence:
    [runs/phase7c/phase7c_demo_summary.csv](/work/vmo703/tristan-mp-v2/runs/phase7c/phase7c_demo_summary.csv),
    [runs/phase7c/phase7c_control_uniformT_summary.csv](/work/vmo703/tristan-mp-v2/runs/phase7c/phase7c_control_uniformT_summary.csv)
  - Interpretation: localized seed plus hot outer wings was still a
    near-global response.
- Phase 7D: stronger null result than 7C.
  - Evidence:
    [runs/phase7d/phase7d_demo_summary.csv](/work/vmo703/tristan-mp-v2/runs/phase7d/phase7d_demo_summary.csv)
  - Interpretation: localizing the drive as well as the seed still did not
    produce delayed outer response.
- Phase 8A: source-plus-buffer control failed to isolate delayed outward
  response.
  - Evidence:
    [runs/phase8/phase8a_control_summary.csv](/work/vmo703/tristan-mp-v2/runs/phase8/phase8a_control_summary.csv)
  - Interpretation: adding absorbing x-edge buffers did not break the
    synchronous behavior.
- Phase 8A-long: same failure at larger `Lx`.
  - Evidence:
    [runs/phase8/phase8a_long_control_summary.csv](/work/vmo703/tristan-mp-v2/runs/phase8/phase8a_long_control_summary.csv),
    [runs/phase8/phase8a_long_control.png](/work/vmo703/tristan-mp-v2/runs/phase8/phase8a_long_control.png)
  - Interpretation: this closes the "just make the box longer" branch.

### What the current repo already supports

- A validated local PIC workflow for MRI-capable shearing-box pair-plasma runs.
- A robust local claim that MRI-like linear growth persists over the tested
  temperature range.
- A weaker local claim that hotter / more marginal runs become less MHD-like in
  coherence and transport.

### What the current repo does not support

- A validated PBH-analogue claim that MRI launched near the source region
  propagates outward into surrounding early-universe plasma.
- A paper-safe claim that the original GRMHD dark-matter argument has already
  been kinetically validated.

## Paper and referee assessment

### Original paper

The proof PDF in
[notebooks/docs/59 CurdAnantuaLujanFowler - Journal of Plasma Physics - PLA-2025-0146_Proof_hi.pdf](/work/vmo703/tristan-mp-v2/notebooks/docs/59%20CurdAnantuaLujanFowler%20-%20Journal%20of%20Plasma%20Physics%20-%20PLA-2025-0146_Proof_hi.pdf)
argues that MRI-driven PBH accretion in the Positronium Era can account for a
significant dark-matter contribution.

That submitted claim is stronger than the current PIC branch can justify.

### Referee report

The referee report in
[notebooks/docs/Review report on JPP paper entitled “MRI driven BH as DM”.pdf](/work/vmo703/tristan-mp-v2/notebooks/docs/Review%20report%20on%20JPP%20paper%20entitled%20%E2%80%9CMRI%20driven%20BH%20as%20DM%E2%80%9D.pdf)
raises five main concerns:

1. Accretion vs capture.
   - This is still not resolved by current PIC work.
   - It remains outside what the present local PIC program can prove.
2. MRI region extending beyond the Bondi radius.
   - This is not resolved by the current 7/8 local-box PIC runs.
   - It remains part of the stronger PBH/macroscopic transport claim.
3. Collisionality vs MRI.
   - This is only partially addressable with current work.
   - It becomes meaningfully addressable only if the local PIC evidence package
     is completed cleanly.
4. Plasma validity / `N_D ~ O(10)`.
   - This is the decisive point the PIC project can answer directly.
   - This should be treated as the primary near-term scientific target.
5. "Positronium Era" terminology.
   - This is an editorial issue, not a gating simulation issue.

### Minimal defensible next manuscript claim

If the local PIC workstream is completed, the narrow defensible claim is:

- **Local MRI-like growth persists in an early-universe-analogue marginal pair
  plasma over part of the relevant parameter space.**

The current repo does **not** yet support the stronger claim:

- **MRI launched near a PBH propagates outward through the surrounding
  early-universe plasma.**

## Stage 1: referee-proof local MRI persistence package

This is the primary objective. It answers whether MRI persists locally in the
early-universe analogue.

### Deliverables

- Close Phase 7B with one final, stable analysis package.
- Add one compact bridge sweep connecting the current temperature sweep to lower
  collectivity.
- Add an explicit cosmology-to-PIC interpretation section.

### Canonical local setup

- Keep the Phase 7 local shearing-box branch as the canonical local-persistence
  setup.
- Use the existing `user_mri_6` / Phase 7 geometry and diagnostics.
- Do not use the propagation branch for the local-persistence claim.

### Required analysis closure

Finish the Phase 7B analysis package with one final figure that combines:

- existing temperature-persistence evidence from
  [runs/phase7/validate_phase7a.png](/work/vmo703/tristan-mp-v2/runs/phase7/validate_phase7a.png)
- field-only coherence or spectral metric from the current dumps
- particle anisotropy from new high-cadence particle reruns

Run two particle-enabled reruns with higher `prtl_write_every` to fix the
current 7B diagnostic gap:

- clean reference case: `ppc = 16`, `Theta = 1`, `B0 = 0.06`
- marginal hot case: `ppc = 8`, `Theta = 4`, `B0 = 0.06`

### Required bridge sweep

Add a compact field-only bridge sweep to close the current
low-collectivity/high-temperature gap:

- `ppc = 4, 8, 16`
- `Theta = 1, 4`
- `B0 = 0.06`
- same Phase 7 box and runtime conventions

This sweep is meant to answer the referee-facing question the current repo still
leaves exposed: whether the local MRI-like response persists once temperature
and collectivity are stressed together, not separately.

### Required interpretation layer

Produce a referee-facing interpretation table in this README that maps:

- physical early-universe quantities from the referee report:
  - `T`
  - `N_D`
  - collisionality concern
- to the PIC surrogate knobs already used in the repo:
  - `Theta`
  - `ppc`
  - coherence diagnostics
  - stress/transport diagnostics

### Stage 1 acceptance criteria

- local MRI-like exponential growth is still present in at least part of the
  bridge sweep
- the result is robust against the already completed Phase 6 `ppc`/`B` scan and
  the existing Phase 7 temperature scan
- the final 7B figure supports the sentence:
  - MRI-like growth persists, but coherence and transport become less MHD-like
    as the plasma becomes more marginal and hotter
- the README can cleanly state what the local PIC program does and does not
  prove

## Stage 2: PBH-propagation program

This is a separate, stronger claim. It should not be mixed into the Stage 1
success criterion.

### Current propagation conclusion

7C, 7D, 8A, and 8A-long all failed because the centered, symmetric source-box
design still excites near-global response.

Do **not** continue the current centered-box escalation beyond `8A-long`. The
next propagation experiment must be a design change, not another box-size
increase.

### Next propagation design

Replace the centered symmetric source with a **one-sided source-to-medium**
experiment:

- source region occupies the left side of the domain, not the center
- only the source region is driven and seeded
- the medium to the right is passive and is the only propagation target
- the far-right boundary is absorbing
- the left side is treated as source support, not as a second passive wing

Diagnostics must also become one-directional:

- right-going `Sx`
- source-to-right lag
- right-side peak delay
- visible front speed in `E_perp(x,t)`

### Run order

1. one-sided uniform-temperature control
2. one-sided unfavorable-medium run with `backgr_T_outer > backgr_T_inner`
3. only after both work, any production scan over outer conditions

### Stage 2 acceptance criteria

- `best_lag_right` positive by more than one output interval
- `peak_delay_right` positive by more than one output interval
- `zero_lag_corr_right < 0.99`
- net positive right-going source-interface flux over the propagation window
- visible source-to-right front in `E_perp(x,t)`

## Manuscript strategy

If Stage 1 succeeds and Stage 2 is still unresolved, the correct near-term
product is a **PIC persistence paper / referee-response evidence package**, not
a full claim that PBH MRI accretion in the early universe has been validated.

The original JPP dark-matter paper can only be revised around the current
evidence if its claims are narrowed.

A full PBH-propagation or dark-matter-growth claim should wait until the
one-sided propagation design passes.

## Test and evidence checklist

Before claiming "MRI persists in the early universe," the project should have:

- existing: Phase 6 robustness figure
  - [runs/phase6/phase6.png](/work/vmo703/tristan-mp-v2/runs/phase6/phase6.png)
- existing: Phase 7A temperature-persistence figure
  - [runs/phase7/validate_phase7a.png](/work/vmo703/tristan-mp-v2/runs/phase7/validate_phase7a.png)
- required: final Phase 7B coherence/anisotropy figure from corrected
  high-cadence particle data
- required: compact bridge sweep at lower collectivity and higher temperature
- required: written cosmology-to-PIC mapping table

Explicitly **not** required for Stage 1:

- a positive PBH outward-propagation result

Explicitly required for the stronger PBH claim:

- successful one-sided propagation control
- successful one-sided unfavorable-medium follow-up

## Defaults and assumptions

- Strategy: `Two-stage plan`
  - Stage 1 = referee-proof local MRI persistence
  - Stage 2 = stronger PBH propagation claim
- Existing 7C/7D/8A/8A-long runs are treated as completed null results, not
  branches to continue scaling
- [notebooks/assets/map.md](/work/vmo703/tristan-mp-v2/notebooks/assets/map.md)
  is retained as an archival reference, not the governing roadmap
