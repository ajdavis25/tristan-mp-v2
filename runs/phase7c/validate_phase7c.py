#!/usr/bin/env python3
"""
Phase 7C offline diagnostics from standard field dumps.

Inputs:
  python validate_phase7c.py <outdir> --input input_phase7c_demo \
    [--save phase7c_demo.png] [--csv phase7c_demo_summary.csv]
"""

import argparse
import csv
import glob
import os
import re
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np


class FieldReadError(RuntimeError):
    pass


def extract_step(path):
    base = os.path.basename(path)
    for tok in reversed(base.split(".")):
        if tok.isdigit():
            return int(tok)
    return None


def _h5_datasets(h5):
    names = []

    def _visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            names.append(name)

    h5.visititems(_visit)
    return sorted(names)


def _read_dataset(h5, name, *, file_path, required=True):
    candidates = [
        name,
        name.lower(),
        name.upper(),
        f"flds/{name}",
        f"flds/{name.lower()}",
        f"flds/{name.upper()}",
        f"fields/{name}",
        f"fields/{name.lower()}",
        f"fields/{name.upper()}",
    ]
    for key in candidates:
        try:
            obj = h5[key]
        except KeyError:
            continue
        if not isinstance(obj, h5py.Dataset):
            continue
        arr = np.asarray(obj[()])
        if arr.size == 0:
            raise FieldReadError(f"{file_path}: dataset '{key}' empty")
        if not np.issubdtype(arr.dtype, np.number):
            raise FieldReadError(f"{file_path}: dataset '{key}' non-numeric dtype {arr.dtype}")
        return arr
    if required:
        avail = ", ".join(_h5_datasets(h5))
        raise FieldReadError(f"{file_path}: missing '{name}'. Available: {avail}")
    return None


def parse_input_file(path):
    values = {}
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.split("#", 1)[0].strip()
            if not line or line.startswith("<"):
                continue
            if "=" not in line:
                continue
            key, value = [part.strip() for part in line.split("=", 1)]
            values[key] = value
    return values


def input_float(values, key, default):
    try:
        return float(values.get(key, default))
    except ValueError as exc:
        raise RuntimeError(f"Could not parse float for '{key}' from input file") from exc


def build_config(input_path):
    values = parse_input_file(input_path)
    inner_xmin = input_float(values, "inner_xmin_frac", 0.30)
    inner_xmax = input_float(values, "inner_xmax_frac", 0.70)
    transition = input_float(values, "transition_width_frac", 0.10)
    if inner_xmin <= 0.0 or inner_xmax >= 1.0 or inner_xmin >= inner_xmax:
        raise RuntimeError("Invalid inner_xmin_frac / inner_xmax_frac in input file")
    if transition <= 0.0:
        raise RuntimeError("transition_width_frac must be > 0")
    if inner_xmin - transition < 0.0 or inner_xmax + transition > 1.0:
        raise RuntimeError("Transition region extends outside the box")
    return {
        "inner_xmin_frac": inner_xmin,
        "inner_xmax_frac": inner_xmax,
        "transition_width_frac": transition,
        "zone_x0_frac": inner_xmin - transition,
        "zone_x1_frac": inner_xmin,
        "zone_x2_frac": inner_xmax,
        "zone_x3_frac": inner_xmax + transition,
    }


def read_profiles(outdir):
    candidates = sorted(
        path
        for path in glob.glob(os.path.join(outdir, "flds", "flds.tot.*"))
        if extract_step(path) is not None and not path.endswith(".xdmf")
    )

    used_files = []
    times = []
    eperp_x = []
    stress_x = []
    sx_x = []

    for path in candidates:
        try:
            with h5py.File(path, "r") as h5:
                bx = _read_dataset(h5, "bx", file_path=path)
                by = _read_dataset(h5, "by", file_path=path)
                bz = _read_dataset(h5, "bz", file_path=path)
                ey = _read_dataset(h5, "ey", file_path=path)
                ez = _read_dataset(h5, "ez", file_path=path)
        except OSError:
            continue

        if bx.ndim != 3:
            raise RuntimeError(f"{path}: expected 3D field arrays, got shape {bx.shape}")

        times.append(float(extract_step(path)))
        eperp_x.append(np.mean(0.5 * (bx * bx + by * by), axis=(1, 2), dtype=np.float64))
        stress_x.append(np.mean(-(bx * by), axis=(1, 2), dtype=np.float64))
        sx_x.append(np.mean(ey * bz - ez * by, axis=(1, 2), dtype=np.float64))
        used_files.append(path)

    if not used_files:
        raise RuntimeError(f"No readable HDF5 field dumps in {outdir}/flds")

    eperp_x = np.asarray(eperp_x, dtype=float)
    nx = eperp_x.shape[1]
    x_frac = (np.arange(nx, dtype=float) + 0.5) / float(nx)

    return {
        "times": np.asarray(times, dtype=float),
        "x_frac": x_frac,
        "eperp_x": eperp_x,
        "stress_x": np.asarray(stress_x, dtype=float),
        "sx_x": np.asarray(sx_x, dtype=float),
        "files": used_files,
    }


def region_mask(x_frac, xmin, xmax, *, include_right=True):
    if include_right:
        return (x_frac >= xmin) & (x_frac <= xmax)
    return (x_frac >= xmin) & (x_frac < xmax)


def interface_mask(nx, center_frac, half_width=2):
    center_idx = int(np.clip(np.rint(center_frac * nx - 0.5), 0, nx - 1))
    lo = max(0, center_idx - half_width)
    hi = min(nx - 1, center_idx + half_width)
    mask = np.zeros(nx, dtype=bool)
    mask[lo : hi + 1] = True
    return mask


def onset_time(times, series, factor=3.0):
    series = np.asarray(series, dtype=float)
    times = np.asarray(times, dtype=float)
    if len(series) < 2 or len(times) < 2:
        return np.nan

    baseline = float(series[0])
    if not np.isfinite(baseline):
        return np.nan

    threshold = factor * baseline
    if threshold <= baseline:
        idx = np.where(series[1:] > baseline)[0] + 1
    else:
        idx = np.where(series[1:] >= threshold)[0] + 1
    if idx.size == 0:
        return np.nan
    return float(times[int(idx[0])])


def peak_time(times, series):
    if not np.isfinite(series).any():
        return np.nan
    return float(times[int(np.nanargmax(series))])


def pearson_corr(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    mask = np.isfinite(a) & np.isfinite(b)
    if np.count_nonzero(mask) < 5:
        return np.nan

    a = a[mask]
    b = b[mask]
    a = a - np.mean(a)
    b = b - np.mean(b)

    denom = np.sqrt(np.sum(a * a) * np.sum(b * b))
    if denom <= 0.0:
        return np.nan
    return float(np.sum(a * b) / denom)


def best_nonnegative_lag(times, lead_series, lag_series, *, min_overlap_frac=0.5, min_overlap_points=8):
    if len(times) < min_overlap_points:
        return np.nan, np.nan

    lead = np.asarray(lead_series, dtype=float)
    lag = np.asarray(lag_series, dtype=float)
    mask = np.isfinite(lead) & np.isfinite(lag)
    if np.count_nonzero(mask) < min_overlap_points:
        return np.nan, np.nan

    lead = lead[mask]
    lag = lag[mask]
    t = np.asarray(times, dtype=float)[mask]
    if len(t) < min_overlap_points:
        return np.nan, np.nan

    dt = float(np.median(np.diff(t))) if len(t) > 1 else 1.0
    min_overlap = max(min_overlap_points, int(np.ceil(min_overlap_frac * len(t))))
    max_lag_idx = len(t) - min_overlap
    if max_lag_idx < 0:
        return np.nan, np.nan

    best_corr = -np.inf
    best_lag = np.nan
    best_idx = None
    tol = 1.0e-6
    for lag_idx in range(max_lag_idx + 1):
        overlap = len(lead) - lag_idx
        corr = pearson_corr(lead[:overlap], lag[lag_idx:lag_idx + overlap])
        if not np.isfinite(corr):
            continue
        if corr > best_corr + tol or (abs(corr - best_corr) <= tol and (best_idx is None or lag_idx < best_idx)):
            best_corr = corr
            best_lag = lag_idx * dt
            best_idx = lag_idx

    if not np.isfinite(best_corr):
        return np.nan, np.nan
    return float(best_lag), float(best_corr)


def safe_mean(values):
    values = np.asarray(values, dtype=float)
    return float(np.nanmean(values)) if np.isfinite(values).any() else np.nan


def summarize_run(data, cfg):
    times = data["times"]
    x_frac = data["x_frac"]
    eperp_x = data["eperp_x"]
    stress_x = data["stress_x"]
    sx_x = data["sx_x"]
    nx = len(x_frac)

    inner_mask = region_mask(x_frac, cfg["inner_xmin_frac"], cfg["inner_xmax_frac"])
    left_outer_mask = region_mask(x_frac, 0.0, cfg["zone_x0_frac"], include_right=False)
    right_outer_mask = region_mask(x_frac, cfg["zone_x3_frac"], 1.0)
    left_interface_mask = interface_mask(nx, cfg["inner_xmin_frac"], half_width=2)
    right_interface_mask = interface_mask(nx, cfg["inner_xmax_frac"], half_width=2)

    if not inner_mask.any() or not left_outer_mask.any() or not right_outer_mask.any():
        raise RuntimeError("One or more region masks are empty; check input fractions and nx")

    e_inner = np.mean(eperp_x[:, inner_mask], axis=1)
    e_outer_left = np.mean(eperp_x[:, left_outer_mask], axis=1)
    e_outer_right = np.mean(eperp_x[:, right_outer_mask], axis=1)

    sx_out_left = -np.mean(sx_x[:, left_interface_mask], axis=1)
    sx_out_right = np.mean(sx_x[:, right_interface_mask], axis=1)
    sx_out_avg = 0.5 * (sx_out_left + sx_out_right)

    t_onset_inner = onset_time(times, e_inner, factor=3.0)
    t_peak_inner = peak_time(times, e_inner)
    t_peak_outer_left = peak_time(times, e_outer_left)
    t_peak_outer_right = peak_time(times, e_outer_right)

    peak_delay_left = t_peak_outer_left - t_peak_inner if np.isfinite(t_peak_inner) and np.isfinite(t_peak_outer_left) else np.nan
    peak_delay_right = t_peak_outer_right - t_peak_inner if np.isfinite(t_peak_inner) and np.isfinite(t_peak_outer_right) else np.nan

    if np.isfinite(t_onset_inner):
        post_mask = times >= t_onset_inner
    else:
        post_mask = np.zeros_like(times, dtype=bool)

    zero_lag_corr_left = pearson_corr(e_inner[post_mask], e_outer_left[post_mask]) if np.any(post_mask) else np.nan
    zero_lag_corr_right = pearson_corr(e_inner[post_mask], e_outer_right[post_mask]) if np.any(post_mask) else np.nan
    best_lag_left, corr_peak_left = best_nonnegative_lag(times[post_mask], e_inner[post_mask], e_outer_left[post_mask])
    best_lag_right, corr_peak_right = best_nonnegative_lag(times[post_mask], e_inner[post_mask], e_outer_right[post_mask])

    mean_outward_flux_left = safe_mean(sx_out_left[post_mask]) if np.any(post_mask) else np.nan
    mean_outward_flux_right = safe_mean(sx_out_right[post_mask]) if np.any(post_mask) else np.nan
    mean_outward_flux_avg = safe_mean(sx_out_avg[post_mask]) if np.any(post_mask) else np.nan

    summary = {
        "t_onset_inner": t_onset_inner,
        "t_peak_inner": t_peak_inner,
        "t_peak_outer_left": t_peak_outer_left,
        "t_peak_outer_right": t_peak_outer_right,
        "peak_delay_left": peak_delay_left,
        "peak_delay_right": peak_delay_right,
        "zero_lag_corr_left": zero_lag_corr_left,
        "zero_lag_corr_right": zero_lag_corr_right,
        "best_lag_left": best_lag_left,
        "best_lag_right": best_lag_right,
        "corr_peak_left": corr_peak_left,
        "corr_peak_right": corr_peak_right,
        "mean_outward_flux_left": mean_outward_flux_left,
        "mean_outward_flux_right": mean_outward_flux_right,
        "mean_outward_flux_avg": mean_outward_flux_avg,
    }

    series = {
        "e_inner": e_inner,
        "e_outer_left": e_outer_left,
        "e_outer_right": e_outer_right,
        "sx_out_left": sx_out_left,
        "sx_out_right": sx_out_right,
        "sx_out_avg": sx_out_avg,
        "left_outer_mask": left_outer_mask,
        "right_outer_mask": right_outer_mask,
        "inner_mask": inner_mask,
    }
    return summary, series


def write_csv(path, summary):
    fieldnames = [
        "t_onset_inner",
        "t_peak_inner",
        "t_peak_outer_left",
        "t_peak_outer_right",
        "peak_delay_left",
        "peak_delay_right",
        "zero_lag_corr_left",
        "zero_lag_corr_right",
        "best_lag_left",
        "best_lag_right",
        "corr_peak_left",
        "corr_peak_right",
        "mean_outward_flux_left",
        "mean_outward_flux_right",
        "mean_outward_flux_avg",
    ]
    with open(path, "w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({name: summary.get(name) for name in fieldnames})


def add_zone_lines(ax, cfg):
    for xpos in (cfg["zone_x0_frac"], cfg["zone_x1_frac"], cfg["zone_x2_frac"], cfg["zone_x3_frac"]):
        ax.axvline(xpos, color="white", linestyle="--", linewidth=1.0, alpha=0.9)


def plot_results(data, cfg, summary, series, savepath):
    times = data["times"]
    x_frac = data["x_frac"]

    fig, axs = plt.subplots(2, 2, figsize=(12, 9))
    ax_e, ax_stress, ax_ts, ax_flux = axs.flatten()

    eps = np.finfo(float).tiny
    im_e = ax_e.imshow(
        np.log10(np.clip(data["eperp_x"], eps, None)),
        origin="lower",
        aspect="auto",
        extent=[0.0, 1.0, times[0], times[-1]],
        cmap="magma",
    )
    add_zone_lines(ax_e, cfg)
    ax_e.set_title(r"$\log_{10} E_\perp(x,t)$")
    ax_e.set_xlabel(r"$x/L_x$")
    ax_e.set_ylabel("time")
    fig.colorbar(im_e, ax=ax_e, fraction=0.046, pad=0.04)

    stress_lim = np.nanmax(np.abs(data["stress_x"]))
    im_s = ax_stress.imshow(
        data["stress_x"],
        origin="lower",
        aspect="auto",
        extent=[0.0, 1.0, times[0], times[-1]],
        cmap="RdBu_r",
        vmin=-stress_lim,
        vmax=stress_lim,
    )
    add_zone_lines(ax_stress, cfg)
    ax_stress.set_title(r"$\langle -B_x B_y \rangle_{y,z}(x,t)$")
    ax_stress.set_xlabel(r"$x/L_x$")
    ax_stress.set_ylabel("time")
    fig.colorbar(im_s, ax=ax_stress, fraction=0.046, pad=0.04)

    ax_ts.plot(times, series["e_inner"], label="inner")
    ax_ts.plot(times, series["e_outer_left"], label="outer left")
    ax_ts.plot(times, series["e_outer_right"], label="outer right")
    if np.isfinite(summary["t_onset_inner"]):
        ax_ts.axvline(summary["t_onset_inner"], color="k", linestyle="--", linewidth=1.0, label="inner onset")
    if np.isfinite(summary["t_peak_inner"]):
        ax_ts.scatter([summary["t_peak_inner"]], [np.nanmax(series["e_inner"])], color="C0", zorder=3)
    if np.isfinite(summary["t_peak_outer_left"]):
        ax_ts.scatter([summary["t_peak_outer_left"]], [np.nanmax(series["e_outer_left"])], color="C1", zorder=3)
    if np.isfinite(summary["t_peak_outer_right"]):
        ax_ts.scatter([summary["t_peak_outer_right"]], [np.nanmax(series["e_outer_right"])], color="C2", zorder=3)
    ax_ts.set_yscale("log")
    ax_ts.set_title(r"region-averaged $E_\perp(t)$")
    ax_ts.set_xlabel("time")
    ax_ts.set_ylabel(r"$E_\perp$")
    ax_ts.legend(fontsize=9)
    ax_ts.grid(True, alpha=0.3)

    ax_flux.plot(times, series["sx_out_left"], label="left outward flux")
    ax_flux.plot(times, series["sx_out_right"], label="right outward flux")
    ax_flux.plot(times, series["sx_out_avg"], label="avg outward flux")
    ax_flux.axhline(0.0, color="k", linewidth=1.0)
    if np.isfinite(summary["t_onset_inner"]):
        ax_flux.axvline(summary["t_onset_inner"], color="k", linestyle="--", linewidth=1.0)
    ax_flux.set_title(r"interface $S_x(t)$ (outward convention)")
    ax_flux.set_xlabel("time")
    ax_flux.set_ylabel(r"$S_x$")
    ax_flux.legend(fontsize=9)
    ax_flux.grid(True, alpha=0.3)

    fig.suptitle("phase 7C: two-zone propagation demonstrator")
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    fig.text(
        0.5,
        0.01,
        "Zone lines mark outer core / transition / inner core boundaries at x/Lx = 0.20, 0.30, 0.70, 0.80.",
        ha="center",
        fontsize=9,
    )

    if savepath:
        fig.savefig(savepath, dpi=160)
        print(f"wrote: {savepath}")
    else:
        plt.show()


def build_parser():
    parser = argparse.ArgumentParser(description="Offline diagnostics for the Phase 7C two-zone MRI run.")
    parser.add_argument("outdir", help="Path to the Phase 7C output directory")
    parser.add_argument("--input", required=True, dest="input_path", help="Input file used for the run")
    parser.add_argument("--save", dest="savepath", default=None, help="Optional PNG output path")
    parser.add_argument("--csv", dest="csvpath", default=None, help="Optional CSV summary path")
    return parser


def main(argv):
    args = build_parser().parse_args(argv)
    cfg = build_config(args.input_path)
    data = read_profiles(args.outdir)
    summary, series = summarize_run(data, cfg)

    print("\nphase 7C summary")
    print(f"  t_onset_inner      = {summary['t_onset_inner']}")
    print(f"  t_peak_inner       = {summary['t_peak_inner']}")
    print(f"  t_peak_outer_left  = {summary['t_peak_outer_left']}")
    print(f"  t_peak_outer_right = {summary['t_peak_outer_right']}")
    print(f"  peak_delay_left    = {summary['peak_delay_left']}")
    print(f"  peak_delay_right   = {summary['peak_delay_right']}")
    print(f"  zero_lag_corr_left = {summary['zero_lag_corr_left']}")
    print(f"  zero_lag_corr_right= {summary['zero_lag_corr_right']}")
    print(f"  best_lag_left      = {summary['best_lag_left']}")
    print(f"  best_lag_right     = {summary['best_lag_right']}")
    print(f"  corr_peak_left     = {summary['corr_peak_left']}")
    print(f"  corr_peak_right    = {summary['corr_peak_right']}")
    print(f"  mean_outward_flux_left  = {summary['mean_outward_flux_left']}")
    print(f"  mean_outward_flux_right = {summary['mean_outward_flux_right']}")
    print(f"  mean_outward_flux_avg   = {summary['mean_outward_flux_avg']}")

    if args.csvpath:
        write_csv(args.csvpath, summary)
        print(f"wrote: {args.csvpath}")

    plot_results(data, cfg, summary, series, args.savepath)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
