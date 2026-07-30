#!/usr/bin/env python3
"""
3D animation of distributed tracking behavior.

What this shows:
- Observer satellites and target objects moving in 3D.
- Active tracking links observer -> target over time.
- Optional event overlays parsed from simulation log:
  Local track init/close, IOD init, DEKF init.

Usage examples:
  python3 scripts/plotting/navigation/animate_results_3d.py \
      --case no-da --show

  python3 scripts/plotting/navigation/animate_results_3d.py \
      --input /path/to/run --save /path/to/run/tracking_animation.mp4
"""

from __future__ import annotations

import argparse
import math
import re
import tempfile
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, writers


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SINGLE_OUTPUT_ROOT = (
    Path(tempfile.gettempdir())
    / "spaceagora_navigation"
    / "runs"
    / "single"
)


def repo_path(value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else REPO_ROOT / path


def case_output_name(case: str) -> str:
    return case.strip().replace("-", "_")


def numeric_col(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.nan, index=df.index)
    s = df[col]
    if pd.api.types.is_bool_dtype(s):
        return s.astype(float)
    return pd.to_numeric(s, errors="coerce")


def discover_spacecraft_ids(df: pd.DataFrame) -> list[int]:
    ids = set()
    rx = re.compile(r"^sc(\d+)_pos_1$")
    for c in df.columns:
        m = rx.match(c)
        if m:
            ids.add(int(m.group(1)))
    return sorted(ids)


def discover_target_ids(df: pd.DataFrame) -> list[int]:
    ids = set()
    rx = re.compile(r"^dekf_t(\d+)_state_obs\d+_1$")
    for c in df.columns:
        m = rx.match(c)
        if m:
            ids.add(int(m.group(1)))
    return sorted(ids)


def discover_observer_ids(df: pd.DataFrame) -> list[int]:
    ids = set()
    rx = re.compile(r"^dekf_t\d+_state_obs(\d+)_1$")
    for c in df.columns:
        m = rx.match(c)
        if m:
            ids.add(int(m.group(1)))
    return sorted(ids)


def find_target_slot_column(df: pd.DataFrame, target: int, obs: int) -> str | None:
    c = f"dekf_t{target}_slot_obs{obs}"
    if c in df.columns:
        return c
    c = f"dekf_t{target}_slot_obs{obs}_1"
    if c in df.columns:
        return c
    rx = re.compile(rf"^dekf_t{target}_slot_obs{obs}(?:_1)?$")
    for cc in df.columns:
        if rx.match(cc):
            return cc
    return None


def resolve_log_path(csv_path: Path, cli_log: str | None) -> Path | None:
    if cli_log:
        p = repo_path(cli_log)
        return p if p.exists() else None
    candidates = [
        csv_path.with_suffix(".log"),
        csv_path.parent / "run.log",
        csv_path.parent / "simulation.log",
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def parse_event_log(log_path: Path | None) -> dict[float, list[str]]:
    events = defaultdict(list)
    if log_path is None:
        return events

    rx_t = re.compile(r"\|\s*t=([0-9]+(?:\.[0-9]+)?)\s*s")
    rx_slot = re.compile(r"slot=(\d+)")
    rx_obs = re.compile(r"obs=(\d+)")

    with log_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            tm = rx_t.search(s)
            if tm is None:
                continue
            t = float(tm.group(1))
            slot = rx_slot.search(s)
            obs = rx_obs.search(s)
            slot_txt = f"slot={slot.group(1)}" if slot else ""
            obs_txt = f"obs={obs.group(1)}" if obs else ""

            if s.startswith("Local track init"):
                events[t].append(f"L-init {obs_txt} {slot_txt}".strip())
            elif s.startswith("Local track closed"):
                events[t].append(f"L-close {obs_txt} {slot_txt}".strip())
            elif s.startswith("IOD init"):
                events[t].append(f"IOD {obs_txt} {slot_txt}".strip())
            elif s.startswith("DEKF init"):
                events[t].append(f"DEKF {obs_txt} {slot_txt}".strip())

    return events


def build_aligned_df(df: pd.DataFrame, nav_dt_sec: float, nav_tol_sec: float, aligned_only: bool) -> pd.DataFrame:
    if not aligned_only:
        return df.copy().reset_index(drop=True)

    t = numeric_col(df, "time").to_numpy(dtype=float)
    phase = np.mod(t, nav_dt_sec)
    dist_to_grid = np.minimum(phase, nav_dt_sec - phase)
    aligned = np.isfinite(t) & (dist_to_grid <= nav_tol_sec)
    out = df.loc[aligned].copy().reset_index(drop=True)
    if len(out) > 1:
        out = out.iloc[:-1].copy().reset_index(drop=True)
    return out


def lerp_xyz(a: np.ndarray, b: np.ndarray, alpha: float) -> np.ndarray:
    """Linear interpolation for arrays with last dimension = 3, robust to NaNs."""
    out = np.full_like(a, np.nan, dtype=float)
    a_ok = np.all(np.isfinite(a), axis=-1)
    b_ok = np.all(np.isfinite(b), axis=-1)
    both = a_ok & b_ok
    only_a = a_ok & (~b_ok)
    only_b = b_ok & (~a_ok)
    out[both] = (1.0 - alpha) * a[both] + alpha * b[both]
    out[only_a] = a[only_a]
    out[only_b] = b[only_b]
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preset", type=str, default="default", choices=["default", "presentation"], help="Rendering preset")
    parser.add_argument("--case", type=str, default="no-da", help="Navigation case used for the default single-run input")
    parser.add_argument("--input", type=str, default=None, help="Relative or absolute folder with simulation_results.csv")
    parser.add_argument("--csv", type=str, default=None, help="Relative or absolute CSV path (overrides --input)")
    parser.add_argument("--log", type=str, default=None, help="Optional simulation log path")
    parser.add_argument("--show", action="store_true", help="Show interactive window")
    parser.add_argument("--save", type=str, default=None, help="Output animation path (.mp4 or .gif)")
    parser.add_argument("--fps", type=int, default=10)
    parser.add_argument("--stride", type=int, default=1, help="Frame decimation factor")
    parser.add_argument("--trail", type=int, default=20, help="Trail length in frames")
    parser.add_argument("--max-links-per-observer", type=int, default=0, help="Max active links shown per observer (0 = all)")
    parser.add_argument("--show-target-trails", action="store_true", default=True, help="Show true target trails (default on)")
    parser.add_argument("--hide-target-trails", action="store_true", help="Disable true target trails")
    parser.add_argument("--hide-est-trails", action="store_true", help="Disable estimated dashed trails")
    parser.add_argument("--speed-factor", type=float, default=1.5, help="Playback speed multiplier (<1 slower, >1 faster)")
    parser.add_argument("--interp-substeps", type=int, default=5, help="Interpolated subframes between samples for smoother motion")
    parser.add_argument("--camera-smooth", type=float, default=0.14, help="Camera smoothing factor in [0,1]")
    parser.add_argument("--event-focus-frames", type=int, default=18, help="Frames of temporary zoom-in after init events")
    parser.add_argument("--nav-dt", type=float, default=5.0)
    parser.add_argument("--nav-tol", type=float, default=0.05)
    parser.add_argument("--aligned-only", action="store_true", help="Use only nav-aligned samples")
    parser.add_argument("--earth", action="store_true", help="Draw Earth sphere")
    parser.add_argument("--visibility-km", type=float, default=500.0, help="Observer visibility radius used for focus/filtering")
    parser.add_argument("--focus-scale", type=float, default=1.25, help="Camera scale multiplier around observer visibility neighborhood")
    parser.add_argument("--near-margin", type=float, default=1.10, help="Margin factor on visibility radius for displaying nearby targets")
    parser.add_argument("--hide-visibility-spheres", action="store_true", help="Disable transparent visibility spheres around observers")
    args = parser.parse_args()

    if args.preset == "presentation":
        args.fps = max(args.fps, 24)
        args.trail = max(args.trail, 32)
        args.speed_factor = max(args.speed_factor, 0.75)
        args.interp_substeps = max(args.interp_substeps, 6)
        args.camera_smooth = max(args.camera_smooth, 0.20)

    input_dir = (
        repo_path(args.input)
        if args.input
        else DEFAULT_SINGLE_OUTPUT_ROOT / case_output_name(args.case)
    )
    csv_path = repo_path(args.csv) if args.csv else input_dir / "simulation_results.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    df = build_aligned_df(df, args.nav_dt, args.nav_tol, args.aligned_only)

    spacecraft = discover_spacecraft_ids(df)
    targets = discover_target_ids(df)
    observers = discover_observer_ids(df)

    if not observers or not targets:
        raise RuntimeError("Could not infer observers/targets from DEKF columns.")

    t = numeric_col(df, "time").to_numpy(dtype=float)
    n = len(df)

    obs_idx = {obs: i for i, obs in enumerate(observers)}
    tgt_idx = {tgt: i for i, tgt in enumerate(targets)}

    obs_pos = np.full((n, len(observers), 3), np.nan, dtype=float)
    tgt_pos = np.full((n, len(targets), 3), np.nan, dtype=float)
    est_pos = np.full((n, len(observers), len(targets), 3), np.nan, dtype=float)

    for obs in observers:
        cols = [f"sc{obs}_pos_1", f"sc{obs}_pos_2", f"sc{obs}_pos_3"]
        if all(c in df.columns for c in cols):
            obs_pos[:, obs_idx[obs], :] = df[cols].to_numpy(dtype=float)

    for tgt in targets:
        cols = [f"sc{tgt}_pos_1", f"sc{tgt}_pos_2", f"sc{tgt}_pos_3"]
        if all(c in df.columns for c in cols):
            tgt_pos[:, tgt_idx[tgt], :] = df[cols].to_numpy(dtype=float)

    visibility_m = max(1.0, args.visibility_km * 1e3)

    # near_mask[k, ti] = target ti is inside at least one observer visibility neighborhood at time k
    # (with configurable margin) and should be shown in the focused animation.
    delta_ot = obs_pos[:, :, None, :] - tgt_pos[:, None, :, :]
    dist_ot = np.linalg.norm(delta_ot, axis=3)
    near_mask = np.any(
        np.isfinite(dist_ot) & (dist_ot <= visibility_m * args.near_margin),
        axis=1,
    )

    # Active links: observer -> target if estimate exists (and slot exists if provided)
    active = np.zeros((n, len(observers), len(targets)), dtype=bool)
    slots = np.full((n, len(observers), len(targets)), np.nan, dtype=float)

    for obs in observers:
        oi = obs_idx[obs]
        for tgt in targets:
            ti = tgt_idx[tgt]
            est_col = f"dekf_t{tgt}_state_obs{obs}_1"
            if est_col in df.columns:
                est_ok = np.isfinite(numeric_col(df, est_col).to_numpy(dtype=float))
            else:
                est_ok = np.zeros(n, dtype=bool)
            est_xyz_cols = [f"dekf_t{tgt}_state_obs{obs}_{k}" for k in range(1, 4)]
            if all(c in df.columns for c in est_xyz_cols):
                est_pos[:, oi, ti, :] = df[est_xyz_cols].to_numpy(dtype=float)

            slot_col = find_target_slot_column(df, tgt, obs)
            if slot_col is not None:
                slot_series = numeric_col(df, slot_col).to_numpy(dtype=float)
                slots[:, oi, ti] = slot_series
                slot_ok = np.isfinite(slot_series) & (slot_series > 0.0)
                active[:, oi, ti] = est_ok & slot_ok
            else:
                active[:, oi, ti] = est_ok

    # Active-track start/stop events (used for HUD/event markers).
    active_start = np.zeros_like(active, dtype=bool)
    active_stop = np.zeros_like(active, dtype=bool)
    active_start[0, :, :] = active[0, :, :]
    active_start[1:, :, :] = active[1:, :, :] & (~active[:-1, :, :])
    active_stop[1:, :, :] = (~active[1:, :, :]) & active[:-1, :, :]

    # Visibility in each observer FoV and its entry/exit events.
    vis_in_fov = np.isfinite(dist_ot) & (dist_ot <= visibility_m)
    vis_entry = np.zeros_like(vis_in_fov, dtype=bool)
    vis_exit = np.zeros_like(vis_in_fov, dtype=bool)
    vis_entry[1:, :, :] = vis_in_fov[1:, :, :] & (~vis_in_fov[:-1, :, :])
    vis_exit[1:, :, :] = (~vis_in_fov[1:, :, :]) & vis_in_fov[:-1, :, :]

    # Frame decimation
    frame_ids = np.arange(0, n, max(1, args.stride), dtype=int)
    if len(frame_ids) == 0 or frame_ids[-1] != n - 1:
        frame_ids = np.append(frame_ids, n - 1)

    # Build interpolated animation states for smoother playback.
    substeps = max(1, int(args.interp_substeps))
    frame_states: list[tuple[int, int, float, float, int]] = []
    if len(frame_ids) == 1:
        k = int(frame_ids[0])
        frame_states.append((k, k, 0.0, float(t[k]), k))
    else:
        for i in range(len(frame_ids) - 1):
            k0 = int(frame_ids[i])
            k1 = int(frame_ids[i + 1])
            t0 = float(t[k0])
            t1 = float(t[k1])
            for s in range(substeps):
                alpha = float(s) / float(substeps)
                t_interp = (1.0 - alpha) * t0 + alpha * t1
                k_disc = k0 if alpha < 0.5 else k1
                frame_states.append((k0, k1, alpha, t_interp, k_disc))
        k = int(frame_ids[-1])
        frame_states.append((k, k, 0.0, float(t[k]), k))

    # View limits from all available points
    all_pts = np.vstack([
        obs_pos.reshape(-1, 3),
        tgt_pos.reshape(-1, 3),
        est_pos.reshape(-1, 3),
    ])
    valid = np.all(np.isfinite(all_pts), axis=1)
    if not np.any(valid):
        raise RuntimeError("No finite positions found for observers/targets.")
    pts = all_pts[valid]
    center = np.nanmean(pts, axis=0)
    radius = np.nanmax(np.linalg.norm(pts - center, axis=1)) * 1.05
    radius = max(radius, 1.0)

    # Parse events from log
    log_path = resolve_log_path(csv_path, args.log)
    events_by_time = parse_event_log(log_path)

    # Build nearest time lookup for events -> frame index
    event_by_frame = defaultdict(list)
    if events_by_time:
        for te, msgs in events_by_time.items():
            idx = int(np.argmin(np.abs(t - te)))
            event_by_frame[idx].extend(msgs)

    fig = plt.figure(figsize=(13.5, 9.0), facecolor="#f8fafc")
    ax = fig.add_subplot(111, projection="3d")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=0.04)
    ax.set_facecolor("#fbfcfe")
    ax.grid(False)
    ax.view_init(elev=24, azim=38)

    obs0 = obs_pos[0]
    obs0_valid = obs0[np.all(np.isfinite(obs0), axis=1)]
    if len(obs0_valid) > 0:
        center0 = np.mean(obs0_valid, axis=0)
        spread0 = np.max(np.linalg.norm(obs0_valid - center0, axis=1)) if len(obs0_valid) > 1 else 0.0
        radius0 = max(spread0 + visibility_m * args.focus_scale, visibility_m * 1.05)
    else:
        center0 = center
        radius0 = max(radius, visibility_m * 1.05)

    ax.set_xlim(center0[0] - radius0, center0[0] + radius0)
    ax.set_ylim(center0[1] - radius0, center0[1] + radius0)
    ax.set_zlim(center0[2] - radius0, center0[2] + radius0)
    ax.set_box_aspect((1.0, 1.0, 1.0))
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        axis.pane.set_edgecolor((0.85, 0.88, 0.92, 0.35))

    if args.earth:
        earth_r = 6_378_136.6
        u = np.linspace(0, 2 * np.pi, 40)
        v = np.linspace(0, np.pi, 20)
        x = earth_r * np.outer(np.cos(u), np.sin(v))
        y = earth_r * np.outer(np.sin(u), np.sin(v))
        z = earth_r * np.outer(np.ones_like(u), np.cos(v))
        ax.plot_surface(x, y, z, color="lightsteelblue", alpha=0.12, linewidth=0)

    # Pre-create artists
    # High-contrast palette for observers (targets stay a single shared color).
    obs_palette = [
        "#1f77b4",  # blue
        "#d62728",  # red
        "#9467bd",  # purple
        "#8c564b",  # brown
        "#e377c2",  # pink
        "#7f7f7f",  # gray
        "#17becf",  # cyan
        "#393b79",  # deep blue
        "#ad494a",  # brick red
    ]
    obs_colors = [obs_palette[i % len(obs_palette)] for i in range(len(observers))]
    obs_markers = []
    for oi, c_obs in enumerate(obs_colors):
        label = "Observers" if oi == 0 else "_nolegend_"
        (mk,) = ax.plot(
            [],
            [],
            [],
            linestyle="None",
            marker="o",
            markersize=11.5,
            markerfacecolor=c_obs,
            markeredgecolor="#0b1220",
            markeredgewidth=0.8,
            alpha=1.0,
            label=label,
        )
        obs_markers.append(mk)
    tgt_scatter = ax.scatter([], [], [], s=18, c="tab:orange", depthshade=False, alpha=0.62, label="Targets")
    est_scatter = ax.scatter([], [], [], s=11, marker="o", c=[], depthshade=False, edgecolors="none", alpha=0.72, label="_nolegend_")

    obs_trails = []
    for _ in observers:
        (ln,) = ax.plot([], [], [], color="#2563eb", alpha=0.58, linewidth=1.15, linestyle="-")
        obs_trails.append(ln)

    tgt_trails = []
    for _ in targets:
        (ln,) = ax.plot([], [], [], color="tab:orange", alpha=0.12, linewidth=0.65, linestyle="--")
        tgt_trails.append(ln)
    tgt_trails_active = []
    for _ in targets:
        (ln,) = ax.plot([], [], [], color="tab:orange", alpha=0.32, linewidth=1.0, linestyle="-")
        tgt_trails_active.append(ln)

    max_links = len(observers) * len(targets)
    link_lines = []
    for _ in range(max_links):
        (ln,) = ax.plot([], [], [], color="tab:green", alpha=0.24, linewidth=0.6, linestyle=":")
        link_lines.append(ln)
    est_trails = []
    for _ in range(max_links):
        (ln,) = ax.plot([], [], [], color="tab:red", alpha=0.42, linewidth=0.85, linestyle="-")
        est_trails.append(ln)
    lock_ring_scatter = ax.scatter(
        [],
        [],
        [],
        s=95,
        marker="o",
        facecolors="none",
        edgecolors="#ef4444",
        linewidths=1.4,
        alpha=0.92,
        depthshade=False,
        label="Locked targets",
    )

    # Event marker: only global track close (target no longer tracked by any observer).
    stop_scatter_global = ax.scatter([], [], [], s=62, marker="x", c="#111111", alpha=0.95, depthshade=False, label="Global close")

    # Visibility spheres (transparent surfaces), updated every frame around each observer.
    sphere_u = np.linspace(0.0, 2.0 * np.pi, 48)
    sphere_v = np.linspace(0.0, np.pi, 24)
    sphere_x = np.outer(np.cos(sphere_u), np.sin(sphere_v))
    sphere_y = np.outer(np.sin(sphere_u), np.sin(sphere_v))
    sphere_z = np.outer(np.ones_like(sphere_u), np.cos(sphere_v))
    visibility_spheres = [None for _ in observers]

    title = ax.set_title("Distributed Tracking Animation")
    ax.legend(loc="upper right")

    # HUD text
    hud = fig.text(
        0.02,
        0.96,
        "",
        va="top",
        ha="left",
        fontsize=10,
        family="monospace",
        bbox=dict(facecolor="white", alpha=0.86, edgecolor="lightgray"),
    )

    def set_scatter_xyz(sc, xyz):
        if xyz.size == 0:
            sc._offsets3d = ([], [], [])
            return
        sc._offsets3d = (xyz[:, 0], xyz[:, 1], xyz[:, 2])

    cam_center = np.array(center0, dtype=float)
    cam_radius = float(radius0)
    camera_smooth = min(max(float(args.camera_smooth), 0.0), 1.0)
    focus_boost = 0
    stop_ttl = np.zeros(len(targets), dtype=int)

    def update(frame_k: int):
        nonlocal cam_center, cam_radius, focus_boost, stop_ttl
        k0, k1, alpha, t_now, k_disc = frame_states[frame_k]

        op = lerp_xyz(obs_pos[k0], obs_pos[k1], alpha)
        tp = lerp_xyz(tgt_pos[k0], tgt_pos[k1], alpha)
        ep_all = lerp_xyz(est_pos[k0], est_pos[k1], alpha)
        active_now = active[k_disc]
        slot_now = slots[k_disc]

        # Recompute near-target mask on interpolated geometry.
        d_interp = np.linalg.norm(op[:, None, :] - tp[None, :, :], axis=2)
        near_t = np.any(np.isfinite(d_interp) & (d_interp <= visibility_m * args.near_margin), axis=0)

        # Per-frame event counts.
        init_count = int(np.count_nonzero(active_start[k_disc]))
        close_count = int(np.count_nonzero(active_stop[k_disc]))
        if init_count > 0:
            focus_boost = max(focus_boost, int(args.event_focus_frames))

        # Update close-event TTL markers on targets.
        stopped_tgts = np.where(np.any(active_stop[k_disc], axis=0))[0].astype(np.intp, copy=False)
        stop_ttl[stopped_tgts] = 8
        stop_ttl[stop_ttl > 0] -= 1

        # Refresh observer visibility spheres.
        if not args.hide_visibility_spheres:
            for oi in range(len(observers)):
                old_surface = visibility_spheres[oi]
                if old_surface is not None:
                    old_surface.remove()
                    visibility_spheres[oi] = None
                c = op[oi]
                if np.all(np.isfinite(c)):
                    pulse = int(np.count_nonzero(vis_entry[k_disc, oi, :]) + np.count_nonzero(vis_exit[k_disc, oi, :]))
                    pulse_alpha = 0.015 + min(0.04, 0.01 * pulse)
                    rgba = mcolors.to_rgba(obs_colors[oi], alpha=pulse_alpha)
                    sx = c[0] + visibility_m * sphere_x
                    sy = c[1] + visibility_m * sphere_y
                    sz = c[2] + visibility_m * sphere_z
                    visibility_spheres[oi] = ax.plot_surface(
                        sx,
                        sy,
                        sz,
                        color=rgba,
                        linewidth=0.0,
                        antialiased=True,
                        shade=False,
                    )

        # Camera follow with smoothing + temporary event zoom.
        op_valid = op[np.all(np.isfinite(op), axis=1)]
        if len(op_valid) > 0:
            c = np.mean(op_valid, axis=0)
            spread = np.max(np.linalg.norm(op_valid - c, axis=1)) if len(op_valid) > 1 else 0.0
            rr = max(spread + visibility_m * args.focus_scale, visibility_m * 1.05)
            if focus_boost > 0:
                zoom_gain = 0.25 * (focus_boost / max(1, int(args.event_focus_frames)))
                rr *= (1.0 - zoom_gain)
                focus_boost -= 1
            cam_center = (1.0 - camera_smooth) * cam_center + camera_smooth * c
            cam_radius = (1.0 - camera_smooth) * cam_radius + camera_smooth * rr
            ax.set_xlim(cam_center[0] - cam_radius, cam_center[0] + cam_radius)
            ax.set_ylim(cam_center[1] - cam_radius, cam_center[1] + cam_radius)
            ax.set_zlim(cam_center[2] - cam_radius, cam_center[2] + cam_radius)

        for oi, mk in enumerate(obs_markers):
            p = op[oi]
            if np.all(np.isfinite(p)):
                mk.set_data([p[0]], [p[1]])
                mk.set_3d_properties([p[2]])
            else:
                mk.set_data([], [])
                mk.set_3d_properties([])
        tp_show = tp[near_t] if len(tp) == len(near_t) else tp
        set_scatter_xyz(tgt_scatter, tp_show)

        # trails
        i0 = max(0, k_disc - args.trail)
        for oi, ln in enumerate(obs_trails):
            seg = obs_pos[i0:k_disc + 1, oi, :]
            ln.set_data(seg[:, 0], seg[:, 1])
            ln.set_3d_properties(seg[:, 2])

        show_target_trails = args.show_target_trails and (not args.hide_target_trails)
        if show_target_trails:
            for ti, ln in enumerate(tgt_trails):
                if not near_t[ti]:
                    ln.set_data([], [])
                    ln.set_3d_properties([])
                    tgt_trails_active[ti].set_data([], [])
                    tgt_trails_active[ti].set_3d_properties([])
                    continue
                seg = tgt_pos[i0:k_disc + 1, ti, :]
                seg_near = near_mask[i0:k_disc + 1, ti]
                seg = np.where(seg_near[:, None], seg, np.nan)
                ln.set_data(seg[:, 0], seg[:, 1])
                ln.set_3d_properties(seg[:, 2])
                if np.any(active_now[:, ti]):
                    tgt_trails_active[ti].set_data(seg[:, 0], seg[:, 1])
                    tgt_trails_active[ti].set_3d_properties(seg[:, 2])
                else:
                    tgt_trails_active[ti].set_data([], [])
                    tgt_trails_active[ti].set_3d_properties([])
        else:
            for ln in tgt_trails:
                ln.set_data([], [])
                ln.set_3d_properties([])
            for ln in tgt_trails_active:
                ln.set_data([], [])
                ln.set_3d_properties([])

        # links + estimated trajectories
        used = 0
        used_est = 0
        active_count = 0
        per_obs = []
        est_points = []
        est_point_colors = []
        lock_points = []
        active_err = []
        for obs in observers:
            oi = obs_idx[obs]
            links_for_obs = []
            candidates = []
            for tgt in targets:
                ti = tgt_idx[tgt]
                if active_now[oi, ti] and near_t[ti]:
                    p1 = op[oi, :]
                    p2 = tp[ti, :]
                    d = float(np.linalg.norm(p1 - p2)) if np.all(np.isfinite(p1)) and np.all(np.isfinite(p2)) else float("inf")
                    candidates.append((d, tgt))

            # By default keep all links; optional cap via CLI.
            candidates.sort(key=lambda x: x[0])
            if args.max_links_per_observer > 0:
                selected_targets = [tgt for _, tgt in candidates[: args.max_links_per_observer]]
            else:
                selected_targets = [tgt for _, tgt in candidates]

            for tgt in selected_targets:
                ti = tgt_idx[tgt]
                if active_now[oi, ti] and near_t[ti]:
                    active_count += 1
                    if used < max_links:
                        p1 = op[oi, :]
                        p2 = tp[ti, :]
                        link_lines[used].set_data([p1[0], p2[0]], [p1[1], p2[1]])
                        link_lines[used].set_3d_properties([p1[2], p2[2]])
                        link_lines[used].set_color(obs_colors[oi])
                        ep = ep_all[oi, ti, :]
                        if np.all(np.isfinite(ep)) and np.all(np.isfinite(p2)):
                            err = float(np.linalg.norm(ep - p2))
                            quality = float(math.exp(-err / 180.0))
                            link_lines[used].set_alpha(0.08 + 0.32 * quality)
                            link_lines[used].set_linewidth(0.35 + 0.50 * quality)
                        else:
                            link_lines[used].set_alpha(0.18)
                            link_lines[used].set_linewidth(0.45)
                        used += 1

                    ep = ep_all[oi, ti, :]
                    if used_est < max_links and np.all(np.isfinite(ep)):
                        # Estimated trajectory (observer-colored, continuous) over dashed target trajectory.
                        if not args.hide_est_trails:
                            seg_est = est_pos[i0:k_disc + 1, oi, ti, :]
                            seg_near = near_mask[i0:k_disc + 1, ti]
                            good = np.all(np.isfinite(seg_est), axis=1) & seg_near
                            seg_est = seg_est[good]
                            if len(seg_est) > 0:
                                est_trails[used_est].set_data(seg_est[:, 0], seg_est[:, 1])
                                est_trails[used_est].set_3d_properties(seg_est[:, 2])
                                est_trails[used_est].set_color(obs_colors[oi])
                                used_est += 1

                        est_points.append(ep)
                        est_point_colors.append(obs_colors[oi])

                        p2 = tp[ti, :]
                        if np.all(np.isfinite(p2)):
                            active_err.append(float(np.linalg.norm(ep - p2)))

                    slot = slot_now[oi, ti]
                    if np.isfinite(slot):
                        links_for_obs.append(f"t{tgt}(s{int(round(slot))})")
                    else:
                        links_for_obs.append(f"t{tgt}")

            if links_for_obs:
                per_obs.append(f"obs{obs}: {len(links_for_obs)}")

        for j in range(used, max_links):
            link_lines[j].set_data([], [])
            link_lines[j].set_3d_properties([])
        for j in range(used_est, max_links):
            est_trails[j].set_data([], [])
            est_trails[j].set_3d_properties([])

        est_arr = np.asarray(est_points, dtype=float) if est_points else np.empty((0, 3), dtype=float)
        set_scatter_xyz(est_scatter, est_arr)
        est_scatter.set_color(est_point_colors if est_point_colors else [])
        # One red ring around each target currently locked by at least one observer.
        for ti in range(len(targets)):
            if near_t[ti] and np.any(active_now[:, ti]) and np.all(np.isfinite(tp[ti])):
                lock_points.append(tp[ti])
        lock_arr = np.asarray(lock_points, dtype=float) if lock_points else np.empty((0, 3), dtype=float)
        set_scatter_xyz(lock_ring_scatter, lock_arr)

        mean_err = float(np.mean(active_err)) if active_err else float("nan")
        nearby_count = int(np.count_nonzero(near_t))
        if np.isfinite(mean_err):
            title.set_text(
                f"Distributed Tracking | t = {t_now:.1f} s | nearby targets = {nearby_count} | active links = {active_count} | mean |est-true| = {mean_err:.1f} m"
            )
        else:
            title.set_text(f"Distributed Tracking | t = {t_now:.1f} s | nearby targets = {nearby_count} | active links = {active_count}")

        ev = event_by_frame.get(k_disc, [])
        ev_txt = " | ".join(ev[:3]) if ev else "-"
        p_ratio = frame_k / max(1, len(frame_states) - 1)
        p_len = 24
        p_fill = int(round(p_ratio * p_len))
        p_bar = "█" * p_fill + "·" * (p_len - p_fill)
        hud.set_text(
            f"[{p_bar}] {100.0 * p_ratio:5.1f}%\n"
            f"t = {t_now:8.1f} s\n"
            f"active links: {active_count:4d} | nearby targets: {nearby_count:3d}\n"
            f"new starts: {init_count:3d} | closed: {close_count:3d}\n"
            f"mean |est-true|: {mean_err:7.2f} m\n"
            f"events: {ev_txt}"
        )

        # Event markers on targets (briefly shown): only global closes.
        stop_global_pts = []
        for ti in range(len(targets)):
            if stop_ttl[ti] > 0 and near_t[ti] and np.all(np.isfinite(tp[ti])):
                if not np.any(active_now[:, ti]):
                    stop_global_pts.append(tp[ti])
        stop_global_arr = np.asarray(stop_global_pts, dtype=float) if stop_global_pts else np.empty((0, 3), dtype=float)
        set_scatter_xyz(stop_scatter_global, stop_global_arr)

        artists = [tgt_scatter, est_scatter, lock_ring_scatter, stop_scatter_global, title, hud]
        artists.extend(obs_markers)
        artists.extend(obs_trails)
        artists.extend(tgt_trails)
        artists.extend(tgt_trails_active)
        artists.extend(link_lines)
        artists.extend(est_trails)
        return artists

    base_interval_ms = 1000.0 / max(1, args.fps)
    interval_ms = base_interval_ms / max(1e-6, args.speed_factor)
    ani = FuncAnimation(fig, update, frames=len(frame_states), interval=interval_ms, blit=False)

    if args.save:
        out = repo_path(args.save)
        out.parent.mkdir(parents=True, exist_ok=True)
        total_frames = len(frame_states)
        progress_step = max(1, total_frames // 10)

        def report_progress(frame_number, _total_frames):
            completed = frame_number + 1
            if completed == 1 or completed == total_frames or completed % progress_step == 0:
                percent = 100.0 * completed / max(1, total_frames)
                print(
                    f"Rendering animation: {completed}/{total_frames} "
                    f"frames ({percent:.0f}%)",
                    flush=True,
                )

        if out.suffix.lower() == ".gif":
            writer = "pillow"
            out_dpi = 100
        elif writers.is_available("ffmpeg"):
            writer = "ffmpeg"
            out_dpi = 170 if args.preset == "presentation" else 140
        else:
            out = out.with_suffix(".gif")
            writer = "pillow"
            out_dpi = 100
            print(
                "ffmpeg is not available; saving a GIF instead: "
                f"{out}",
                flush=True,
            )
        ani.save(
            out,
            writer=writer,
            fps=max(1, args.fps),
            dpi=out_dpi,
            progress_callback=report_progress,
        )
        print(f"Saved animation to: {out}")

    if args.show or not args.save:
        plt.show()


if __name__ == "__main__":
    main()
