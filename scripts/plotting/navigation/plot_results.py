#!/usr/bin/env python3
"""Compact plots for a saved navigation run."""

import argparse
import re
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SINGLE_OUTPUT_ROOT = (
    Path(tempfile.gettempdir())
    / "spaceagora_navigation"
    / "runs"
    / "single"
)
DEFAULT_PLOT_OUTPUT_ROOT = (
    Path(tempfile.gettempdir())
    / "spaceagora_navigation"
    / "plots"
)


def repo_path(value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else REPO_ROOT / path


def case_output_name(case: str) -> str:
    return case.strip().replace("-", "_")


def case_display_name(case: str) -> str:
    normalized = case_output_name(case)
    return "No DA" if normalized == "no_da" else normalized.replace("_", " ").title()


parser = argparse.ArgumentParser()
parser.add_argument(
    "--case",
    type=str,
    default="no-da",
    help="Navigation case used for the default single-run input"
)
parser.add_argument(
    "--input",
    type=str,
    default=None,
    help="Relative or absolute folder containing simulation_results.csv"
)
parser.add_argument(
    "--csv",
    type=str,
    default=None,
    help="Relative or absolute CSV path (overrides --input)"
)
parser.add_argument(
    "--output",
    type=str,
    default=None,
    help="Relative or absolute directory for saved plots"
)
parser.add_argument(
    "--max-target-plots",
    type=int,
    default=3,
    help="Maximum number of representative target plots (0 = all)."
)
args = parser.parse_args()

if args.csv is not None:
    CSV_PATH = repo_path(args.csv)
else:
    input_dir = (
        repo_path(args.input)
        if args.input
        else DEFAULT_SINGLE_OUTPUT_ROOT / case_output_name(args.case)
    )
    CSV_PATH = input_dir / "simulation_results.csv"

if not CSV_PATH.exists():
    raise FileNotFoundError(f"File not found: {CSV_PATH}")

CASE_LABEL = case_display_name(args.case)
PLOT_OUTPUT_DIR = (
    repo_path(args.output)
    if args.output
    else DEFAULT_PLOT_OUTPUT_ROOT / case_output_name(args.case)
)
PLOT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
SAVED_PLOT_PATHS = []

DF = pd.read_csv(CSV_PATH)
print(f'Loaded {CSV_PATH}')
print(f'Samples: {len(DF)}')
print(f'Saved plots directory: {PLOT_OUTPUT_DIR}')


def save_figure(fig, relative_path: str):
    output_path = PLOT_OUTPUT_DIR / relative_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    SAVED_PLOT_PATHS.append(output_path)
    print(f"Saved plot: {output_path}")


def print_view_command():
    if not SAVED_PLOT_PATHS:
        return
    quoted_paths = " ".join(f'"{path}"' for path in SAVED_PLOT_PATHS)
    print("\nCopy and paste to view the generated plots:")
    print(f"code {quoted_paths}")


def numeric_col(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.nan, index=df.index)
    s = df[col]
    if pd.api.types.is_bool_dtype(s):
        return s.astype(float)
    return pd.to_numeric(s, errors='coerce')


def discover_spacecraft_ids(df: pd.DataFrame):
    ids = set()
    rx = re.compile(r'^sc(\d+)_pos_1$')
    for c in df.columns:
        m = rx.match(c)
        if m:
            ids.add(int(m.group(1)))
    return sorted(ids)


def discover_target_ids(df: pd.DataFrame):
    ids = set()
    rx = re.compile(r'^dekf_t(\d+)_state_obs\d+_1$')
    for c in df.columns:
        m = rx.match(c)
        if m:
            ids.add(int(m.group(1)))
    return sorted(ids)


def discover_observer_ids(df: pd.DataFrame):
    ids = set()
    rx = re.compile(r'^dekf_t\d+_state_obs(\d+)_1$')
    for c in df.columns:
        m = rx.match(c)
        if m:
            ids.add(int(m.group(1)))
    return sorted(ids)


def find_target_slot_column(df: pd.DataFrame, target: int, obs: int):
    candidates = [
        f"dekf_t{target}_slot_obs{obs}",
        f"dekf_t{target}_slot_obs{obs}_1",
    ]
    for c in candidates:
        if c in df.columns:
            return c
    rx = re.compile(rf"^dekf_t{target}_slot_obs{obs}(?:_1)?$")
    for c in df.columns:
        if rx.match(c):
            return c
    return None


spacecraft = discover_spacecraft_ids(DF)
targets = discover_target_ids(DF)
observers = discover_observer_ids(DF)

if not spacecraft:
    raise RuntimeError('No spacecraft state columns found (expected sc{n}_pos_1).')
if not targets:
    raise RuntimeError('No DEKF target state columns found (expected dekf_t{target}_state_obs{obs}_{k}).')
if not observers:
    raise RuntimeError('No observer IDs found in DEKF columns.')

print(
    f'Scenario: {len(spacecraft)} spacecraft '
    f'({len(observers)} observers, {len(targets)} targets)'
)
print(f'Observer IDs: {observers}')

plot_df = DF.copy().reset_index(drop=True)
NAV_DT_SEC = 5.0
NAV_ALIGN_TOL_SEC = 0.05
t_all = pd.to_numeric(plot_df['time'], errors='coerce').to_numpy(dtype=float)
phase = np.mod(t_all, NAV_DT_SEC)
dist_to_grid = np.minimum(phase, NAV_DT_SEC - phase)
aligned_mask = np.isfinite(t_all) & (dist_to_grid <= NAV_ALIGN_TOL_SEC)
plot_df_est = plot_df.loc[aligned_mask].copy().reset_index(drop=True)
if len(plot_df_est) > 1:
    plot_df_est = plot_df_est.iloc[:-1].copy().reset_index(drop=True)
print(
    "Estimation plots: using nav-aligned samples and dropping final sample "
    f"(dt={NAV_DT_SEC:g}s, tol={NAV_ALIGN_TOL_SEC:g}s): {len(plot_df_est)}/{len(plot_df)} rows."
)
# -----------------------------------------------------------------------------
# Initial observer-only reference geometry
# -----------------------------------------------------------------------------
COMMUNICATION_RANGE_M = 300_000.0  # 300 km
VISIBILITY_RADIUS_M = 500_000.0  # 500 km
obs_color = 'tab:blue'

init_row = plot_df.iloc[0]


def kepler_orbit_from_state(position, velocity, samples=400):
    mu_earth = 3.986004418e14
    r0 = np.asarray(position, dtype=float)
    v0 = np.asarray(velocity, dtype=float)
    if not np.all(np.isfinite(r0)) or not np.all(np.isfinite(v0)):
        return None

    r_norm = np.linalg.norm(r0)
    angular_momentum = np.cross(r0, v0)
    h_norm = np.linalg.norm(angular_momentum)
    if r_norm <= 0.0 or h_norm <= 0.0:
        return None

    h_hat = angular_momentum / h_norm
    eccentricity_vector = np.cross(v0, angular_momentum) / mu_earth - r0 / r_norm
    eccentricity = np.linalg.norm(eccentricity_vector)
    if eccentricity >= 1.0:
        return None

    if eccentricity > 1.0e-10:
        periapsis_hat = eccentricity_vector / eccentricity
    else:
        periapsis_hat = r0 / r_norm
    transverse_hat = np.cross(h_hat, periapsis_hat)
    semilatus_rectum = h_norm**2 / mu_earth
    anomaly = np.linspace(0.0, 2.0 * np.pi, samples)
    radius = semilatus_rectum / (1.0 + eccentricity * np.cos(anomaly))
    return radius[:, None] * (
        np.cos(anomaly)[:, None] * periapsis_hat
        + np.sin(anomaly)[:, None] * transverse_hat
    )


def save_initial_target_scenario():
    earth_radius_km = 6378.137
    target_color = "0.45"
    observer_color = "navy"
    fig, ax = plt.subplots(figsize=(9.5, 9.0), constrained_layout=True)

    earth = plt.Circle(
        (0.0, 0.0),
        earth_radius_km,
        facecolor="#dff1fb",
        edgecolor="#7ec8f5",
        linewidth=1.0,
        alpha=0.65,
        zorder=0,
    )
    ax.add_patch(earth)

    target_positions = []
    observer_positions = []
    plotted_orbits = []
    for spacecraft_id in targets + observers:
        position_cols = [f"sc{spacecraft_id}_pos_{j}" for j in range(1, 4)]
        velocity_cols = [f"sc{spacecraft_id}_vel_{j}" for j in range(1, 4)]
        if not all(
            column in plot_df.columns
            for column in position_cols + velocity_cols
        ):
            continue

        position = np.array(
            [float(init_row[column]) for column in position_cols], dtype=float
        )
        velocity = np.array(
            [float(init_row[column]) for column in velocity_cols], dtype=float
        )
        orbit = kepler_orbit_from_state(position, velocity)
        if orbit is None:
            continue
        orbit_km = orbit / 1000.0
        plotted_orbits.append(orbit_km)

        if spacecraft_id in observers:
            ax.plot(
                orbit_km[:, 0],
                orbit_km[:, 1],
                color=observer_color,
                linewidth=1.5,
                alpha=0.8,
                zorder=3,
            )
            observer_positions.append(position / 1000.0)
        else:
            ax.plot(
                orbit_km[:, 0],
                orbit_km[:, 1],
                color=target_color,
                linewidth=0.55,
                alpha=0.12,
                zorder=1,
            )
            target_positions.append(position / 1000.0)

    if target_positions:
        target_positions_array = np.asarray(target_positions)
        ax.scatter(
            target_positions_array[:, 0],
            target_positions_array[:, 1],
            color=target_color,
            s=9,
            alpha=0.85,
            linewidths=0,
            zorder=2,
        )
    if observer_positions:
        observer_positions_array = np.asarray(observer_positions)
        ax.scatter(
            observer_positions_array[:, 0],
            observer_positions_array[:, 1],
            color=observer_color,
            s=35,
            edgecolors="white",
            linewidths=1.0,
            zorder=4,
        )

    if plotted_orbits:
        xy_values = np.concatenate([orbit[:, :2] for orbit in plotted_orbits])
        extent = max(
            earth_radius_km * 1.08,
            1.05 * float(np.nanmax(np.abs(xy_values))),
        )
    else:
        extent = earth_radius_km * 1.15

    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("ECI x [km]")
    ax.set_ylabel("ECI y [km]")
    ax.set_title("Target Population and Initial Orbital Geometry")
    ax.grid(True, alpha=0.18)
    ax.legend(
        handles=[
            Patch(
                facecolor="#dff1fb",
                edgecolor="#7ec8f5",
                alpha=0.65,
                label="Earth",
            ),
            Line2D(
                [0], [0], color=target_color, linewidth=1.0,
                alpha=0.35, label="Target reference orbits"
            ),
            Line2D(
                [0], [0], marker="o", linestyle="none", color=target_color,
                markersize=6, label="Target initial positions"
            ),
            Line2D(
                [0], [0], color=observer_color, linewidth=2.5,
                label="Observer reference orbits"
            ),
            Line2D(
                [0], [0], marker="o", linestyle="none", color=observer_color,
                markeredgecolor="white", markersize=7,
                label="Observer initial positions"
            ),
        ],
        loc="upper left",
    )
    save_figure(fig, "target_scenario.png")


save_initial_target_scenario()

obs_points = []
obs_ids_with_pos = []
for obs in observers:
    cols = [f'sc{obs}_pos_{j}' for j in range(1, 4)]
    if all(c in plot_df.columns for c in cols):
        obs_points.append([float(init_row[c]) for c in cols])
        obs_ids_with_pos.append(obs)
obs_points = np.asarray(obs_points, dtype=float) if obs_points else np.empty((0, 3), dtype=float)

if len(obs_points) > 0:
    fig_obs = plt.figure(figsize=(9.5, 8.0))
    ax_obs = fig_obs.add_subplot(111, projection='3d')
    obs_cmap = plt.get_cmap('tab20')

    for oi, obs in enumerate(obs_ids_with_pos):
        cols = [f'sc{obs}_pos_{j}' for j in range(1, 4)]
        traj = plot_df[cols].apply(pd.to_numeric, errors='coerce').to_numpy(dtype=float)
        finite = np.all(np.isfinite(traj), axis=1)
        color = obs_cmap(oi % obs_cmap.N)
        ax_obs.plot(
            np.where(finite, traj[:, 0], np.nan),
            np.where(finite, traj[:, 1], np.nan),
            np.where(finite, traj[:, 2], np.nan),
            color=color, linewidth=1.5, alpha=0.8
        )
        ax_obs.scatter(
            obs_points[oi, 0], obs_points[oi, 1], obs_points[oi, 2],
            color=color, s=58, edgecolors='black', linewidths=0.35,
            depthshade=False
        )
        ax_obs.text(
            obs_points[oi, 0], obs_points[oi, 1], obs_points[oi, 2],
            f'  obs{obs}', fontsize=8, color='0.20'
        )

    # Initial communication links: sparse local graph induced by R_com.
    for i in range(len(obs_points)):
        for j in range(i + 1, len(obs_points)):
            d = np.linalg.norm(obs_points[i] - obs_points[j])
            if np.isfinite(d) and d <= COMMUNICATION_RANGE_M:
                ax_obs.plot(
                    [obs_points[i, 0], obs_points[j, 0]],
                    [obs_points[i, 1], obs_points[j, 1]],
                    [obs_points[i, 2], obs_points[j, 2]],
                    color='0.20', linewidth=0.8, alpha=0.22
                )

    center = np.mean(obs_points, axis=0)
    formation_extent = np.max(np.linalg.norm(obs_points - center, axis=1))
    radius = max(float(formation_extent) * 1.35, COMMUNICATION_RANGE_M * 1.6)
    ax_obs.set_xlim(center[0] - radius, center[0] + radius)
    ax_obs.set_ylim(center[1] - radius, center[1] + radius)
    ax_obs.set_zlim(center[2] - radius, center[2] + radius)
    ax_obs.set_box_aspect((1, 1, 1))
    ax_obs.set_xlabel('X [m]')
    ax_obs.set_ylabel('Y [m]')
    ax_obs.set_zlabel('Z [m]')
    ax_obs.set_title('Observer Initial Configuration and Reference Orbits')
    ax_obs.legend(
        handles=[
            Line2D([0], [0], color=obs_color, linewidth=1.8, label='Observer reference orbits'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor=obs_color,
                   markeredgecolor='black', markersize=8, label='Initial observer positions'),
            Line2D([0], [0], color='0.20', linewidth=1.0, alpha=0.45,
                   label='Initial communication links (<= 300 km)')
        ],
        loc='upper left'
    )
    observer_plot_path = (
        PLOT_OUTPUT_DIR / 'observer_initial_configuration_reference_orbits.png'
    )
    fig_obs.savefig(observer_plot_path, dpi=180, bbox_inches='tight')
    SAVED_PLOT_PATHS.append(observer_plot_path)
    print(f'Saved observer-only reference plot: {observer_plot_path}')
    plt.close(fig_obs)
else:
    print('Observer-only reference plot skipped: observer position columns not found.')


def load_run_metrics(csv_path: Path):
    metrics_path = csv_path.parent / "association_quality_table.csv"
    if not metrics_path.exists():
        return {}, None
    table = pd.read_csv(metrics_path)
    required = {"section", "metric", "value"}
    if not required.issubset(table.columns):
        return {}, metrics_path
    metrics = {}
    for row in table.itertuples(index=False):
        try:
            metrics[f"{row.section}.{row.metric}"] = float(row.value)
        except (TypeError, ValueError):
            continue
    return metrics, metrics_path


def metric_text(metrics, key, unit="", digits=1):
    value = metrics.get(key, np.nan)
    if not np.isfinite(value):
        return "n/a"
    return f"{value:.{digits}f}{unit}"


def print_compact_tracking_summary():
    n_samples = len(plot_df_est)
    if n_samples == 0:
        print("Tracking summary skipped: no navigation-aligned samples.")
        return

    has_estimate_data = False
    active_target_count = np.zeros(n_samples, dtype=int)

    for target in targets:
        true_cols = [f"sc{target}_pos_{j}" for j in range(1, 4)]
        if not all(c in plot_df_est.columns for c in true_cols):
            continue
        true_position = plot_df_est[true_cols].to_numpy(dtype=float)
        target_active = np.zeros(n_samples, dtype=bool)

        for obs in observers:
            estimate_cols = [
                f"dekf_t{target}_state_obs{obs}_{j}" for j in range(1, 4)
            ]
            if not all(c in plot_df_est.columns for c in estimate_cols):
                continue
            estimate = plot_df_est[estimate_cols].to_numpy(dtype=float)
            valid = (
                np.all(np.isfinite(true_position), axis=1)
                & np.all(np.isfinite(estimate), axis=1)
            )

            has_estimate_data |= bool(np.any(valid))
            target_active |= valid

        active_target_count += target_active.astype(int)

    if not has_estimate_data:
        print("No valid target estimates found in the saved time interval.")

    time = pd.to_numeric(
        plot_df_est["time"], errors="coerce"
    ).to_numpy(dtype=float)
    metrics, metrics_path = load_run_metrics(CSV_PATH)
    print("\nRun summary")
    print(f"  duration: {time[-1] - time[0]:.0f} s")
    print(f"  targets with an estimate at final epoch: {active_target_count[-1]}/{len(targets)}")
    if metrics_path is None:
        print("  association metrics: n/a (association_quality_table.csv not found)")
        return

    print(f"  M2T precision / recall: "
          f"{metric_text(metrics, 'meas_assoc.commit_accuracy_pct', '%')} / "
          f"{metric_text(metrics, 'meas_assoc.recall_pct', '%')}")
    print(f"  T2T precision / recall: "
          f"{metric_text(metrics, 'track_assoc.tt_accuracy_pct_known_only', '%')} / "
          f"{metric_text(metrics, 'track_assoc.tt_recall_pct', '%')}")
    print(f"  tracking success: "
          f"{metric_text(metrics, 'tracking.success_rate_possible_pct', '%')}")
    print(f"  converged median error: "
          f"{metric_text(metrics, 'tracking.converged_median_error_m', ' m')}")


print_compact_tracking_summary()


def estimate_sample_count(target):
    count = 0
    for obs in observers:
        cols = [f"dekf_t{target}_state_obs{obs}_{j}" for j in range(1, 7)]
        if all(c in plot_df_est.columns for c in cols):
            estimate = plot_df_est[cols].to_numpy(dtype=float)
            count += int(np.count_nonzero(np.all(np.isfinite(estimate), axis=1)))
    return count


target_scores = {target: estimate_sample_count(target) for target in targets}
ranked_targets = sorted(
    targets, key=lambda target: (-target_scores[target], target)
)
if args.max_target_plots < 0:
    raise ValueError("--max-target-plots must be zero or positive")
representative_targets = (
    ranked_targets
    if args.max_target_plots == 0
    else ranked_targets[:args.max_target_plots]
)
print(f"Representative target plots: {representative_targets}")

# -----------------------------------------------------------------------------
# DEKF error norms (all observers in same plot, for each target)
# -----------------------------------------------------------------------------
t_est = pd.to_numeric(plot_df_est['time'], errors='coerce').to_numpy(dtype=float)
line_styles = ['-', '--', '-.', ':']
marker_styles = ['o', 's', '^', 'D', 'v', 'P', 'X', '*', '<', '>']
obs_to_color = {obs: f"C{i % 10}" for i, obs in enumerate(sorted(observers))}

for target in representative_targets:
    true_pos_cols = [f'sc{target}_pos_{j}' for j in range(1, 4)]
    true_vel_cols = [f'sc{target}_vel_{j}' for j in range(1, 4)]
    if not all(c in plot_df_est.columns for c in true_pos_cols + true_vel_cols):
        print(f'Skipping target SC{target}: missing true-state columns.')
        continue

    true_pos = plot_df_est[true_pos_cols].to_numpy(dtype=float)
    true_vel = plot_df_est[true_vel_cols].to_numpy(dtype=float)
    pos_norm_by_obs = {}
    vel_norm_by_obs = {}
    visibility_by_obs = {}
    slot_series_by_obs = {}

    for obs in observers:
        est_cols = [f'dekf_t{target}_state_obs{obs}_{j}' for j in range(1, 7)]
        if not all(c in plot_df_est.columns for c in est_cols):
            continue

        est = plot_df_est[est_cols].to_numpy(dtype=float)
        err_pos = est[:, :3] - true_pos
        err_vel = est[:, 3:] - true_vel

        pos_norm_by_obs[obs] = np.linalg.norm(err_pos, axis=1)
        vel_norm_by_obs[obs] = np.linalg.norm(err_vel, axis=1)

        slot_col = find_target_slot_column(plot_df_est, target, obs)
        if slot_col is not None:
            slot_series_by_obs[obs] = numeric_col(plot_df_est, slot_col).to_numpy(dtype=float)
        else:
            slot_series_by_obs[obs] = np.full(len(plot_df_est), np.nan, dtype=float)

        obs_pos_cols = [f'sc{obs}_pos_{j}' for j in range(1, 4)]
        if all(c in plot_df_est.columns for c in obs_pos_cols):
            obs_pos = plot_df_est[obs_pos_cols].to_numpy(dtype=float)
            dist = np.linalg.norm(true_pos - obs_pos, axis=1)
            visible = np.isfinite(dist) & (dist <= VISIBILITY_RADIUS_M)
            visibility_by_obs[obs] = visible
        else:
            visibility_by_obs[obs] = np.zeros(len(plot_df_est), dtype=bool)

    if not pos_norm_by_obs:
        print(f'Skipping target SC{target}: no DEKF observer states available.')
        continue

    fig, axes = plt.subplots(
        3, 1, figsize=(12, 8), sharex=True,
        gridspec_kw={"height_ratios": [3.0, 3.0, 1.4]}
    )
    ax_pos, ax_vel, ax_vis = axes

    for oi, obs in enumerate(sorted(pos_norm_by_obs)):
        color = obs_to_color[obs]
        linestyle = line_styles[oi % len(line_styles)]
        marker = marker_styles[oi % len(marker_styles)]
        markevery = (oi * 9, 45)

        y_pos = pos_norm_by_obs[obs]
        y_vel = vel_norm_by_obs[obs]
        y_pos_plot = y_pos.copy()
        y_vel_plot = y_vel.copy()

        # If saved slot mapping is available, break the line on slot switches to
        # avoid visual spikes from plotting two different tracks as one curve.
        slot_series = slot_series_by_obs.get(obs, np.full(len(plot_df_est), np.nan, dtype=float))
        if np.any(np.isfinite(slot_series)):
            slot_ids = np.where(np.isfinite(slot_series), np.rint(slot_series).astype(int), 0)
            switch = np.zeros(len(slot_ids), dtype=bool)
            switch[1:] = (
                (slot_ids[1:] > 0) & (slot_ids[:-1] > 0)
                & (slot_ids[1:] != slot_ids[:-1])
            )
            y_pos_plot[switch] = np.nan
            y_vel_plot[switch] = np.nan

        # White underlay keeps coincident curves visually separable.
        ax_pos.plot(t_est, y_pos_plot, color='white', linewidth=4.2, alpha=0.95, zorder=1)
        ax_vel.plot(t_est, y_vel_plot, color='white', linewidth=4.2, alpha=0.95, zorder=1)

        ax_pos.plot(
            t_est, y_pos_plot,
            color=color, linestyle=linestyle, linewidth=2.0,
            marker=marker, markersize=3.5, markevery=markevery,
            alpha=0.95, label=f'obs{obs}', zorder=2 + oi
        )
        ax_vel.plot(
            t_est, y_vel_plot,
            color=color, linestyle=linestyle, linewidth=2.0,
            marker=marker, markersize=3.5, markevery=markevery,
            alpha=0.95, label=f'obs{obs}', zorder=2 + oi
        )

        vis = visibility_by_obs.get(obs, np.zeros(len(plot_df_est), dtype=bool))
        y0 = oi
        ax_vis.fill_between(
            t_est, y0, y0 + 0.8, where=vis, step='post',
            color=color, alpha=0.35
        )
        ax_vis.step(
            t_est, y0 + 0.8 * vis.astype(float), where='post',
            color=color, linewidth=1.4
        )

    ax_pos.set_ylabel('||position error|| [m]')
    ax_pos.set_title(
        f'Target SC{target} - estimation errors + visibility ({CASE_LABEL})'
    )
    ax_pos.grid(True, alpha=0.35)
    ax_pos.legend(ncol=3)

    ax_vel.set_ylabel('||velocity error|| [m/s]')
    ax_vel.grid(True, alpha=0.35)

    ax_vis.set_ylabel('visibility')
    ax_vis.set_xlabel('Time [s]')
    ax_vis.set_yticks(np.arange(len(pos_norm_by_obs)) + 0.4)
    ax_vis.set_yticklabels([f'obs{obs}' for obs in sorted(pos_norm_by_obs)])
    ax_vis.set_ylim(0, max(1, len(pos_norm_by_obs)))
    ax_vis.grid(True, axis='x', alpha=0.25)
    ax_vis.text(
        0.01, 1.02, f'visible if distance <= {VISIBILITY_RADIUS_M/1000:.0f} km',
        transform=ax_vis.transAxes, fontsize=8, color='0.35'
    )

    plt.tight_layout()
    save_figure(
        fig,
        f"target_errors/target_sc{target}_estimation_errors.png"
    )


# -----------------------------------------------------------------------------
# RMSE per target over time (aggregated across all observers)
# -----------------------------------------------------------------------------
t_rmse = pd.to_numeric(plot_df_est["time"], errors="coerce").to_numpy(dtype=float)
rmse_ts_pos = {}
rmse_ts_vel = {}
rmse_rows = []

for target in targets:
    true_pos_cols = [f"sc{target}_pos_{j}" for j in range(1, 4)]
    true_vel_cols = [f"sc{target}_vel_{j}" for j in range(1, 4)]
    if not all(c in plot_df_est.columns for c in true_pos_cols + true_vel_cols):
        continue

    true_pos = plot_df_est[true_pos_cols].to_numpy(dtype=float)
    true_vel = plot_df_est[true_vel_cols].to_numpy(dtype=float)
    n = len(plot_df_est)
    sum_pos_sq = np.zeros(n, dtype=float)
    sum_vel_sq = np.zeros(n, dtype=float)
    count = np.zeros(n, dtype=float)

    for obs in observers:
        est_cols = [f"dekf_t{target}_state_obs{obs}_{j}" for j in range(1, 7)]
        if not all(c in plot_df_est.columns for c in est_cols):
            continue
        est = plot_df_est[est_cols].to_numpy(dtype=float)
        err_pos = est[:, :3] - true_pos
        err_vel = est[:, 3:] - true_vel
        valid = (
            np.isfinite(err_pos).all(axis=1)
            & np.isfinite(err_vel).all(axis=1)
            & np.isfinite(true_pos).all(axis=1)
            & np.isfinite(true_vel).all(axis=1)
        )
        if not np.any(valid):
            continue
        pos_sq = np.sum(err_pos ** 2, axis=1)
        vel_sq = np.sum(err_vel ** 2, axis=1)
        sum_pos_sq[valid] += pos_sq[valid]
        sum_vel_sq[valid] += vel_sq[valid]
        count[valid] += 1.0

    has = count > 0
    if not np.any(has):
        continue

    rmse_pos_t = np.full(n, np.nan, dtype=float)
    rmse_vel_t = np.full(n, np.nan, dtype=float)
    rmse_pos_t[has] = np.sqrt(sum_pos_sq[has] / count[has])
    rmse_vel_t[has] = np.sqrt(sum_vel_sq[has] / count[has])
    rmse_ts_pos[target] = rmse_pos_t
    rmse_ts_vel[target] = rmse_vel_t

    rmse_rows.append({
        "target": target,
        "overall_rmse_pos_m": float(np.sqrt(np.nanmean(sum_pos_sq[has] / count[has]))),
        "overall_rmse_vel_mps": float(np.sqrt(np.nanmean(sum_vel_sq[has] / count[has]))),
        "valid_time_samples": int(np.sum(has)),
    })

if rmse_rows:
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    for target in sorted(rmse_ts_pos):
        y = rmse_ts_pos[target]
        y_plot = np.where(y > 0.0, y, np.nan)
        axes[0].plot(t_rmse, y_plot, linewidth=1.8, label=f"T{target}")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("RMSE position [m] (log)")
    axes[0].set_title("RMSE over time by target (all observers pooled)")
    axes[0].grid(True, which="both", alpha=0.35)
    if len(rmse_ts_pos) <= 30:
        axes[0].legend(ncol=4, fontsize=8)

    for target in sorted(rmse_ts_vel):
        y = rmse_ts_vel[target]
        y_plot = np.where(y > 0.0, y, np.nan)
        axes[1].plot(t_rmse, y_plot, linewidth=1.8, label=f"T{target}")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Time [s]")
    axes[1].set_ylabel("RMSE velocity [m/s] (log)")
    axes[1].grid(True, which="both", alpha=0.35)
    plt.tight_layout()
    save_figure(fig, "rmse_all_tracks.png")
else:
    print("RMSE-over-time plots skipped: no valid target/observer estimate samples.")

print(
    f"\nPlots completed: {PLOT_OUTPUT_DIR}\n"
    "  target population and orbital geometry: target_scenario.png\n"
    "  initial orbital scenario: "
    "observer_initial_configuration_reference_orbits.png\n"
    "  representative target plots: target_errors/\n"
    "  all-track logarithmic errors: rmse_all_tracks.png"
)
print_view_command()
