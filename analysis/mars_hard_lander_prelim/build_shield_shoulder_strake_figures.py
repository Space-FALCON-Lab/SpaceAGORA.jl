from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


TRAJ_ROOT = Path(
    "/Users/josephine/Research/Space-FALCON/_Proposals/[FY27]SURP/SHIELD_ShoulderStrake_Trajectory_Study"
)
STAB_ROOT = Path(
    "/Users/josephine/Research/Space-FALCON/_Proposals/[FY27]SURP/SHIELD_ShoulderStrake_Study"
)
PLOTS = TRAJ_ROOT / "plots"

NAVY = "#0b1f3a"
SKY = "#6eaef7"
TEAL = "#2da7b8"
ORANGE = "#f28c18"
MAIZE = "#ffcb05"
RED = "#c2410c"
GRAY = "#475569"
LIGHT = "#f8fafc"
GREEN = "#15803d"


def _style():
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.titlesize": 15,
            "axes.labelsize": 12,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "legend.frameon": False,
            "axes.edgecolor": "#cbd5e1",
            "axes.titlecolor": NAVY,
            "axes.labelcolor": GRAY,
            "xtick.color": GRAY,
            "ytick.color": GRAY,
            "grid.color": "#e2e8f0",
        }
    )


def _load():
    guided = pd.read_csv(TRAJ_ROOT / "shield_shoulder_strake_guided_summary.csv")
    pitch = pd.read_csv(TRAJ_ROOT / "shield_shoulder_strake_pitch_sweep.csv")
    yaw = pd.read_csv(TRAJ_ROOT / "shield_shoulder_strake_yaw_sweep.csv")
    stab = pd.read_csv(STAB_ROOT / "shield_shoulder_strake_stability_screen.csv")
    return guided, pitch, yaw, stab


def build_trim_capability_plot(stab: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(10.5, 6.8))
    ax.set_facecolor(LIGHT)
    ax.grid(True, linestyle=":", linewidth=0.8)

    nom_cg = 0.60
    if "strake_count" not in stab.columns:
        stab = stab.copy()
        stab["strake_count"] = 4
    subset = stab[(stab["cg_axial_fraction_of_body_length"].round(2) == nom_cg) & (stab["strake_count"] == 4)].copy()
    subset["area_pct"] = 100.0 * subset["strake_area_fraction_of_stowed_ref"]
    subset = subset.sort_values(["static_margin_fraction_of_diameter", "area_pct"])

    margin_styles = {
        0.05: (SKY, "5% static margin"),
        0.10: (TEAL, "10% static margin (nominal guided case)"),
        0.15: (ORANGE, "15% static margin"),
    }

    for margin, (color, label) in margin_styles.items():
        df = subset[subset["static_margin_fraction_of_diameter"].round(2) == margin]
        ax.plot(
            df["area_pct"],
            df["pitch_trim_cap_deg"],
            color=color,
            linewidth=3,
            marker="o",
            markersize=7,
            label=label,
            zorder=3,
        )

    ax.axhline(5.0, color=RED, linestyle="--", linewidth=2, zorder=2)
    ax.text(
        0.98,
        0.12,
        "Trim gate = 5 deg",
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=11,
        color=RED,
    )

    nominal = subset[
        (subset["static_margin_fraction_of_diameter"].round(2) == 0.10)
        & (subset["area_pct"].round(2) == 5.0)
    ].iloc[0]
    ax.scatter(
        [nominal["area_pct"]],
        [nominal["pitch_trim_cap_deg"]],
        color=GREEN,
        s=100,
        marker="*",
        zorder=4,
    )
    ax.annotate(
        "Nominal guided case\n5% area (~6.0 kg)\ntrim cap = 5.51 deg",
        xy=(nominal["area_pct"], nominal["pitch_trim_cap_deg"]),
        xytext=(5.45, 7.3),
        fontsize=11,
        color=GREEN,
        arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.8),
        bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="none", alpha=0.95),
    )

    ax.text(
        0.60,
        0.89,
        "Pitch/yaw trim caps are identical in this first-cut surrogate.",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=11,
        color=GRAY,
        bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="none", alpha=0.92),
    )

    eight = stab[
        (stab["cg_axial_fraction_of_body_length"].round(2) == nom_cg)
        & (stab["static_margin_fraction_of_diameter"].round(2) == 0.10)
        & (stab["strake_count"] == 8)
    ].copy()
    if not eight.empty:
        eight["area_pct"] = 100.0 * eight["strake_area_fraction_of_stowed_ref"]
        eight["each_area_pct"] = eight["area_pct"] / eight["strake_count"]
        nominal_four = subset[subset["static_margin_fraction_of_diameter"].round(2) == 0.10].copy()
        nominal_four["each_area_pct"] = nominal_four["area_pct"] / 4.0
        target_each = set(np.round(nominal_four["each_area_pct"], 4))
        eight_same = eight[np.round(eight["each_area_pct"], 4).isin(target_each)].sort_values("area_pct")
        if not eight_same.empty:
            ax.plot(
                eight_same["area_pct"],
                eight_same["pitch_trim_cap_deg"],
                color=NAVY,
                linewidth=2.0,
                linestyle="--",
                marker="D",
                markersize=6,
                label="8 strakes, same individual size (10% static margin)",
                zorder=3,
            )
            last = eight_same.iloc[-1]
            ax.annotate(
                "8 strakes, same individual size\nmeans more total area and mass",
                xy=(last["area_pct"], last["pitch_trim_cap_deg"]),
                xytext=(7.1, last["pitch_trim_cap_deg"] - 0.7),
                fontsize=10.5,
                color=NAVY,
                arrowprops=dict(arrowstyle="->", color=NAVY, lw=1.6),
                bbox=dict(boxstyle="round,pad=0.26", facecolor="white", edgecolor="none", alpha=0.94),
            )

    ax.set_title("Shoulder-Strake Trim Capability On Frozen SHIELD Baseline")
    ax.set_xlabel("Total shoulder-strake area (% of stowed reference area)")
    ax.set_ylabel("Maximum trimmed attitude offset (deg)")
    x_max = 10.5 if "strake_count" in stab.columns and (stab["strake_count"] == 8).any() else 5.25
    ax.set_xlim(0.75, x_max)
    ax.set_ylim(0.0, max(14.5, subset["pitch_trim_cap_deg"].max() + 0.8))
    ax.legend(loc="upper left", bbox_to_anchor=(0.01, 1.0))
    fig.tight_layout()
    fig.savefig(PLOTS / "shoulder_strake_trim_capability.png", dpi=180)
    plt.close(fig)


def build_downrange_authority_plot(guided: pd.DataFrame, pitch: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(10.8, 6.9))
    ax.set_facecolor(LIGHT)
    ax.grid(True, linestyle=":", linewidth=0.8)

    body = guided[guided["case"] == "body_only"].iloc[0]
    passive = guided[guided["case"] == "passive_skirt_subsonic"].iloc[0]
    targeted = guided[guided["case"] == "pitch_targeted"].iloc[0]
    min_row = guided[guided["case"] == "pitch_min_range"].iloc[0]
    max_row = guided[guided["case"] == "pitch_max_range"].iloc[0]
    authority = float(max_row["impact_downrange_km"] - min_row["impact_downrange_km"])
    pitch_x = pitch["pitch_deflection_deg"]
    max_abs_deflection = float(pitch_x.abs().max())

    ax.plot(
        pitch_x,
        pitch["impact_downrange_km"],
        color=SKY,
        linewidth=3,
        marker="o",
        markersize=5,
        label="Pitch-command sweep",
        zorder=3,
    )
    ax.axhline(
        body["impact_downrange_km"],
        color=NAVY,
        linestyle="--",
        linewidth=2.2,
        label=f"Body only ({body['impact_downrange_km']:.1f} km)",
        zorder=2,
    )
    ax.axhline(
        passive["impact_downrange_km"],
        color=ORANGE,
        linestyle="--",
        linewidth=2.2,
        label=f"Passive skirt at 6 km ({passive['impact_downrange_km']:.1f} km)",
        zorder=2,
    )
    ax.scatter(
        [15.0 * targeted["command_value"]],
        [targeted["impact_downrange_km"]],
        color=RED,
        s=110,
        marker="*",
        zorder=4,
        label=f"Targeted case ({targeted['impact_downrange_km']:.1f} km)",
    )

    x_bracket = 15.9
    ax.annotate(
        "",
        xy=(x_bracket, max_row["impact_downrange_km"]),
        xytext=(x_bracket, min_row["impact_downrange_km"]),
        arrowprops=dict(arrowstyle="<->", color=TEAL, lw=2.2),
        annotation_clip=False,
    )
    ax.text(
        x_bracket + 0.05,
        0.5 * (min_row["impact_downrange_km"] + max_row["impact_downrange_km"]),
        f"{authority:.1f} km\nreachable corridor",
        color=TEAL,
        fontsize=12,
        fontweight="bold",
        va="center",
    )

    ax.text(
        0.02,
        0.95,
        f"Pitch mode uses equal-and-opposite top/bottom strake deflections up to ±{max_abs_deflection:.0f} deg.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11,
        color=GRAY,
        bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="none", alpha=0.92),
    )
    ax.set_title("Downrange Authority From Pre-Skirt Shoulder-Strake Guidance")
    ax.set_xlabel("Pitch strake deflection δ (deg)")
    ax.set_ylabel("Impact downrange (km)")
    ax.set_xlim(-1.05 * max_abs_deflection, 1.23 * max_abs_deflection)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(PLOTS / "shoulder_strake_downrange_authority_vs_pitch_command.png", dpi=180)
    plt.close(fig)


def build_impact_trade_plot(guided: pd.DataFrame, pitch: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax in axes:
        ax.set_facecolor(LIGHT)
        ax.grid(True, linestyle=":", linewidth=0.8)

    body = guided[guided["case"] == "body_only"].iloc[0]
    passive = guided[guided["case"] == "passive_skirt_subsonic"].iloc[0]
    targeted = guided[guided["case"] == "pitch_targeted"].iloc[0]
    pitch_x = pitch["pitch_deflection_deg"]
    targeted_deflection = 15.0 * targeted["command_value"]

    axes[0].plot(
        pitch_x,
        pitch["impact_velocity_mps"],
        color=SKY,
        linewidth=3,
        marker="o",
        markersize=5,
        label="Pitch-command sweep",
    )
    axes[0].axhline(body["impact_velocity_mps"], color=NAVY, linestyle="--", linewidth=2.0, label="Body only")
    axes[0].axhline(passive["impact_velocity_mps"], color=ORANGE, linestyle="--", linewidth=2.0, label="Passive skirt at 6 km")
    axes[0].scatter([targeted_deflection], [targeted["impact_velocity_mps"]], color=RED, s=90, marker="*")
    axes[0].set_title("Impact Velocity vs Pitch-Strake Deflection")
    axes[0].set_xlabel("Pitch strake deflection δ (deg)")
    axes[0].set_ylabel("Impact velocity (m/s)")
    axes[0].legend(loc="center right", bbox_to_anchor=(0.98, 0.76))

    axes[1].plot(
        pitch_x,
        pitch["peak_total_decel_earth_g"],
        color=TEAL,
        linewidth=3,
        marker="o",
        markersize=5,
        label="Pitch-command sweep",
    )
    axes[1].axhline(body["peak_total_decel_earth_g"], color=NAVY, linestyle="--", linewidth=2.0, label="Body only")
    axes[1].axhline(passive["peak_total_decel_earth_g"], color=ORANGE, linestyle="--", linewidth=2.0, label="Passive skirt at 6 km")
    axes[1].scatter([targeted_deflection], [targeted["peak_total_decel_earth_g"]], color=RED, s=90, marker="*")
    axes[1].set_title("Peak Aero Load vs Pitch-Strake Deflection")
    axes[1].set_xlabel("Pitch strake deflection δ (deg)")
    axes[1].set_ylabel("Peak aero load (Earth-g)")
    axes[1].legend(loc="lower left")

    fig.suptitle("Impact Conditions Trade Against Shoulder-Strake Deflection", fontsize=16, color=NAVY, y=1.02)
    fig.tight_layout()
    fig.savefig(PLOTS / "shoulder_strake_impact_trade_vs_pitch_command.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def build_impact_footprint_plot(guided: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(10.2, 6.7))
    ax.set_facecolor(LIGHT)
    ax.grid(True, linestyle=":", linewidth=0.8)

    style = {
        "body_only": (NAVY, "o", "Body only"),
        "passive_skirt_subsonic": (ORANGE, "s", "Passive skirt at 6 km"),
        "pitch_min_range": (SKY, "v", "Pitch min-range"),
        "pitch_targeted": (RED, "*", "Pitch targeted"),
        "pitch_max_range": (TEAL, "^", "Pitch max-range"),
        "yaw_left": (GRAY, "<", "Yaw left"),
        "yaw_right": (GREEN, ">", "Yaw right"),
    }

    for _, row in guided.iterrows():
        color, marker, label = style[row["case"]]
        size = 170 if marker == "*" else 85
        ax.scatter(
            row["impact_downrange_km"],
            row["impact_crossrange_km"],
            color=color,
            marker=marker,
            s=size,
            zorder=4,
            label=label,
        )

    pitch_min = guided[guided["case"] == "pitch_min_range"].iloc[0]
    pitch_max = guided[guided["case"] == "pitch_max_range"].iloc[0]
    yaw_left = guided[guided["case"] == "yaw_left"].iloc[0]
    yaw_right = guided[guided["case"] == "yaw_right"].iloc[0]
    targeted = guided[guided["case"] == "pitch_targeted"].iloc[0]

    ax.plot(
        [pitch_min["impact_downrange_km"], pitch_max["impact_downrange_km"]],
        [pitch_min["impact_crossrange_km"], pitch_max["impact_crossrange_km"]],
        color=SKY,
        linewidth=2.0,
        alpha=0.85,
        zorder=2,
    )
    ax.plot(
        [yaw_left["impact_downrange_km"], yaw_right["impact_downrange_km"]],
        [yaw_left["impact_crossrange_km"], yaw_right["impact_crossrange_km"]],
        color=GREEN,
        linewidth=2.0,
        alpha=0.85,
        zorder=2,
    )

    ax.annotate(
        f"{pitch_max['impact_downrange_km'] - pitch_min['impact_downrange_km']:.1f} km downrange span",
        xy=(targeted["impact_downrange_km"], 0.6),
        xytext=(664.0, 2.4),
        fontsize=11,
        color=SKY,
        arrowprops=dict(arrowstyle="->", color=SKY, lw=1.8),
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="none", alpha=0.95),
    )
    ax.annotate(
        f"{yaw_right['impact_crossrange_km'] - yaw_left['impact_crossrange_km']:.1f} km crossrange span",
        xy=(671.25, 4.1),
        xytext=(678.2, 4.55),
        fontsize=11,
        color=GREEN,
        arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.8),
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="none", alpha=0.95),
    )

    ax.set_title("Nominal Deterministic Endpoint Map From Shoulder-Strake Commands")
    ax.set_xlabel("Impact downrange (km)")
    ax.set_ylabel("Impact crossrange (km)")
    ax.set_xlim(658.0, 685.5)
    ax.set_ylim(-5.3, 5.5)
    ax.legend(loc="lower right", ncol=2)
    fig.tight_layout()
    fig.savefig(PLOTS / "shoulder_strake_impact_footprint.png", dpi=180)
    plt.close(fig)


def _clean_terminal_artifact(df: pd.DataFrame) -> pd.DataFrame:
    if len(df) < 2:
        return df
    last = df.iloc[-1]
    prev = df.iloc[-2]
    deployed_col = "deployed" in df.columns
    state_col = "state_label" in df.columns
    if (
        abs(float(last["altitude_km"])) < 1e-9
        and deployed_col
        and state_col
        and bool(last["deployed"]) is False
        and bool(prev["deployed"]) is True
        and str(last["state_label"]) != str(prev["state_label"])
    ):
        return df.iloc[:-1].copy()
    return df


def _interp_against_altitude(df: pd.DataFrame, y_column: str, altitude_grid: np.ndarray) -> np.ndarray:
    clean = _clean_terminal_artifact(df)
    alt = clean["altitude_km"].to_numpy()
    y = clean[y_column].to_numpy()
    order = np.argsort(alt)
    return np.interp(altitude_grid, alt[order], y[order])


def _build_altitude_envelope_plot(
    *,
    y_column: str,
    x_label: str,
    title: str,
    outfile: str,
    xlim_left: float | None = None,
    xlim_right: float | None = None,
    add_mach1: bool = False,
    note: str = "Shaded envelope spans all tried pitch and yaw command cases.",
):
    traj_dir = TRAJ_ROOT / "trajectories"
    pitch_paths = sorted(traj_dir.glob("pitch_cmd_*.csv"))
    yaw_paths = sorted(traj_dir.glob("yaw_cmd_*.csv"))
    command_paths = pitch_paths + yaw_paths
    if not command_paths:
        return

    command_trajs = [_clean_terminal_artifact(pd.read_csv(path)) for path in command_paths]
    body = _clean_terminal_artifact(pd.read_csv(traj_dir / "body_only.csv"))
    passive = _clean_terminal_artifact(pd.read_csv(traj_dir / "passive_skirt_subsonic.csv"))
    targeted = _clean_terminal_artifact(pd.read_csv(traj_dir / "pitch_targeted.csv"))

    altitude_grid = np.linspace(0.0, 125.0, 500)

    command_stack = np.vstack([_interp_against_altitude(df, y_column, altitude_grid) for df in command_trajs])
    y_min = command_stack.min(axis=0)
    y_max = command_stack.max(axis=0)
    body_y = _interp_against_altitude(body, y_column, altitude_grid)
    passive_y = _interp_against_altitude(passive, y_column, altitude_grid)
    targeted_y = _interp_against_altitude(targeted, y_column, altitude_grid)

    fig, ax = plt.subplots(figsize=(10.5, 6.8))
    ax.set_facecolor(LIGHT)
    ax.grid(True, linestyle=":", linewidth=0.8)

    ax.fill_betweenx(
        altitude_grid,
        y_min,
        y_max,
        color=SKY,
        alpha=0.22,
        label="All pitch/yaw command cases",
    )
    ax.plot(body_y, altitude_grid, color=NAVY, linestyle="--", linewidth=2.4, label="Body only")
    ax.plot(passive_y, altitude_grid, color=ORANGE, linestyle="--", linewidth=2.4, label="Passive skirt at 6 km")
    ax.plot(targeted_y, altitude_grid, color=TEAL, linewidth=2.8, label="Targeted guided case")

    if add_mach1:
        ax.axvline(1.0, color=RED, linestyle=":", linewidth=2.0)
        ax.text(
            1.02,
            0.10,
            "Mach 1",
            transform=ax.get_xaxis_transform(),
            color=RED,
            fontsize=11,
            ha="left",
            va="bottom",
        )

    ax.text(
        0.02,
        0.95,
        note,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11,
        color=GRAY,
        bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="none", alpha=0.92),
    )

    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Altitude (km)")
    if xlim_left is not None or xlim_right is not None:
        ax.set_xlim(left=xlim_left, right=xlim_right)
    else:
        ax.set_xlim(left=0.0)
    ax.set_ylim(0.0, 125.0)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(PLOTS / outfile, dpi=180)
    plt.close(fig)


def build_mach_altitude_envelope_plot():
    _build_altitude_envelope_plot(
        y_column="mach",
        x_label="Mach number",
        title="Mach History Envelope Across Shoulder-Strake Pitch Cases",
        outfile="shoulder_strake_mach_vs_altitude_envelope.png",
        xlim_left=0.0,
        add_mach1=True,
        note="Shaded envelope spans all tried pitch and yaw command cases.\nFinal surface bookkeeping artifact removed from the plotted histories.",
    )


def build_velocity_altitude_envelope_plot():
    _build_altitude_envelope_plot(
        y_column="velocity_mps",
        x_label="Velocity (m/s)",
        title="Velocity History Envelope Across Shoulder-Strake Pitch Cases",
        outfile="shoulder_strake_velocity_vs_altitude_envelope.png",
        xlim_left=0.0,
    )


def build_dynamic_pressure_altitude_envelope_plot():
    _build_altitude_envelope_plot(
        y_column="q_pa",
        x_label="Dynamic pressure q (Pa)",
        title="Dynamic Pressure Envelope Across Shoulder-Strake Pitch Cases",
        outfile="shoulder_strake_dynamic_pressure_vs_altitude_envelope.png",
        xlim_left=0.0,
    )


def main():
    _style()
    PLOTS.mkdir(parents=True, exist_ok=True)
    guided, pitch, yaw, stab = _load()
    build_trim_capability_plot(stab)
    build_downrange_authority_plot(guided, pitch)
    build_impact_trade_plot(guided, pitch)
    build_impact_footprint_plot(guided)
    build_mach_altitude_envelope_plot()
    build_velocity_altitude_envelope_plot()
    build_dynamic_pressure_altitude_envelope_plot()
    print(PLOTS / "shoulder_strake_trim_capability.png")
    print(PLOTS / "shoulder_strake_downrange_authority_vs_pitch_command.png")
    print(PLOTS / "shoulder_strake_impact_trade_vs_pitch_command.png")
    print(PLOTS / "shoulder_strake_impact_footprint.png")
    print(PLOTS / "shoulder_strake_mach_vs_altitude_envelope.png")
    print(PLOTS / "shoulder_strake_velocity_vs_altitude_envelope.png")
    print(PLOTS / "shoulder_strake_dynamic_pressure_vs_altitude_envelope.png")


if __name__ == "__main__":
    main()
