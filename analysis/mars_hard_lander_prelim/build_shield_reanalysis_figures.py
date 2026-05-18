from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(
    "/Users/josephine/Research/Space-FALCON/_Proposals/[FY27]SURP/SHIELD_Reanalysis_Preliminary_Results"
)
PLOTS = ROOT / "plots"

NAVY = "#0b1f3a"
SKY = "#6eaef7"
TEAL = "#2da7b8"
ORANGE = "#f28c18"
MAIZE = "#ffcb05"
RED = "#c2410c"
GRAY = "#475569"
LIGHT = "#f8fafc"
PANEL = "#eaf2ff"


def _load_data():
    summary = pd.read_csv(ROOT / "prelim_summary.csv")
    authority = pd.read_csv(ROOT / "authority_summary.csv")
    target = pd.read_csv(ROOT / "target_guidance_summary.csv")
    return summary, authority, target


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


def _fmt_km(value: float) -> str:
    return f"{value:.0f}" if abs(value - round(value)) < 1e-9 else f"{value:.1f}"


def build_downrange_authority_plot(summary: pd.DataFrame, authority: pd.DataFrame, target: pd.DataFrame):
    bang = summary[summary["policy"] == "bang_bang"].sort_values("h_switch_km")
    body = summary[summary["policy"] == "body_only"].iloc[0]
    deployed = summary[summary["policy"] == "fixed_deployed"].iloc[0]
    auth = authority.iloc[0]
    target_row = target.iloc[0]
    switch_max_km = float(bang["h_switch_km"].max())

    fig, ax = plt.subplots(figsize=(11, 7))
    ax.set_facecolor(LIGHT)

    ax.plot(
        bang["h_switch_km"],
        bang["impact_downrange_km"],
        color=SKY,
        linewidth=3,
        marker="o",
        markersize=5,
        label="Deploy below h_switch",
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
        deployed["impact_downrange_km"],
        color=ORANGE,
        linestyle="--",
        linewidth=2.2,
        label=f"Always deployed ({deployed['impact_downrange_km']:.1f} km)",
        zorder=2,
    )
    ax.scatter(
        [target_row["h_switch_km"]],
        [target_row["achieved_range_km"]],
        color=RED,
        s=90,
        marker="*",
        label=f"Targeted solve ({target_row['h_switch_km']:.2f} km)",
        zorder=4,
    )

    x_bracket = bang["h_switch_km"].max() + 2.2
    y_lo = auth["min_impact_downrange_km"]
    y_hi = auth["max_impact_downrange_km"]
    ax.annotate(
        "",
        xy=(x_bracket, y_hi),
        xytext=(x_bracket, y_lo),
        arrowprops=dict(arrowstyle="<->", color=TEAL, lw=2.2),
        annotation_clip=False,
    )
    ax.text(
        x_bracket + 0.5,
        0.5 * (y_lo + y_hi),
        f"{auth['downrange_authority_km']:.1f} km\nreachable corridor",
        color=TEAL,
        fontsize=12,
        fontweight="bold",
        va="center",
    )

    ax.text(
        0.02,
        0.88,
        f"Published-size SHIELD surrogate: 120 kg entry mass, β_high = {auth['target_beta_high_kg_m2']:.1f} kg/m², β_low = {auth['target_beta_low_kg_m2']:.1f} kg/m²",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11,
        color=GRAY,
        bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="none", alpha=0.9),
    )
    ax.set_title("Downrange Authority From Subsonic Deployment Timing")
    ax.set_xlabel("Switch altitude h_switch (km)")
    ax.set_ylabel("Impact downrange (km)")
    ax.set_xlim(-0.3, switch_max_km + 3.2)
    ax.grid(True, linestyle=":", linewidth=0.8)
    ax.legend(loc="lower left", bbox_to_anchor=(0.02, 0.02))
    fig.tight_layout()
    fig.savefig(PLOTS / "downrange_authority_vs_switch_altitude.png", dpi=180)
    plt.close(fig)


def build_terminal_trade_plot(summary: pd.DataFrame):
    bang = summary[summary["policy"] == "bang_bang"].sort_values("h_switch_km")
    body = summary[summary["policy"] == "body_only"].iloc[0]
    deployed = summary[summary["policy"] == "fixed_deployed"].iloc[0]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax in axes:
        ax.set_facecolor(LIGHT)
        ax.grid(True, linestyle=":", linewidth=0.8)

    axes[0].plot(
        bang["h_switch_km"],
        bang["impact_velocity_mps"],
        color=SKY,
        linewidth=3,
        marker="o",
        markersize=5,
        label="Deploy below h_switch",
    )
    axes[0].axhline(body["impact_velocity_mps"], color=NAVY, linestyle="--", linewidth=2, label="Body only")
    axes[0].axhline(deployed["impact_velocity_mps"], color=ORANGE, linestyle="--", linewidth=2, label="Always deployed")
    axes[0].set_title("Impact Velocity vs Switch Altitude")
    axes[0].set_xlabel("Switch altitude h_switch (km)")
    axes[0].set_ylabel("Impact velocity (m/s)")
    axes[0].legend(loc="center right", bbox_to_anchor=(1.0, 0.84))

    axes[1].plot(
        bang["h_switch_km"],
        bang["peak_total_decel_earth_g"],
        color=TEAL,
        linewidth=3,
        marker="o",
        markersize=5,
        label="Deploy below h_switch",
    )
    axes[1].axhline(body["peak_total_decel_earth_g"], color=NAVY, linestyle="--", linewidth=2, label="Body only")
    axes[1].axhline(deployed["peak_total_decel_earth_g"], color=ORANGE, linestyle="--", linewidth=2, label="Always deployed")
    axes[1].set_title("Peak Aero Load vs Switch Altitude")
    axes[1].set_xlabel("Switch altitude h_switch (km)")
    axes[1].set_ylabel("Peak aero load (Earth-g)")
    axes[1].legend(loc="lower left")

    fig.suptitle("Impact Conditions Trade Against Deployment Timing", fontsize=16, color=NAVY, y=1.02)
    fig.tight_layout()
    fig.savefig(PLOTS / "terminal_velocity_and_peak_g_vs_switch_altitude.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def build_clean_drag_plot(summary: pd.DataFrame):
    body = pd.read_csv(ROOT / "trajectories" / "body_only_swnan_jetnan.csv")
    fixed = pd.read_csv(ROOT / "trajectories" / "fixed_deployed_swnan_jet000.csv")
    targeted = pd.read_csv(ROOT / "trajectories" / "targeted_nominal.csv")
    bang = summary[summary["policy"] == "bang_bang"].sort_values("h_switch_km")
    switch_max_km = float(bang["h_switch_km"].max())

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.set_facecolor(LIGHT)
    ax.grid(True, linestyle=":", linewidth=0.8)

    ax.plot(body["drag_accel_mps2"], body["altitude_km"], color=NAVY, linewidth=2.6, label="Body only")
    ax.plot(fixed["drag_accel_mps2"], fixed["altitude_km"], color=ORANGE, linewidth=2.6, label="Always deployed to impact")
    ax.plot(targeted["drag_accel_mps2"], targeted["altitude_km"], color=SKY, linewidth=3.0, label="Targeted switch solution")

    ax.axhline(switch_max_km, color=TEAL, linestyle="--", linewidth=1.8)
    ax.text(
        0.98,
        0.66,
        f"{_fmt_km(switch_max_km)} km\nsubsonic deployment\nband starts",
        transform=ax.transAxes,
        ha="right",
        va="top",
        color=TEAL,
        fontsize=11,
    )

    ax.set_title("Representative Drag Histories")
    ax.set_xlabel("Drag acceleration (m/s²)")
    ax.set_ylabel("Altitude (km)")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(PLOTS / "drag_accel_vs_altitude_clean.png", dpi=180)
    plt.close(fig)


def main():
    _style()
    PLOTS.mkdir(parents=True, exist_ok=True)
    summary, authority, target = _load_data()
    build_downrange_authority_plot(summary, authority, target)
    build_terminal_trade_plot(summary)
    build_clean_drag_plot(summary)
    print(PLOTS / "downrange_authority_vs_switch_altitude.png")
    print(PLOTS / "terminal_velocity_and_peak_g_vs_switch_altitude.png")
    print(PLOTS / "drag_accel_vs_altitude_clean.png")


if __name__ == "__main__":
    main()
