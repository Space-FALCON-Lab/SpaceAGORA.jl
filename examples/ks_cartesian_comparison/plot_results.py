"""Render the native-coordinate comparison CSVs without requiring Julia plotting packages."""
from pathlib import Path
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullFormatter
import pandas as pd

OUT = Path(__file__).resolve().parents[2] / "outputs/ks_cartesian_comparison"
plt.rcParams.update({"font.size": 10, "pdf.fonttype": 42, "ps.fonttype": 42})
LABELS = {"a_m": r"$\Delta a$ (m)", "i_deg": r"$\Delta i$ (deg)",
          "omega_deg": "Delta argument of\nperiapsis (deg)"}


def tick(value, _):
    if value == 0 or .01 <= abs(value) <= 100:
        return f"{value:.5g}"
    exponent = math.floor(math.log10(abs(value)))
    return rf"${value / 10**exponent:.3g}\cdot 10^{{{exponent}}}$"


def save(fig, name):
    for ax in fig.axes:
        if ax.get_yscale() == "linear":
            ax.yaxis.set_major_formatter(FuncFormatter(tick))
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    data = pd.read_csv(OUT / "nonlinear.csv")
    truth = pd.read_csv(OUT / "truth_refinement.csv")
    long = OUT / "extended_horizons"
    if (long / "nonlinear.csv").exists():
        extended = pd.read_csv(long / "nonlinear.csv")
        extended_truth = pd.read_csv(long / "truth_refinement.csv")
        assert len(extended) == 240 and len(extended_truth) == 10, "Extended run is incomplete"
        assert extended.groupby(["e", "orbits"]).size().eq(24).all()
        data = extended
        truth = extended_truth
    data.to_csv(OUT / "nonlinear_all_horizons.csv", index=False)
    truth.to_csv(OUT / "truth_refinement_all_horizons.csv", index=False)
    horizons = sorted(data.orbits.unique())
    for e in data.e.unique():
        nrows = math.ceil(len(horizons) / 2)
        fig, grid = plt.subplots(nrows, 2, figsize=(7.1, 3.1 * nrows), squeeze=False)
        axes = grid.ravel()
        for ax, horizon in zip(axes, horizons):
            subset = data[(data.e == e) & (data.orbits == horizon)]
            for route, color in (("cartesian", "#176b9b"), ("ks_target_s", "#b33a3a")):
                for solver, marker in (("Tsit5", "o"), ("DP5", "s")):
                    for mode, fill in (("adaptive", color), ("fixed", "none")):
                        rows = subset[(subset.route == route) & (subset.integrator == solver) & (subset["mode"] == mode)]
                        ax.scatter(rows.total_s, rows.max_position_m, marker=marker,
                                   facecolors=fill, edgecolors=color, s=28,
                                   label=f"{'Cartesian' if route == 'cartesian' else 'KS (s)'} {solver} {mode}")
            ax.set(xscale="log", yscale="log", xlabel="Solve + sampling + metrics (s)",
                   ylabel="Maximum position error (m)", title=f"e = {e}; {horizon} {'orbit' if horizon == 1 else 'orbits'}")
            lo, hi = subset.total_s.min(), subset.total_s.max()
            ax.set_xticks([lo, math.sqrt(lo * hi), hi])
            ax.xaxis.set_major_formatter(FuncFormatter(tick))
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.grid(alpha=.2)
        if len(axes) > len(horizons):
            for ax in axes[len(horizons):]:
                ax.set_axis_off()
            handles, labels = axes[0].get_legend_handles_labels()
            axes[-1].legend(handles, labels, fontsize=8, loc="center")
        else:
            axes[-1].legend(fontsize=6, loc="best")
        save(fig, f"nonlinear_work_precision_e{e}")
    terminal = pd.read_csv(OUT / "terminal.csv")
    for objective in terminal.objective.unique():
        history = pd.read_csv(OUT / f"maximize_{objective}/history.csv")
        fig, axes = plt.subplots(3, 1, figsize=(7.1, 6.0), sharex=True)
        for ax, element in zip(axes, ("a_m", "i_deg", "omega_deg")):
            h = history[history.element == element]
            for key, label, style in (("ks_nonlinear", "Nonlinear KS", "-"),
                                      ("cartesian_nonlinear", "Nonlinear Cartesian", "--"),
                                      ("ks_linear", "KS-derived linear", ":"),
                                      ("cartesian_linear", "Direct Cartesian linear", "-.")):
                ax.plot(h.time_orbits, h[key], style, label=label, linewidth=1.2)
            ax.set_ylabel(LABELS[element])
            ax.grid(alpha=.2)
        axes[0].legend(fontsize=8)
        axes[-1].set_xlabel("Time (orbits)")
        save(fig, f"scheduling_{objective}")
        native = pd.read_csv(OUT / f"maximize_{objective}/inertial_zoh.csv")
        fig, axes = plt.subplots(3, 1, figsize=(7.1, 6.0), sharex=True)
        for ax, element in zip(axes, ("a_m", "i_deg", "omega_deg")):
            h = history[history.element == element]
            n = native[native.element == element]
            ax.plot(h.time_orbits, h.ks_nonlinear, label="Nonlinear KS", color="black")
            ax.plot(n.time_orbits, n.ks_native_s, "--", label="Native-s KS inertial ZOH")
            ax.plot(n.time_orbits, n.cartesian_inertial, ":", label="Cartesian inertial ZOH")
            ax.set_ylabel(LABELS[element])
            ax.grid(alpha=.2)
        axes[0].legend(fontsize=8)
        axes[-1].set_xlabel("Time (orbits)")
        save(fig, f"inertial_zoh_{objective}")
        refinement = pd.read_csv(OUT / f"maximize_{objective}/inertial_refinement.csv")
        fig, axes = plt.subplots(3, 1, figsize=(7.1, 6.0), sharex=True)
        for ax, element in zip(axes, ("a_m", "i_deg", "omega_deg")):
            r = refinement[refinement.element == element]
            ax.plot(r.substeps, r.ks_native_s, "o-", label="Native-s KS")
            ax.plot(r.substeps, r.cartesian_inertial, "s--", label="Cartesian inertial")
            h = history[history.element == element].iloc[-1]
            ax.axhline(h.ks_nonlinear, color="black", linestyle=":", label="Nonlinear KS")
            ax.set_ylabel(LABELS[element])
            ax.set_xscale("log", base=2)
            ax.set_xticks([1, 4, 16], ["1", "4", "16"])
            ax.grid(alpha=.2)
        axes[0].legend(fontsize=8)
        axes[-1].set_xlabel("Linear substeps per scheduling interval (schedule unchanged)")
        save(fig, f"inertial_refinement_{objective}")
        fig, axes = plt.subplots(2, 1, figsize=(7.1, 4.5), sharex=True)
        for ax, route in zip(axes, ("ks", "cartesian")):
            schedule = pd.read_csv(OUT / f"maximize_{objective}/schedule_{route}.csv")
            period = schedule.end_s.iloc[-1] / history.time_orbits.max()
            for column in schedule.columns[2:]:
                helper = int(column.split("_")[0])
                for row in schedule[schedule[column] == 1].itertuples():
                    ax.plot([row.start_s/period, row.end_s/period], [helper, helper], color="#176b9b", lw=2)
            ax.set(ylabel="Emitter ID", title=f"{objective}: {route}", ylim=(1, 28))
        axes[-1].set_xlabel("Time (orbits)")
        save(fig, f"schedules_{objective}")
    text = ["# Local native-coordinate comparison results", "",
            "Times are recorded measurements from the validation machine; see docs/ks_cartesian_comparison.md for methodology and limitations.", "",
            "## Fastest passing nonlinear settings tested", "",
            "Limits: 0.01 m maximum position error and 0.00001 m/s maximum velocity error, adding the measured truth-refinement difference as a safety margin (not a rigorous bound on reference error). This selects among the tested settings only; it is not a global matched-accuracy search. One warmed timing sample per configuration.", "",
            "| e | Orbits | Route | Solver | Mode | Setting | Solve (s) | Total (s) | Max position (m) |",
            "| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: |"]
    for (e, orbits, route), group in data.groupby(["e", "orbits", "route"]):
        uncertainty = truth[(truth.e == e) & (truth.orbits == orbits)].iloc[0]
        passing = group[(group.max_position_m + uncertainty.max_position_m <= .01) &
                        (group.max_velocity_m_s + uncertainty.max_velocity_m_s <= .00001)]
        if passing.empty:
            text.append(f"| {e} | {orbits} | {route} | No passing setting | | | | | |")
        else:
            best = passing.loc[passing.total_s.idxmin()]
            text.append(f"| {e} | {orbits} | {route} | {best.integrator} | {best['mode']} | {best.setting:g} | {best.solve_s:.4g} | {best.total_s:.4g} | {best.max_position_m:.4g} |")
    text += ["",
            "## Five-orbit scheduling results", "",
            "Both nonlinear routes replay the KS-selected schedule. Changes subtract each route's own no-laser trajectory.", "",
            "| Objective | Element | KS linear | Cartesian linear | KS nonlinear | Cartesian nonlinear | Different intervals |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: |"]
    for r in terminal.itertuples():
        text.append(f"| {r.objective} | {r.element} | {r.ks_linear:.10g} | {r.cartesian_linear:.10g} | {r.ks_nonlinear:.10g} | {r.cartesian_nonlinear:.10g} | {r.different_intervals} |")
    text += ["", "## Interpretation", "",
             "The scheduler uses the original frozen-interval physical-time model in RTN. Its KS route maps SpaceAGORA's analytical KS tangent to Cartesian before RTN/ZOH. It is not a native fictitious-time STM scheduler.",
             "The nonlinear benchmark is inertial: no RTN conversion, scheduling, or Jacobian assembly enters its solve timing. KS uses target fictitious time; fixed KS steps are nominal physical-second equivalents.",
             "Equal numerical tolerances on differently scaled states are not a matched-accuracy result. Use work-precision plots and truth_refinement.csv, not one runtime alone, to judge tradeoffs.",
             "A separate inertial_zoh.csv per objective propagates the selected schedule with native-s KS and inertial Cartesian matrices, with no RTN conversion. KS applies a fixed-physical-time clock correction after each interval. These finite-step approximations need not coincide with the RTN ZOH or with each other; they are not used to replace the original optimizer.",
             f"Nonlinear horizons completed: {', '.join(str(h) for h in horizons)} orbits; {len(data)} configurations. The scheduling and native inertial linear comparisons remain five-orbit runs. Historical affine/GVE/PPTX experiments have not been re-created by this driver."]
    (OUT / "README.md").write_text("\n".join(text) + "\n")


if __name__ == "__main__":
    main()
