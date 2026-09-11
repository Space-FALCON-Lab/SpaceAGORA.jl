"""
SectionF_plots.jl

Plotting functions for Section F conservation-law validation.

Reproduces the layout of Figures 2 and 3 in the paper:
  Figure 2 — Net OCL force, net OCL torque, and angular-momentum residual ε_H^(OCL)
              vs number of target orbits (panels: J2-inactive and J2-active)
  Figure 3 — ΔE_orb, W_OCL, and energy residual ε_E vs number of target orbits
              (panels: J2-inactive and J2-active)

Usage:
  plot_figure2(diag_no_j2, diag_j2; save_dir=".")
  plot_figure3(diag_no_j2, diag_j2; save_dir=".")

where each `diag_*` is the NamedTuple returned by `compute_sectionF_diagnostics`.
"""

using Plots, Printf

"""
    plot_figure2(diag_no_j2, diag_j2; save_dir=".")

Reproduce Figure 2 of the paper: net OCL force, net OCL torque, and
OCL-induced angular-momentum residual ε_H^(OCL) vs orbit count.

Left panel  (a): J2 inactive
Right panel (b): J2 active  (also shows H_z vs orbits as 4th subplot)
"""
function plot_figure2(diag_no_j2, diag_j2;
                      Δh_label = "Δh = 100 km",
                      save_dir  = @__DIR__,
                      filename  = "SectionF_Figure2_validation_force_torque_angmom.png")

    label = Δh_label
    color = :blue

    # ── Panel (a): J2 inactive ────────────────────────────────────────────────
    oc_a = diag_no_j2.orbit_count

    p1a = plot(oc_a, diag_no_j2.F_net_mag;
               label  = label,
               color  = color,
               lw     = 1.5,
               xlabel = "",
               ylabel = "F⁽ᴼᶜᴸ⁾_net  (N)",
               title  = "(a)  J₂ inactive",
               legend = :topright,
               titlefontsize  = 10,
               guidefontsize  = 9,
               tickfontsize   = 8,
               legendfontsize = 8,
               bottom_margin  = 0Plots.mm)

    p2a = plot(oc_a, diag_no_j2.tau_net_mag;
               label  = label,
               color  = color,
               lw     = 1.5,
               xlabel = "",
               ylabel = "τ⁽ᴼᶜᴸ⁾_net  (N·m)",
               legend = false,
               guidefontsize  = 9,
               tickfontsize   = 8,
               top_margin     = 0Plots.mm,
               bottom_margin  = 0Plots.mm)

    p3a = plot(oc_a, diag_no_j2.eps_H_mag;
               label  = label,
               color  = color,
               lw     = 1.5,
               xlabel = "Time  (number of target orbits)",
               ylabel = "ε⁽ᴼᶜᴸ⁾_H  (kg·m²/s)",
               legend = false,
               guidefontsize  = 9,
               tickfontsize   = 8,
               top_margin     = 0Plots.mm)

    panel_a = plot(p1a, p2a, p3a;
                   layout    = (3, 1),
                   size      = (420, 540),
                   left_margin = 10Plots.mm)

    # ── Panel (b): J2 active ──────────────────────────────────────────────────
    oc_b = diag_j2.orbit_count

    p1b = plot(oc_b, diag_j2.F_net_mag;
               label  = label,
               color  = color,
               lw     = 1.5,
               xlabel = "",
               ylabel = "F⁽ᴼᶜᴸ⁾_net  (N)",
               title  = "(b)  J₂ active",
               legend = :topright,
               titlefontsize  = 10,
               guidefontsize  = 9,
               tickfontsize   = 8,
               legendfontsize = 8,
               bottom_margin  = 0Plots.mm)

    p2b = plot(oc_b, diag_j2.tau_net_mag;
               label  = label,
               color  = color,
               lw     = 1.5,
               xlabel = "",
               ylabel = "τ⁽ᴼᶜᴸ⁾_net  (N·m)",
               legend = false,
               guidefontsize  = 9,
               tickfontsize   = 8,
               top_margin     = 0Plots.mm,
               bottom_margin  = 0Plots.mm)

    p3b = plot(oc_b, diag_j2.eps_H_mag;
               label  = label,
               color  = color,
               lw     = 1.5,
               xlabel = "",
               ylabel = "ε⁽ᴼᶜᴸ⁾_H  (kg·m²/s)",
               legend = false,
               guidefontsize  = 9,
               tickfontsize   = 8,
               top_margin     = 0Plots.mm,
               bottom_margin  = 0Plots.mm)

    # 4th subplot: H_z vs orbits (only in J2 panel)
    Hz0  = diag_j2.Hz_total[1]
    ΔHz  = diag_j2.Hz_total .- Hz0

    p4b = plot(oc_b, ΔHz;
               label  = label,
               color  = color,
               lw     = 1.5,
               xlabel = "Time  (number of target orbits)",
               ylabel = "ΔH_z  (kg·m²/s)",
               legend = false,
               guidefontsize  = 9,
               tickfontsize   = 8,
               top_margin     = 0Plots.mm)

    panel_b = plot(p1b, p2b, p3b, p4b;
                   layout      = (4, 1),
                   size        = (420, 720),
                   left_margin = 10Plots.mm)

    # ── Combine panels side-by-side ───────────────────────────────────────────
    fig = plot(panel_a, panel_b;
               layout = (1, 2),
               size   = (860, 720))

    save_path = joinpath(save_dir, filename)
    savefig(fig, save_path)
    println("Saved Figure 2 → ", save_path)
    return fig
end


"""
    plot_figure3(diag_no_j2, diag_j2; save_dir=".")

Reproduce Figure 3 of the paper: ΔE_orb, W_OCL, and energy residual ε_E vs orbit count.

Left panel  (a): J2 inactive
Right panel (b): J2 active
"""
function plot_figure3(diag_no_j2, diag_j2;
                      Δh_label = "Δh = 100 km",
                      save_dir  = @__DIR__,
                      filename  = "SectionF_Figure3_validation_energy.png")

    label = Δh_label
    color = :blue

    function _make_panel(diag, title_str, show_xlabel)
        oc = diag.orbit_count

        p1 = plot(oc, diag.delta_Eorb;
                  label  = label,
                  color  = color,
                  lw     = 1.5,
                  xlabel = "",
                  ylabel = "ΔE_orb  (J)",
                  title  = title_str,
                  legend = :topleft,
                  titlefontsize  = 10,
                  guidefontsize  = 9,
                  tickfontsize   = 8,
                  legendfontsize = 8,
                  bottom_margin  = 0Plots.mm)

        p2 = plot(oc, diag.W_OCL;
                  label  = label,
                  color  = color,
                  lw     = 1.5,
                  xlabel = "",
                  ylabel = "W_OCL  (J)",
                  legend = false,
                  guidefontsize  = 9,
                  tickfontsize   = 8,
                  top_margin     = 0Plots.mm,
                  bottom_margin  = 0Plots.mm)

        p3 = plot(oc, diag.eps_E;
                  label  = label,
                  color  = color,
                  lw     = 1.5,
                  xlabel = show_xlabel ? "Time  (number of target orbits)" : "",
                  ylabel = "ε_E  (J)",
                  legend = false,
                  guidefontsize  = 9,
                  tickfontsize   = 8,
                  top_margin     = 0Plots.mm)

        return plot(p1, p2, p3; layout=(3,1), size=(420,540), left_margin=10Plots.mm)
    end

    panel_a = _make_panel(diag_no_j2, "(a)  J₂ inactive", true)
    panel_b = _make_panel(diag_j2,    "(b)  J₂ active",   true)

    fig = plot(panel_a, panel_b; layout=(1,2), size=(860, 540))

    save_path = joinpath(save_dir, filename)
    savefig(fig, save_path)
    println("Saved Figure 3 → ", save_path)
    return fig
end


"""
    print_diagnostic_summary(diag; label="")

Print a concise table of max absolute values for each validation residual.
"""
function print_diagnostic_summary(diag; label::String = "")
    tag = isempty(label) ? "" : "[$label]  "
    @printf("\n%sValidation residual summary (max |value| over all time steps):\n", tag)
    @printf("  max |F_net^(OCL)|      = %.4e N\n",      maximum(diag.F_net_mag))
    @printf("  max |τ_net^(OCL)|      = %.4e N·m\n",    maximum(diag.tau_net_mag))
    @printf("  max |ε_E|              = %.4e J\n",       maximum(abs.(diag.eps_E)))
    @printf("  max |ε_H^(OCL)|        = %.4e kg·m²/s\n", maximum(diag.eps_H_mag))
    @printf("  ΔE_orb at t_f         = %.4e J\n",       diag.delta_Eorb[end])
    @printf("  W_OCL  at t_f         = %.4e J\n",       diag.W_OCL[end])
end
