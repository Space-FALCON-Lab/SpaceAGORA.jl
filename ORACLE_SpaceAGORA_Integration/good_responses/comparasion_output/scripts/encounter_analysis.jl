"""
Shared, code-agnostic post-processing/plotting functions for the ORACLE-SpaceAGORA
laser-link encounter comparison (30 orbits, gve_sma schedule).

Consumes plain Vectors of time / per-satellite ECI position & velocity (SVector{3,Float64},
meters & m/s) plus a per-sample "which helper is actively firing" index (0 = none), and
produces the following deliverables for one target-vs-helper pair:
  1) relative-speed magnitude |v_target - v_helper| vs time, one plot per range-feasible encounter
  2) signed range-rate vs time (with range + max-range band on a twin axis)
  3) geometric contact-window duration vs real laser-on duration, per encounter
plus, across all pairs of one code run:
  4) combined relative-speed vs time-since-encounter-start (all encounters overlaid)
  5) combined range & range-rate vs time-since-encounter-start (all encounters overlaid)
"""

using Plots, Printf, LinearAlgebra
gr()
default(titlefontsize=10)

const Encounter = NamedTuple{(:i0, :i1, :t0, :t1), Tuple{Int,Int,Float64,Float64}}

"""
    short_label(pair_label, encounter) -> String

"helper02_vs_target", 1  ->  "h02_vs_t01"
"""
function short_label(pair_label::String, encounter::Int)
    m = match(r"helper(\d+)_vs_target", pair_label)
    hnum = m === nothing ? pair_label : m.captures[1]
    return "h$(hnum)_vs_t$(lpad(encounter, 2, '0'))"
end

"""
    find_encounters(t, rng, max_range) -> Vector{Encounter}

Maximal contiguous index intervals [i0,i1] (into `t`/`rng`) where rng[k] <= max_range.
"""
function find_encounters(t::AbstractVector{<:Real}, rng::AbstractVector{<:Real}, max_range::Real)
    encounters = Encounter[]
    n = length(rng)
    i = 1
    while i <= n
        if rng[i] <= max_range
            j = i
            while j + 1 <= n && rng[j+1] <= max_range
                j += 1
            end
            push!(encounters, (i0=i, i1=j, t0=t[i], t1=t[j]))
            i = j + 1
        else
            i += 1
        end
    end
    return encounters
end

"""
    range_rate_relspeed(t, rt, rh, vt, vh) -> (rng, rate, relspeed)

rng[k]      = |r_target - r_helper|                       [m]
rate[k]     = d(rng)/dt = (Δr·Δv)/|Δr|  (signed; <0 closing, >0 opening)  [m/s]
relspeed[k] = |v_target - v_helper|                        [m/s]
"""
function range_rate_relspeed(rt, rh, vt, vh)
    n = length(rt)
    rng  = Vector{Float64}(undef, n)
    rate = Vector{Float64}(undef, n)
    relspeed = Vector{Float64}(undef, n)
    @inbounds for k in 1:n
        dr = rt[k] - rh[k]
        dv = vt[k] - vh[k]
        r = norm(dr)
        rng[k] = r
        rate[k] = r > 0 ? dot(dr, dv) / r : 0.0
        relspeed[k] = norm(dv)
    end
    return rng, rate, relspeed
end

"""
    plot_relspeed_per_encounter(t, relspeed, encounters, pair_label, outdir, code_label, scenario_label)

Deliverable (1): one PNG per encounter, |v_target - v_helper| vs time-since-encounter-start.
"""
function plot_relspeed_per_encounter(t, relspeed, encounters::Vector{Encounter}, pair_label::String,
                                      outdir::String, code_label::String, scenario_label::String)
    mkpath(outdir)
    n = length(encounters)
    for (k, enc) in enumerate(encounters)
        tt = t[enc.i0:enc.i1] .- t[enc.i0]
        plt = plot(tt, relspeed[enc.i0:enc.i1];
                    xlabel="Time since encounter start (s)",
                    ylabel="Relative speed |v_target - v_helper| (m/s)",
                    title=@sprintf("%s — encounter %d/%d (t=%.0f–%.0f s)\n%s\n%s", pair_label, k, n, enc.t0, enc.t1, code_label, scenario_label),
                    legend=false, lw=2, color=:steelblue, size=(900, 550))
        savefig(plt, joinpath(outdir, @sprintf("%s_encounter%02d_relspeed.png", pair_label, k)))
    end
    return n
end

"""
    plot_range_and_rangerate(t, rng, rate, max_range, pair_label, outdir, code_label, scenario_label)

Deliverable (2): one PNG per pair (full timeline) — range (km, left axis) with two
horizontal reference lines at 0 and +max_range (the within-range band), overlaid with
signed range-rate (m/s, right axis).
"""
function plot_range_and_rangerate(t, rng, rate, max_range::Real, pair_label::String, outdir::String,
                                   code_label::String, scenario_label::String)
    mkpath(outdir)
    rng_km = rng ./ 1000.0
    max_range_km = max_range / 1000.0
    plt = plot(t, rng_km;
               xlabel="Time (s)", ylabel="Range (km)",
               label="range", color=:blue, lw=1.3, size=(900, 600),
               title="$pair_label — range & signed range-rate\n$code_label\n$scenario_label",
               legend=:outerbottom, legend_column=-1)
    hline!(plt, [max_range_km]; label=@sprintf("max range (%.0f km)", max_range_km), ls=:dash, color=:red)
    hline!(plt, [0.0]; label="0 km", ls=:dash, color=:gray)
    plot!(plt, [NaN], [NaN]; label="range-rate", color=:darkgreen, lw=0.9)   # phantom entry for the shared legend row
    plt2 = twinx(plt)
    plot!(plt2, t, rate; ylabel="Signed range-rate (m/s)", label="",
          color=:darkgreen, lw=0.9, legend=false)
    savefig(plt, joinpath(outdir, "$(pair_label)_range_rangerate.png"))
    return plt
end

"""
    contact_window_rows(t, encounters, pair_label, active_mask)

Deliverable (3) data: for each encounter, the geometric contact-window duration (time the
pair stays range-feasible) vs the real laser-on duration (time `active_mask` is true for
this helper, integrated with a left-endpoint rectangle rule over the same samples).
Returns a Vector of NamedTuples (one row per encounter).
"""
function contact_window_rows(t, encounters::Vector{Encounter}, pair_label::String, active_mask::AbstractVector{Bool})
    rows = NamedTuple[]
    for (k, enc) in enumerate(encounters)
        geom_dur = enc.t1 - enc.t0
        laser_dur = 0.0
        for idxk in enc.i0:(enc.i1 - 1)
            dt = t[idxk+1] - t[idxk]
            active_mask[idxk] && (laser_dur += dt)
        end
        push!(rows, (pair=pair_label, encounter=k, t_start_s=enc.t0, t_end_s=enc.t1,
                      geometric_window_s=geom_dur, laser_on_s=laser_dur,
                      duty_cycle=geom_dur > 0 ? laser_dur / geom_dur : 0.0))
    end
    return rows
end

"""
    write_contact_window_table(rows, csv_path, md_path, png_path, code_label, scenario_label)

Writes deliverable (3) as CSV + a small Markdown table, plus a grouped bar chart PNG
(geometric window duration vs real laser-on duration for each encounter). Bar x-axis
labels use the short "h02_vs_t01" form.
"""
function write_contact_window_table(rows::Vector{<:NamedTuple}, csv_path::String, md_path::String,
                                     png_path::String, code_label::String, scenario_label::String)
    mkpath(dirname(csv_path))
    open(csv_path, "w") do io
        println(io, "pair,encounter,t_start_s,t_end_s,geometric_window_s,laser_on_s,duty_cycle")
        for r in rows
            @printf(io, "%s,%d,%.3f,%.3f,%.3f,%.3f,%.4f\n",
                    r.pair, r.encounter, r.t_start_s, r.t_end_s, r.geometric_window_s, r.laser_on_s, r.duty_cycle)
        end
    end
    open(md_path, "w") do io
        println(io, "| Pair | Encounter | t_start (s) | t_end (s) | Geometric window (s) | Laser-on (s) | Duty cycle |")
        println(io, "|---|---|---|---|---|---|---|")
        for r in rows
            @printf(io, "| %s | %d | %.1f | %.1f | %.1f | %.1f | %.1f%% |\n",
                    r.pair, r.encounter, r.t_start_s, r.t_end_s, r.geometric_window_s, r.laser_on_s, 100r.duty_cycle)
        end
    end
    if !isempty(rows)
        labels = [short_label(r.pair, r.encounter) for r in rows]
        geom = [r.geometric_window_s for r in rows]
        laser = [r.laser_on_s for r in rows]
        n = length(rows)
        xs = collect(1:n)
        w = 0.35
        plt = bar(xs .- w/2, geom; bar_width=w, label="geometric window (s)", color=:steelblue,
                  ylabel="Duration (s)", title="Contact-window vs laser-on duration per encounter\n$code_label\n$scenario_label",
                  legend=:outerbottom, legend_column=-1, xticks=(xs, labels), xrotation=60,
                  size=(max(900, 90*n), 600))
        bar!(plt, xs .+ w/2, laser; bar_width=w, label="laser-on (s)", color=:orange)
        mkpath(dirname(png_path))
        savefig(plt, png_path)
    end
    return nothing
end

"""
    plot_relspeed_combined(results, code_label, scenario_label, outdir)

Extra deliverable (4): overlays every range-feasible encounter's relative-speed curve
(one line per encounter, x-axis reset to time-since-encounter-start) on a single plot.
`results` is a Vector of the NamedTuples returned by `run_pair_analysis`.
"""
function plot_relspeed_combined(results, code_label::String, scenario_label::String, outdir::String)
    mkpath(outdir)
    plt = plot(xlabel="Time since encounter start (s)",
               ylabel="Relative speed |v_target - v_helper| (m/s)",
               title="Combined relative speed — all encounters\n$code_label\n$scenario_label",
               legend=:outerbottom, legend_column=-1, size=(950, 600))
    for res in results, (k, enc) in enumerate(res.encounters)
        tt = res.t[enc.i0:enc.i1] .- res.t[enc.i0]
        plot!(plt, tt, res.relspeed[enc.i0:enc.i1]; label=short_label(res.pair_label, k), lw=2)
    end
    savefig(plt, joinpath(outdir, "combined_relspeed_vs_time_since_encounter.png"))
    return plt
end

"""
    plot_range_rangerate_combined(results, max_range, code_label, scenario_label, outdir)

Extra deliverable (5): two stacked panels (range on top, signed range-rate on bottom),
x-axis reset to time-since-encounter-start, one color per encounter, plus a third
(axis-less) strip holding the single shared legend row below both panels. `results` is
a Vector of the NamedTuples returned by `run_pair_analysis`. The 3-line title (main +
code label + property list) is drawn once via `plot_title`, spanning all panels.
"""
function plot_range_rangerate_combined(results, max_range::Real, code_label::String, scenario_label::String, outdir::String)
    mkpath(outdir)
    max_range_km = max_range / 1000.0
    colors = palette(:tab10)
    p1 = plot(ylabel="Range (km)", legend=false)
    hline!(p1, [max_range_km]; label="max range", ls=:dash, color=:red)
    hline!(p1, [0.0]; label="0 km", ls=:dash, color=:gray)
    p2 = plot(xlabel="Time since encounter start (s)", ylabel="Signed range-rate (m/s)", legend=false)
    legend_plot = plot(; framestyle=:none, legend=:inside, legend_column=-1, grid=false, showaxis=false, ticks=false)
    plot!(legend_plot, [NaN], [NaN]; label="max range", ls=:dash, color=:red)
    plot!(legend_plot, [NaN], [NaN]; label="0 km", ls=:dash, color=:gray)
    c_i = 1
    for res in results, (k, enc) in enumerate(res.encounters)
        tt = res.t[enc.i0:enc.i1] .- res.t[enc.i0]
        lbl = short_label(res.pair_label, k)
        c = colors[mod1(c_i, length(colors))]
        plot!(p1, tt, res.rng[enc.i0:enc.i1] ./ 1000.0; label=lbl, color=c, lw=1.6)
        plot!(p2, tt, res.rate[enc.i0:enc.i1]; label=lbl, color=c, lw=1.6)
        plot!(legend_plot, [NaN], [NaN]; label=lbl, color=c, lw=1.6)
        c_i += 1
    end
    plt = plot(p1, p2, legend_plot; layout=grid(3, 1; heights=[0.45, 0.45, 0.1]), size=(950, 800),
               plot_title="Combined range & range-rate — all encounters\n$code_label\n$scenario_label", plot_titlefontsize=10)
    savefig(plt, joinpath(outdir, "combined_range_rangerate_vs_time_since_encounter.png"))
    return plt
end

"""
    run_pair_analysis(t, rt, vt, rh, vh, max_range, pair_label, active_mask, outdir, code_label, scenario_label)

Runs the full per-pair pipeline (deliverables 1 & 2) and returns a NamedTuple carrying
the contact-window rows plus the raw series, so the caller can build the cross-pair
combined plots (deliverables 4 & 5) afterwards.
"""
function run_pair_analysis(t, rt, vt, rh, vh, max_range::Real, pair_label::String,
                            active_mask::AbstractVector{Bool}, outdir::String, code_label::String, scenario_label::String)
    rng, rate, relspeed = range_rate_relspeed(rt, rh, vt, vh)
    encounters = find_encounters(t, rng, max_range)
    plot_relspeed_per_encounter(t, relspeed, encounters, pair_label, joinpath(outdir, "plot1_relative_speed"), code_label, scenario_label)
    plot_range_and_rangerate(t, rng, rate, max_range, pair_label, joinpath(outdir, "plot2_range_and_rangerate"), code_label, scenario_label)
    rows = contact_window_rows(t, encounters, pair_label, active_mask)
    return (rows=rows, t=t, rng=rng, rate=rate, relspeed=relspeed, encounters=encounters, pair_label=pair_label)
end
