module FiringLatitude

include("ComparisonResults.jl")
const CR = ComparisonResults
const SC = CR.SC
using Arrow, CSV, DataFrames, Dates, LinearAlgebra, OrdinaryDiffEq, StaticArrays, TOML
using Plots
include("gravity_kernel.jl")
include(joinpath(CR.CampaignPaths.BASE_CODE, "EncounterTimingRecorder.jl"))

const U_EDGES = deg2rad.(collect(0.0:10.0:360.0))
const RELATIVE_EDGES = deg2rad.(collect(-2.0:0.1:2.0))
const OUTPUT = joinpath(CR.CampaignPaths.OUTPUT, "images", "firing_latitude_h1100_i1_range_only")
const COLORS = ("#652884", "#D9A421")
const NODE_LABEL_TOLERANCE_DEG = 0.05

function geometry(state)
    helper = SVector{3}(state[1:3])
    r, v = SVector{3}(state[7:9]), SVector{3}(state[10:12])
    h = cross(r, v)
    norm(h) > 0 && norm(r) > 0 || error("Degenerate target orbit")
    normal, radial = normalize(h), normalize(r)
    transverse = cross(normal, radial)
    node = cross(SVector(0.0, 0.0, 1.0), normal)
    norm(node) > 1e-10 || error("Argument of latitude requires a nonequatorial target")
    node = normalize(node)
    u = mod(atan(dot(radial, cross(normal, node)), dot(radial, node)), 2pi)
    relative_phase = atan(dot(helper, transverse), dot(helper, radial))
    distance = norm(r - helper)
    distance > 0 || error("Coincident spacecraft in a firing interval")
    (; helper, r, v, h, normal, radial, transverse, node, u, relative_phase, distance,
        helper_normal_m=dot(helper, normal))
end

function mutual_node_u(state)
    g = geometry(state)
    helper_h = cross(g.helper, SVector{3}(state[4:6]))
    norm(helper_h) > 0 || error("Degenerate helper orbit")
    ascending = cross(normalize(helper_h), g.normal)
    norm(ascending) > 1e-10 || error("Coplanar orbits have no unique mutual node")
    atan(dot(ascending, cross(g.normal, g.node)), dot(ascending, g.node))
end

function held_link_dynamics!(du, state, parameters, time)
    g = geometry(state)
    force = parameters.force_N * (g.r - g.helper) / g.distance
    helper_gravity = _inverse_squared_j2_gravity_accel(g.helper, parameters.planet)
    target_gravity = _inverse_squared_j2_gravity_accel(g.r, parameters.planet)
    du[1:3] .= state[4:6]
    du[4:6] .= helper_gravity - force / parameters.mass
    du[7:9] .= g.v
    du[10:12] .= target_gravity + force / parameters.mass
    normal_acceleration = dot(force, g.normal) / parameters.mass
    j2 = target_gravity - _inverse_squared_gravity_accel(g.r, parameters.planet)
    factor = norm(g.r) / norm(g.h) * cos(g.u)
    rate = factor * normal_acceleration
    du[13] = rate
    du[14] = factor * dot(j2, g.normal)
    du[15] = max(rate, 0.0)
    du[16] = min(rate, 0.0)
    du[17] = normal_acceleration
    du[18] = max(normal_acceleration, 0.0)
    du[19] = min(normal_acceleration, 0.0)
    du[20] = cos(g.u) * normal_acceleration
    du[21] = g.relative_phase
    du[22] = g.distance
    du[23] = g.helper_normal_m
    nothing
end

function boundary_state(row, prefix, body)
    [Float64(row[Symbol("$(prefix)_sc_$(body)_$(field)_$axis")])
        for field in ("r", "v") for axis in ("x", "y", "z")]
end

function inclination(state)
    h = cross(SVector{3}(state[1:3]), SVector{3}(state[4:6]))
    atan(hypot(h[1], h[2]), h[3])
end

function bin_index(edges, value)
    isfinite(value) && first(edges) <= value <= last(edges) ||
        error("Angle $value is outside the histogram domain")
    min(searchsortedlast(edges, value), length(edges) - 1)
end

function breakpoints(angles, left, right)
    start, stop = angles(left), angles(right)
    advance = CR.OE.signed_angle(stop.u - start.u)
    0 < advance < U_EDGES[2] - U_EDGES[1] ||
        error("A reconstruction step must advance by less than one latitude bin")
    unwrapped(t) = start.u + CR.OE.signed_angle(angles(t).u - start.u)
    middle = angles((left + right) / 2)
    unwrapped(left) <= unwrapped((left + right) / 2) <= start.u + advance ||
        error("Latitude is not monotone within a reconstruction step")
    lo, hi = extrema((start.relative_phase, stop.relative_phase))
    lo - 1e-10 <= middle.relative_phase <= hi + 1e-10 ||
        error("Relative phase reverses inside a reconstruction step")
    first(RELATIVE_EDGES) <= lo <= hi <= last(RELATIVE_EDGES) ||
        error("Selected helper is outside the relative-phase histogram")
    times = [left, right]
    for boundary in U_EDGES
        start.u < boundary < start.u + advance || continue
        # Scale angular residuals to metres for the shared range-root helper.
        gap(t) = SC.EARTH_RADIUS_M * (unwrapped(t) - boundary)
        push!(times, EncounterTimingRecorder.bisect_crossing(gap, left, right))
    end
    for boundary in RELATIVE_EDGES
        lo < boundary < hi || continue
        gap(t) = SC.EARTH_RADIUS_M * (angles(t).relative_phase - boundary)
        push!(times, EncounterTimingRecorder.bisect_crossing(gap, left, right))
    end
    sort!(unique(times))
end

function accumulate!(bins, joint, joint_di, solution)
    angles(t) = geometry(solution(t))
    for k in 1:length(solution.t)-1
        times = breakpoints(angles, solution.t[k], solution.t[k + 1])
        for j in 1:length(times)-1
            left, right = times[j], times[j + 1]
            g = angles((left + right) / 2)
            latitude_bin = bin_index(U_EDGES, g.u)
            relative_bin = bin_index(RELATIVE_EDGES, g.relative_phase)
            change = solution(right)[13:23] - solution(left)[13:23]
            bins[latitude_bin, 1] += right - left
            bins[latitude_bin, 2:end] .+= change
            joint[latitude_bin, relative_bin] += right - left
            joint_di[latitude_bin, relative_bin] += change[1]
        end
    end
end

function analyze(helpers; max_step_s=10.0)
    case = SC.StudyCase(1100, 1, helpers, :naive_next_entering)
    track = CR.read_target(SC.case_directory(CR.CampaignPaths.OUTPUT, case), case)
    info, settings = track.info, track.info["settings"]
    path = joinpath(track.directory, "laser_firing.feather")
    SC.filehash(path) == info["checksums"]["laser_firing.feather"] ||
        error("Firing table checksum mismatch")
    firing = DataFrame(Arrow.Table(path))
    issorted(firing.start_time_s) || error("Firing intervals are not ordered")
    parameters = (; planet=(μ=SC.MU, J2=SC.J2, Rp_e=SC.EARTH_RADIUS_M),
        mass=SC.MASS_KG, force_N=settings["force_N"])
    bins = zeros(length(U_EDGES) - 1, 12)
    joint = zeros(length(U_EDGES) - 1, length(RELATIVE_EDGES) - 1)
    joint_di = zero(joint)
    total_integrals = zeros(11)
    initial = info["initial_states"][helpers + 1]
    initial_i = inclination(vcat(initial["position_m"], initial["velocity_mps"]))
    previous_i = initial_i
    previous_time, gap_j2 = 0.0, 0.0
    max_position_error, max_velocity_error = 0.0, 0.0
    min_node_u, max_node_u = Inf, -Inf
    for event in eachrow(firing)
        event.sc_b == helpers + 1 && 1 <= event.sc_a <= helpers ||
            error("Unexpected firing endpoints")
        event.start_time_s >= previous_time && event.duration_s > 0 ||
            error("Overlapping or empty firing interval")
        event.force_N == parameters.force_N || error("Firing force differs from run metadata")
        helper = boundary_state(event, "start", "a")
        target = boundary_state(event, "start", "b")
        gap_j2 += inclination(target) - previous_i
        state = vcat(helper, target, zeros(11))
        solution = solve(ODEProblem(held_link_dynamics!, state, (0.0, event.duration_s), parameters),
            Vern9(); reltol=1e-12, abstol=1e-12, dtmax=max_step_s,
            save_everystep=true, dense=true)
        string(solution.retcode) == "Success" || error("Firing reconstruction failed")
        expected = vcat(boundary_state(event, "end", "a"), boundary_state(event, "end", "b"))
        final = last(solution.u)
        for sc in 1:2
            max_position_error = max(max_position_error,
                norm(final[6sc-5:6sc-3] - expected[6sc-5:6sc-3]))
            max_velocity_error = max(max_velocity_error,
                norm(final[6sc-2:6sc] - expected[6sc-2:6sc]))
        end
        node_lo, node_hi = extrema(mutual_node_u, solution.u)
        min_node_u, max_node_u = min(min_node_u, node_lo), max(max_node_u, node_hi)
        accumulate!(bins, joint, joint_di, solution)
        total_integrals .+= final[13:23]
        previous_i = inclination(expected[7:12])
        previous_time = event.end_time_s
    end
    final_i = last(track.elements).i
    gap_j2 += final_i - previous_i
    total_time = sum(bins[:, 1])
    time_error = abs(total_time - info["actual_firing_s"])
    time_error < 1e-5 || error("Latitude histogram does not conserve laser-on time")
    abs(sum(joint) - total_time) < 1e-5 ||
        error("Joint histogram does not conserve laser-on time")
    all(abs.(vec(sum(bins[:, 2:end]; dims=1)) - total_integrals) .<=
        1e-8 .+ 1e-10abs.(total_integrals)) || error("Binning lost diagnostic integrals")
    abs(sum(joint_di) - total_integrals[1]) < 1e-10 ||
        error("Joint histogram lost inclination contributions")
    inclination_residual = rad2deg(total_integrals[1] + total_integrals[2] + gap_j2 -
        (final_i - initial_i))
    abs(inclination_residual) < 1e-7 || error("Laser/J2 inclination attribution does not close")
    max_position_error < 0.01 && max_velocity_error < 1e-5 ||
        error("Reconstructed firing endpoints do not match saved states")
    normal_index = findfirst(==(:dv_N_mps), CR.RTN_FIELDS)
    abs(total_integrals[5] - track.rtn[normal_index, end]) < 1e-6 ||
        error("Reconstructed normal impulse differs from the recorded diagnostic")
    all(>(0), bins[:, 1]) || error("Mean-per-bin diagnostics require nonempty latitude bins")
    max(abs(min_node_u), abs(max_node_u)) < deg2rad(NODE_LABEL_TOLERANCE_DEG) ||
        error("Mutual nodes shifted too far to label them approximately 0 and 180 degrees")
    summary = (; helpers, target_altitude_km=1100, initial_inclination_deg=1.0,
        policy="naive_next_entering", orbits=100,
        firing_intervals=nrow(firing), duration_s=last(track.time_s),
        firing_s=total_time, firing_hours=total_time / 3600,
        duty_fraction=total_time / last(track.time_s),
        outside_range_firing_s=info["outside_range_firing_s"],
        total_di_deg=rad2deg(final_i - initial_i),
        laser_di_deg=rad2deg(total_integrals[1]),
        J2_di_deg=rad2deg(total_integrals[2] + gap_j2),
        positive_laser_di_deg=rad2deg(total_integrals[3]),
        negative_laser_di_deg=rad2deg(total_integrals[4]),
        normal_dv_mps=total_integrals[5],
        phase_weighted_normal_dv_mps=total_integrals[8],
        time_conservation_error_s=time_error,
        inclination_reconstruction_error_deg=inclination_residual,
        max_boundary_position_error_m=max_position_error,
        max_boundary_velocity_error_mps=max_velocity_error,
        mutual_ascending_node_min_deg=rad2deg(min_node_u),
        mutual_ascending_node_max_deg=rad2deg(max_node_u),
        source_trajectory=track.source_trajectory,
        source_trajectory_sha256=info["checksums"]["trajectory.feather"],
        source_firing=path, source_firing_sha256=info["checksums"]["laser_firing.feather"])
    (; helpers, bins, joint, joint_di, summary)
end

function latitude_rows(result)
    b = result.bins
    [(; helpers=result.helpers, latitude_start_deg=rad2deg(U_EDGES[k]),
        latitude_end_deg=rad2deg(U_EDGES[k + 1]),
        latitude_center_deg=rad2deg((U_EDGES[k] + U_EDGES[k + 1]) / 2),
        laser_on_s=b[k, 1], laser_on_hours=b[k, 1] / 3600,
        fraction_of_all_firing=b[k, 1] / result.summary.firing_s,
        laser_di_deg=rad2deg(b[k, 2]), J2_di_while_firing_deg=rad2deg(b[k, 3]),
        positive_laser_di_deg=rad2deg(b[k, 4]), negative_laser_di_deg=rad2deg(b[k, 5]),
        normal_dv_mps=b[k, 6], normal_dv_positive_mps=b[k, 7], normal_dv_negative_mps=b[k, 8],
        phase_weighted_normal_dv_mps=b[k, 9],
        mean_inclination_rate_deg_per_s=rad2deg(b[k, 2]) / b[k, 1],
        mean_normal_acceleration_mps2=b[k, 6] / b[k, 1],
        mean_relative_helper_phase_deg=rad2deg(b[k, 10] / b[k, 1]),
        mean_link_distance_m=b[k, 11] / b[k, 1],
        mean_helper_cross_track_m=b[k, 12] / b[k, 1]) for k in axes(b, 1)]
end

function joint_rows(result)
    [(; helpers=result.helpers,
        latitude_start_deg=rad2deg(U_EDGES[i]), latitude_end_deg=rad2deg(U_EDGES[i + 1]),
        relative_phase_start_deg=rad2deg(RELATIVE_EDGES[j]),
        relative_phase_end_deg=rad2deg(RELATIVE_EDGES[j + 1]),
        laser_on_s=result.joint[i, j], laser_di_deg=rad2deg(result.joint_di[i, j]))
        for i in axes(result.joint, 1) for j in axes(result.joint, 2)]
end

function render(results; output=OUTPUT)
    mkpath(output)
    CSV.write(joinpath(output, "latitude_bins.csv"), vcat(latitude_rows.(results)...))
    CSV.write(joinpath(output, "joint_latitude_relative_phase.csv"), vcat(joint_rows.(results)...))
    CSV.write(joinpath(output, "case_summary.csv"), [result.summary for result in results])
    style = (; fontfamily=CR.PP.FONT_FAMILY, framestyle=:box, gridalpha=0.2,
        tickfontsize=10, guidefontsize=12, titlefontsize=12, legendfontsize=10)
    ticks = collect(0:45:360)
    tick_labels = [angle in (0, 360) ? "$angle\nAsc. node" :
        angle == 180 ? "$angle\nDesc. node" : string(angle) for angle in ticks]
    latitude_axis = (; xlims=(0, 360), xticks=(ticks, tick_labels),
        xlabel="Target argument of latitude u (deg)\nHelper-plane nodes are approximate")
    description = "Target 1100 km, initial i = 1 deg; range-only; 100 reference orbits"
    exposure = plot(; style..., latitude_axis..., size=(1000, 560), dpi=200,
        ylabel="Laser-on time per 10-deg bin (hours)", title=description, legend=:bottomright,
        left_margin=12Plots.mm, bottom_margin=10Plots.mm, top_margin=5Plots.mm,
        right_margin=12Plots.mm)
    contributions = plot(; style..., latitude_axis...,
        ylabel="Laser inclination change (1e-3 deg)", legend=:topright,
        title="Signed inclination contributions: negative values lower inclination")
    efficiency = plot(; style..., latitude_axis..., legend=false,
        ylabel="Mean laser di/dt (1e-8 deg/s)",
        title="Inclination effect per second of firing in each latitude bin")
    centers = rad2deg.((U_EDGES[1:end-1] + U_EDGES[2:end]) / 2)
    for (index, result) in enumerate(results)
        rows = latitude_rows(result)
        hours = getproperty.(rows, :laser_on_hours)
        plot!(exposure, rad2deg.(U_EDGES), vcat(hours, last(hours));
            seriestype=:steppost, lw=2.5, color=COLORS[index], label="$(result.helpers) helpers")
        plot!(contributions, centers, getproperty.(rows, :laser_di_deg) .* 1000;
            lw=2, marker=:circle, ms=3, markerstrokewidth=0,
            color=COLORS[index], label="$(result.helpers) helpers")
        plot!(efficiency, centers, getproperty.(rows, :mean_inclination_rate_deg_per_s) .* 1e8;
            lw=2, marker=:circle, ms=3, markerstrokewidth=0, color=COLORS[index], label="")
    end
    hline!(contributions, [0.0]; color=:gray, ls=:dash, label="")
    hline!(efficiency, [0.0]; color=:gray, ls=:dash, label="")
    phase_centers = rad2deg.((RELATIVE_EDGES[1:end-1] + RELATIVE_EDGES[2:end]) / 2)
    upper = maximum(maximum(result.joint) for result in results) / 3600
    panels = [heatmap(centers, phase_centers, result.joint' / 3600;
        style..., latitude_axis..., ylabel="Relative helper phase (deg)",
        ylims=rad2deg.((first(RELATIVE_EDGES), last(RELATIVE_EDGES))),
        clims=(0, upper), color=:viridis,
        colorbar_title="Laser-on hours", title="$(result.helpers) helpers")
        for result in results]
    clean_exposure = deepcopy(exposure)
    plot!(clean_exposure; title="", xticks=ticks, xlabel="Target argument of latitude u (deg)")
    for panel in (exposure, contributions, efficiency, panels...)
        vline!(panel, [0.0, 180.0, 360.0]; color=:gray, ls=:dash, lw=1.5, label="")
    end
    diagnostic = plot(contributions, efficiency; layout=(2, 1), size=(950, 860), dpi=200,
        margin=6Plots.mm, right_margin=12Plots.mm)
    geometry_plot = plot(panels...; layout=(1, 2), size=(1150, 530), dpi=200,
        margin=5Plots.mm, bottom_margin=12Plots.mm)
    combined_figures = map((exposure, clean_exposure)) do source
        combined = deepcopy(source)
        plot!(combined; size=(1150, 660), legend=false,
            title=source === exposure ?
                description * "\nSolid steps: laser-on time (left); dashed circles: inclination change (right)" : "",
            right_margin=22Plots.mm,
            ylims=(0, ceil(1.1 * maximum(maximum(result.bins[:, 1]) for result in results) / 3600)))
        inclination_axis = twinx(combined)
        inclination_limit = ceil(1.1 * maximum(maximum(abs, rad2deg.(result.bins[:, 2])) for result in results) * 1000;
            digits=1)
        plot!(inclination_axis; fontfamily=CR.PP.FONT_FAMILY, tickfontsize=10, guidefontsize=12,
            xlims=(0, 360), ylabel="Laser inclination change per 10-deg bin (1e-3 deg)",
            ylims=(-inclination_limit, inclination_limit), legend=false, grid=false)
        hline!(inclination_axis, [0.0]; color=:gray, ls=:dot, label="")
        for (index, result) in enumerate(results)
            plot!(inclination_axis, centers, getproperty.(latitude_rows(result), :laser_di_deg) .* 1000;
                lw=2, ls=:dash, marker=:circle, ms=3, markerstrokewidth=0,
                color=COLORS[index], label="")
        end
        # Draw actual series samples: GR's automatic legend enlarges markers and lines.
        legend_panel = plot(; fontfamily=CR.PP.FONT_FAMILY,
            framestyle=:none, showaxis=false, grid=false, xlims=(0, 1), ylims=(0, 1),
            legend=false, margin=0Plots.mm)
        plot!(legend_panel, [0.005, 0.995, 0.995, 0.005, 0.005],
            [0.25, 0.25, 0.8, 0.8, 0.25]; color=:black, lw=1, label="")
        for (index, result) in enumerate(results)
            x, y = (0.01, 0.24)[index], 0.525
            plot!(legend_panel, x .+ [0.015, 0.06], [y, y]; lw=2.5, color=COLORS[index], label="")
            annotate!(legend_panel, x + 0.075, y, text("$(result.helpers) helpers: laser-on time",
                10, :left, CR.PP.FONT_FAMILY))
        end
        for (index, result) in enumerate(results)
            x, y = (0.47, 0.73)[index], 0.525
            plot!(legend_panel, x .+ [0.015, 0.06], [y, y]; lw=2, ls=:dash, color=COLORS[index], label="")
            scatter!(legend_panel, [x + 0.0375], [y]; marker=:circle, ms=3, markerstrokewidth=0,
                color=COLORS[index], label="")
            annotate!(legend_panel, x + 0.075, y, text("$(result.helpers) helpers: inclination change",
                10, :left, CR.PP.FONT_FAMILY))
        end
        plot(combined, legend_panel; layout=grid(2, 1; heights=[0.9, 0.1]),
            size=(1400, 760), dpi=200)
    end
    figures = (("laser_on_time_vs_argument_of_latitude", exposure),
        ("inclination_contribution_vs_argument_of_latitude", diagnostic),
        ("joint_target_latitude_helper_relative_phase", geometry_plot),
        ("laser_on_time_and_inclination_vs_argument_of_latitude", combined_figures[1]),
        ("laser_on_time_and_inclination_vs_argument_of_latitude_clean", combined_figures[2]))
    paths = String[]
    for (name, figure) in figures, extension in ("png", "pdf")
        path = joinpath(output, "$name.$extension")
        savefig(figure, path)
        filesize(path) > 10_000 || error("Empty firing-latitude figure: $path")
        push!(paths, path)
        println("Saved $path")
    end
    paths
end

function main()
    isempty(ARGS) || error("Usage: plot_firing_latitude.jl")
    mkpath(OUTPUT)
    SC.atomic_toml(joinpath(OUTPUT, "validation.toml"),
        Dict("passed" => false, "status" => "running", "started_utc" => string(now(UTC))))
    results, convergence = NamedTuple[], Dict{String,Any}()
    for helpers in (300, 800)
        println("Reconstructing $helpers-helper firing intervals with 10 s step cap")
        flush(stdout)
        coarse = analyze(helpers)
        println("Checking $helpers-helper exposure with 5 s step cap")
        flush(stdout)
        fine = analyze(helpers; max_step_s=5.0)
        time_error = maximum(abs.(coarse.bins[:, 1] - fine.bins[:, 1]))
        joint_error = maximum(abs.(coarse.joint - fine.joint))
        di_error = rad2deg(maximum(abs.(coarse.bins[:, 2] - fine.bins[:, 2])))
        node_error = max(
            abs(coarse.summary.mutual_ascending_node_min_deg - fine.summary.mutual_ascending_node_min_deg),
            abs(coarse.summary.mutual_ascending_node_max_deg - fine.summary.mutual_ascending_node_max_deg))
        time_error < 0.01 && joint_error < 0.01 && di_error < 1e-8 ||
            error("Firing-latitude reconstruction did not converge")
        node_error < 1e-4 || error("Sampled mutual-node extrema did not converge")
        convergence[string(helpers)] = Dict("max_bin_time_change_s" => time_error,
            "max_joint_bin_time_change_s" => joint_error, "max_bin_inclination_change_deg" => di_error,
            "max_node_extrema_change_deg" => node_error)
        push!(results, fine)
        println(fine.summary)
        flush(stdout)
    end
    paths = render(results)
    SC.atomic_toml(joinpath(OUTPUT, "validation.toml"), Dict(
        "passed" => true, "completed_utc" => string(now(UTC)), "orbits" => 100,
        "target_altitude_km" => 1100, "target_initial_inclination_deg" => 1,
        "helper_counts" => [300, 800], "policy" => "naive_next_entering",
        "latitude_bin_width_deg" => 10, "relative_phase_bin_width_deg" => 0.1,
        "method" => "Replay each recorded held-link interval from its saved boundary states; split time at dense-output angular-bin crossings; integrate on-time and Gauss inclination rate.",
        "latitude_definition" => "u = atan2(r_hat dot (N_hat cross node_hat), r_hat dot node_hat), modulo 360 deg; instantaneous target orbital plane",
        "relative_phase_definition" => "atan2(r_helper dot T_target, r_helper dot R_target); positive is helper ahead; helper projected into target plane",
        "orbit_node_definition" => "Ascending mutual node is parallel to N_helper cross N_target; target crosses toward the positive helper-plane normal. Descending node is opposite.",
        "orbit_node_markers_deg" => [0, 180, 360],
        "orbit_node_marker_tolerance_deg" => NODE_LABEL_TOLERANCE_DEG,
        "orbit_node_sampling" => "Selected-helper and target osculating planes at every replay solver state during firing; plotted case summaries use the 5 s step cap. Signed ascending-node offsets are near 0 deg.",
        "force_model" => "Recorded BP/c force, equal/opposite on selected helper and target; includes held out-of-range portions",
        "physical_sources" => ["gravity_kernel.jl", "recorded run_info.toml and laser_firing.feather"],
        "not_event_start_histogram" => true,
        "solver" => "Vern9, reltol=abstol=1e-12; 10 s versus 5 s convergence; plotted data use 5 s",
        "convergence" => convergence,
        "source_code_checksums" => Dict(path => SC.filehash(path) for path in
            (@__FILE__, joinpath(@__DIR__, "gravity_kernel.jl"),
                joinpath(CR.CampaignPaths.BASE_CODE, "EncounterTimingRecorder.jl"))),
        "figure_checksums" => Dict(basename(path) => SC.filehash(path) for path in paths),
        "csv_checksums" => Dict(file => SC.filehash(joinpath(OUTPUT, file))
            for file in ("latitude_bins.csv", "joint_latitude_relative_phase.csv", "case_summary.csv"))))
    println("Firing-latitude plots and numerical checks completed.")
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    FiringLatitude.main()
end
