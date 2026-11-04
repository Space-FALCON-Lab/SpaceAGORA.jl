using LinearAlgebra
using Statistics
using TOML
using SHA
using Test

const INTEGRATION_ROOT = normpath(joinpath(@__DIR__, "..", "9_SpaceAGORA.jl-main_full integration"))
include(joinpath(INTEGRATION_ROOT, "examples", "common.jl"))
using Arrow
using CSV
using DataFrames
using StaticArrays
using DiffEqBase

const REFERENCE_DIR = normpath(joinpath(@__DIR__, "..", "..", "SpaceAGORA.jl",
    "output", "paper_plot_mode", "h1000km_t1050km", "ih0.0deg_it1.0deg",
    "N50", "T509707s", "gve_sma_e0.0000_nu1.0000"))
const RESULT_DIR = joinpath(@__DIR__, "results", "saved_paper_h1000_t1050_i1_N50_80orbits")
const REFERENCE_FILE = joinpath(REFERENCE_DIR, "simulation_results.feather")
const MASS = 227.0
const SPACECRAFT_COUNT = 51

readvec(table, sc, field, row) =
    SVector{3, Float64}(ntuple(k -> table[Symbol("sc$(sc)_$(field)_$(k)")][row], 3))
filehash(path) = open(io -> bytes2hex(sha256(io)), path)

function elements(r, v, mu)
    h = cross(r, v)
    eccentricity = norm(cross(v, h) / mu - r / norm(r))
    return (sma_m=inv(2 / norm(r) - dot(v, v) / mu),
        eccentricity=eccentricity,
        inclination_deg=rad2deg(atan(hypot(h[1], h[2]), h[3])),
        raan_deg=rad2deg(atan(h[1], -h[2])))
end

function reference()
    table = Arrow.Table(REFERENCE_FILE)
    manifest = TOML.parsefile(joinpath(REFERENCE_DIR, "simulation_results.manifest.toml"))
    @test filehash(REFERENCE_FILE) == manifest["files"]["feather"]["sha256"]
    @test length(table.time) == manifest["steps"] == 1001
    @test first(table.time) == 0
    @test last(table.time) == manifest["mission_time_s"]
    @test manifest["spacecraft_count"] == SPACECRAFT_COUNT
    @test all(table[Symbol("sc$(i)_mass")][1] == MASS for i in 1:SPACECRAFT_COUNT)
    return table
end

function configuration(table, dtmax, duration, output; legacy_radius=false)
    planet = make_no_gram_planet(:earth)
    if legacy_radius
        # Reproduce the historical J2*R_mean^2 term without changing library code.
        planet = typeof(planet)(; J2=planet.J2 * (planet.Rp_m / planet.Rp_e)^2)
    end
    spacecraft = map(1:SPACECRAFT_COUNT) do i
        bus = Link(root=true, m=MASS)
        initial = CartesianInitialCondition(readvec(table, i, "pos", 1),
            readvec(table, i, "vel", 1))
        SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia,
            initial_condition=initial, id=i, n_terminal=1)
    end
    model = ConstellationInterLinkModel()
    for helper in 2:SPACECRAFT_COUNT
        register_allowed_link!(model, spacecraft, (1, 1), (helper, 1);
            parameters=InterLinkParameters(P=10_000.0, B=100.0, range=200e3))
    end
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=true, generate_plots=false,
            results_directory=output, save_csv=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=duration, orientation_sim=false, data_rate=10.0),
        environment_model=make_no_gram_environment(planet=planet),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredJ2GravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2026, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-12,
            abstol_orbit=1e-12, dt_max_orbit=dtmax),
        solver_config=SM.SolverConfig(solver_mode=:tsit5),
        interlink_model=model,
        scheduling_policy_model=SchedulingPolicyModel(:gve_sma; target_idx=1))
end

function diagnostic_callback(table)
    previous_r = readvec(table, 1, "pos", 1)
    previous_v = [readvec(table, i, "vel", 1) for i in 1:SPACECRAFT_COUNT]
    previous_dv = fill(zero(SVector{3, Float64}), SPACECRAFT_COUNT)
    totals = zeros(7)
    rows = DataFrame(time_s=Float64[], dv_R_mps=Float64[], dv_T_mps=Float64[],
        dv_N_mps=Float64[], dv_T_positive_mps=Float64[], dv_T_negative_mps=Float64[],
        target_work_J=Float64[], helpers_work_J=Float64[])
    push!(rows, (0.0, totals...))
    function affect!(integrator)
        state = integrator.u.sc
        r = SVector{3, Float64}(state[1].pos)
        v = SVector{3, Float64}(state[1].vel)
        midpoint_r = (previous_r + r) / 2
        midpoint_v = (previous_v[1] + v) / 2
        radial = normalize(midpoint_r)
        normal = normalize(cross(midpoint_r, midpoint_v))
        tangent = cross(normal, radial)
        increment = SVector{3, Float64}(state[1].laser_dv) - previous_dv[1]
        dR, dT, dN = dot(increment, radial), dot(increment, tangent), dot(increment, normal)
        totals[1] += dR
        totals[2] += dT
        totals[3] += dN
        totals[4] += max(dT, 0)
        totals[5] += min(dT, 0)
        for i in 1:SPACECRAFT_COUNT
            dv = SVector{3, Float64}(state[i].laser_dv)
            vel = SVector{3, Float64}(state[i].vel)
            work = MASS * dot(dv - previous_dv[i], (vel + previous_v[i]) / 2)
            totals[i == 1 ? 6 : 7] += work
            previous_dv[i] = dv
            previous_v[i] = vel
        end
        previous_r = r
        push!(rows, (Float64(integrator.t), totals...))
        return nothing
    end
    # Project solver-integrated ECI impulse increments at the step midpoint.
    # These are quadrature estimates, not the planned RTN ODE states.
    callback = DiffEqBase.DiscreteCallback((u, t, integrator) -> true, affect!;
        save_positions=(false, false))
    return callback, rows
end

function run_case(dtmax; smoke=false, legacy_radius=false)
    old = reference()
    duration = smoke ? 120.0 : last(old.time)
    label = legacy_radius ? "legacy_mean_radius_dt$(dtmax)s" : "dt$(dtmax)s"
    output = joinpath(RESULT_DIR, smoke ? "smoke" : label)
    isdir(output) && !isempty(readdir(output)) &&
        error("Refusing to overwrite existing results: $output")
    mkpath(output)
    args = configuration(old, dtmax, duration, output; legacy_radius)
    planet = args.environment_model.planet
    initial_elements = elements(readvec(old, 1, "pos", 1), readvec(old, 1, "vel", 1), planet.μ)
    @test isapprox(initial_elements.sma_m, planet.Rp_e + 1050e3; atol=1e-6)
    @test isapprox(initial_elements.inclination_deg, 1.0; atol=1e-12)
    @test isapprox(last(old.time) / (2pi * sqrt(initial_elements.sma_m^3 / planet.μ)),
        80.0; atol=1e-10)
    callback, diagnostics = diagnostic_callback(old)
    fields = filter(f -> f.name in (:position, :velocity, :laser_dv, :laser_delta_sma),
        SM.default_save_fields(args))
    println("Running ", smoke ? "smoke test" : "full saved case",
        ": dtmax=", dtmax, ", duration=", duration, ", J2=", planet.J2)
    flush(stdout)
    solution = withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false") do
        run_simulation(args; return_solution=true, save_fields=fields,
            extra_callbacks=(callback,))
    end
    @test string(solution.retcode) == "Success"
    @test solution.t[end] == duration
    @test unique(solution.t) == [0.0, duration]
    @test all(isfinite, solution.u[end])
    for i in 1:SPACECRAFT_COUNT
        @test solution.u[1].sc[i].pos == readvec(old, i, "pos", 1)
        @test solution.u[1].sc[i].vel == readvec(old, i, "vel", 1)
    end
    history = solution.prob.p.args.interlink_model.history
    @test length(history) == solution.stats.naccept + 1
    @test length(diagnostics.time_s) == length(history)
    @test diagnostics.time_s == [sample.time for sample in history]
    @test maximum(diff(diagnostics.time_s)) <= dtmax + 2eps(duration)
    @test all(length(sample.active) <= 1 for sample in history)
    @test isapprox(diagnostics.dv_T_mps[end],
        diagnostics.dv_T_positive_mps[end] + diagnostics.dv_T_negative_mps[end]; atol=1e-12)
    @test norm(sum(MASS * SVector{3, Float64}(s.laser_dv) for s in solution.u[end].sc)) < 1e-8
    @test solution.u[end].sc[1].laser_delta_sma > 0
    Arrow.write(joinpath(output, "accepted_step_diagnostics.feather"), diagnostics)
    schedule = DataFrame(time_s=[s.time for s in history],
        active_helper=[isempty(s.active) ? 0 : s.active[1][2][1] for s in history])
    Arrow.write(joinpath(output, "accepted_step_schedule.feather"), schedule)
    final = solution.u[end].sc[1]
    final_elements = elements(final.pos, final.vel, planet.μ)
    summary = Dict{String, Any}(
        "retcode" => string(solution.retcode), "duration_s" => duration,
        "dtmax_s" => dtmax, "mu_m3_s2" => planet.μ, "earth_radius_m" => planet.Rp_e,
        "J2" => planet.J2, "helper_count" => 50, "mass_kg" => MASS,
        "legacy_mean_radius_diagnostic_only" => legacy_radius,
        "earth_mean_radius_m" => planet.Rp_m,
        "power_W" => 10_000.0, "B" => 100.0, "range_m" => 200e3,
        "target_initial_phase_deg" => 1.0, "scheduler" => "gve_sma",
        "epoch" => "2026-01-01T00:00:00", "reltol" => 1e-12, "abstol" => 1e-12,
        "accepted_steps" => solution.stats.naccept, "rejected_steps" => solution.stats.nreject,
        "solver_saved_times_s" => collect(solution.t),
        "maximum_step_s" => maximum(diff(schedule.time_s)),
        "median_step_s" => median(diff(schedule.time_s)),
        "applied_firing_time_s" => sum(diff(schedule.time_s) .* (schedule.active_helper[1:end-1] .> 0)),
        "target_laser_delta_sma_m" => final.laser_delta_sma,
        "target_laser_dv_ECI_mps" => collect(final.laser_dv),
        "reference_sha256" => filehash(REFERENCE_FILE),
        "test_script_sha256" => filehash(@__FILE__),
        "diagnostic_method" => "Accepted-step ECI-impulse midpoint projection; approximate RTN/work, not RTN ODE states",
        "reference_limit" => "Initial states and duration verified from file; laser settings, epoch, gravity and tolerances inferred from current old runner defaults, not archived run settings")
    for key in keys(initial_elements)
        summary["initial_$(key)"] = initial_elements[key]
        summary["final_$(key)"] = final_elements[key]
        delta = final_elements[key] - initial_elements[key]
        summary["delta_$(key)"] = key == :raan_deg ? mod(delta + 180, 360) - 180 : delta
    end
    for name in names(diagnostics)[2:end]
        summary[name] = diagnostics[end, name]
    end
    open(joinpath(output, "summary.toml"), "w") do io
        TOML.print(io, summary; sorted=true)
    end
    validate_output(output, old, duration)
    println("Saved outputs to ", output, "; see the test summary for validation status.")
    return nothing
end

function validate_output(output, old, duration)
    table = Arrow.Table(joinpath(output, "simulation_results.feather"))
    expected = collect(0.0:10.0:duration)
    last(expected) < duration && push!(expected, duration)
    @test collect(table.time) == expected
    @test all(all(isfinite, table[n]) for n in propertynames(table))
    for i in 1:SPACECRAFT_COUNT
        @test readvec(table, i, "pos", 1) == readvec(old, i, "pos", 1)
        @test readvec(table, i, "vel", 1) == readvec(old, i, "vel", 1)
    end
end

function target_series(table, mu)
    values = [elements(readvec(table, 1, "pos", i), readvec(table, 1, "vel", i), mu)
        for i in eachindex(table.time)]
    return DataFrame(time_s=Float64.(table.time),
        sma_m=getproperty.(values, :sma_m), eccentricity=getproperty.(values, :eccentricity),
        inclination_deg=getproperty.(values, :inclination_deg), raan_deg=getproperty.(values, :raan_deg))
end

function compare()
    old = reference()
    labels = ("dt10.0s", "dt5.0s", "dt2.5s", "legacy_mean_radius_dt5.0s")
    summaries = Dict(label => TOML.parsefile(joinpath(RESULT_DIR, label, "summary.toml"))
        for label in labels)
    mu = summaries["dt10.0s"]["mu_m3_s2"]
    old_start = elements(readvec(old, 1, "pos", 1), readvec(old, 1, "vel", 1), mu)
    old_end = elements(readvec(old, 1, "pos", length(old.time)),
        readvec(old, 1, "vel", length(old.time)), mu)
    CSV.write(joinpath(RESULT_DIR, "reference_target_elements.csv"), target_series(old, mu))
    metrics = DataFrame(metric=String[], old_saved=Float64[], new_dt10=Float64[],
        new_dt5=Float64[], new_dt2_5=Float64[], legacy_radius_dt5=Float64[])
    for key in keys(old_start)
        change = old_end[key] - old_start[key]
        key == :raan_deg && (change = mod(change + 180, 360) - 180)
        push!(metrics, ("delta_$(key)", change,
            (summaries[label]["delta_$(key)"] for label in labels)...))
    end
    for (axis, oldcolumn) in zip(("R", "T", "N"),
        (:dv_r_accumulated, :dv_t_accumulated, :dv_n_accumulated))
        key = "dv_$(axis)_mps"
        push!(metrics, (key, old[oldcolumn][end],
            (summaries[label][key] for label in labels)...))
    end
    metrics.new10_minus_old = metrics.new_dt10 .- metrics.old_saved
    metrics.new5_minus_new10 = metrics.new_dt5 .- metrics.new_dt10
    metrics.new2_5_minus_new5 = metrics.new_dt2_5 .- metrics.new_dt5
    CSV.write(joinpath(RESULT_DIR, "comparison.csv"), metrics)
    tables = Dict(label => Arrow.Table(joinpath(RESULT_DIR, label, "simulation_results.feather"))
        for label in keys(summaries))
    differences = Dict{String, Any}("reference_rows" => length(old.time),
        "reference_spacing_s" => old.time[2] - old.time[1],
        "comparison_is_validation_of_one_case_only" => true)
    for label in sort!(collect(keys(tables)))
        table = tables[label]
        validate_output(joinpath(RESULT_DIR, label), old, last(old.time))
        @test summaries[label]["retcode"] == "Success"
        @test summaries[label]["duration_s"] == last(old.time)
        diagnostics = Arrow.Table(joinpath(RESULT_DIR, label, "accepted_step_diagnostics.feather"))
        schedule = Arrow.Table(joinpath(RESULT_DIR, label, "accepted_step_schedule.feather"))
        @test diagnostics.time_s == schedule.time_s
        @test first(schedule.time_s) == 0.0
        @test last(schedule.time_s) == last(old.time)
        @test length(schedule.time_s) == summaries[label]["accepted_steps"] + 1
        @test all(>(0), diff(schedule.time_s))
        @test maximum(diff(schedule.time_s)) <= summaries[label]["dtmax_s"] + 2eps(last(old.time))
        @test all(i == 0 || 2 <= i <= SPACECRAFT_COUNT for i in schedule.active_helper)
        for name in propertynames(diagnostics)[2:end]
            @test all(isfinite, diagnostics[name])
            @test first(diagnostics[name]) == 0.0
            @test last(diagnostics[name]) == summaries[label][string(name)]
        end
        @test isapprox(sum(diff(schedule.time_s) .* (schedule.active_helper[1:end-1] .> 0)),
            summaries[label]["applied_firing_time_s"]; atol=1e-9)
        CSV.write(joinpath(RESULT_DIR, "$(label)_target_elements.csv"), target_series(table, mu))
        differences["$(label)_final_target_position_difference_from_old_m"] =
            norm(readvec(table, 1, "pos", length(table.time)) - readvec(old, 1, "pos", length(old.time)))
        differences["$(label)_final_target_velocity_difference_from_old_mps"] =
            norm(readvec(table, 1, "vel", length(table.time)) - readvec(old, 1, "vel", length(old.time)))
    end
    a, b = tables["dt10.0s"], tables["dt5.0s"]
    @test a.time == b.time
    differences["dt5_vs_dt10_final_target_position_difference_m"] =
        norm(readvec(a, 1, "pos", length(a.time)) - readvec(b, 1, "pos", length(b.time)))
    differences["dt5_vs_dt10_final_target_velocity_difference_mps"] =
        norm(readvec(a, 1, "vel", length(a.time)) - readvec(b, 1, "vel", length(b.time)))
    differences["dt5_vs_dt10_max_sampled_target_position_difference_m"] =
        maximum(norm(readvec(a, 1, "pos", i) - readvec(b, 1, "pos", i)) for i in eachindex(a.time))
    differences["dt5_vs_dt10_laser_sma_difference_m"] =
        summaries["dt5.0s"]["target_laser_delta_sma_m"] - summaries["dt10.0s"]["target_laser_delta_sma_m"]
    c = tables["dt2.5s"]
    @test b.time == c.time
    differences["dt2_5_vs_dt5_final_target_position_difference_m"] =
        norm(readvec(b, 1, "pos", length(b.time)) - readvec(c, 1, "pos", length(c.time)))
    differences["dt2_5_vs_dt5_final_target_velocity_difference_mps"] =
        norm(readvec(b, 1, "vel", length(b.time)) - readvec(c, 1, "vel", length(c.time)))
    differences["dt2_5_vs_dt5_laser_sma_difference_m"] =
        summaries["dt2.5s"]["target_laser_delta_sma_m"] - summaries["dt5.0s"]["target_laser_delta_sma_m"]
    open(joinpath(RESULT_DIR, "comparison.toml"), "w") do io
        TOML.print(io, differences; sorted=true)
    end
    show(stdout, metrics; allrows=true, allcols=true)
    println()
    TOML.print(stdout, differences; sorted=true)
end

@testset "Saved paper case comparison" begin
    if ARGS == ["smoke"]
        run_case(10.0; smoke=true)
    elseif length(ARGS) == 2 && ARGS[1] == "run"
        dtmax = parse(Float64, ARGS[2])
        isfinite(dtmax) && 0 < dtmax <= 10 || error("Step cap must be in (0, 10] seconds.")
        run_case(dtmax)
    elseif ARGS == ["legacy-radius"]
        run_case(5.0; legacy_radius=true)
    elseif ARGS == ["compare"]
        compare()
    else
        error("Usage: compare_saved_paper_case.jl smoke | run DTMAX_SECONDS | legacy-radius | compare")
    end
end
