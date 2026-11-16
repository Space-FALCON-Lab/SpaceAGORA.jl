#!/usr/bin/env julia
# From ORACLE_SciTech_2027:
# julia --project="9_SpaceAGORA.jl-main_full integration" --threads=1 Tests/prototype_validation/compare_outputs.jl audit
# julia --project="9_SpaceAGORA.jl-main_full integration" --threads=1 Tests/prototype_validation/compare_outputs.jl replay
# Reads production outputs; writes only this test's results directory.
# Replay uses recorded firing times and matched Earth constants, not prototype
# defaults. Native tests independently gate forces; ver3 also needs half power
# to match production's BP/c rather than the prototype's 2BP/c convention.
using Arrow, CSV, DataFrames, Dates, DiffEqCallbacks, LinearAlgebra
using OrdinaryDiffEq, Printf, SHA, StaticArrays, TOML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PROTOTYPE = normpath(joinpath(ROOT, "..", "Kuang's Prototype Code"))
const RESULTS = joinpath(@__DIR__, "results")
const OUTPUT3 = joinpath(ROOT, "Output_3_rangeOnly_vs_positiveAlongTrack", "Output")
const OUTPUT4 = joinpath(ROOT, "Output_4_no_laser", "Output")
const OUTPUT5 = joinpath(ROOT, "Output_5_Enable_Disable_Laser_compare", "Output")

module Prototype
using LinearAlgebra, StaticArrays
const MU = 3.98600436233e14
const R_ATMDEF = 6378136.6 + 100_000
const Ẑ = SVector(0.0, 0.0, 1.0)
const Ẑ = Ẑ
idx(i, off) = 6(i - 1) + off
include(joinpath(Main.PROTOTYPE, "functions", "1_LOS_Metrics.jl"))
include(joinpath(Main.PROTOTYPE, "functions", "2_Laser_Forces_ver2.jl"))
include(joinpath(Main.PROTOTYPE, "functions", "3_Dynamics.jl"))
include(joinpath(Main.PROTOTYPE, "functions", "5_OE_Converters.jl"))
end

module ForwardPrototype
using LinearAlgebra, StaticArrays
const R_ATMDEF = 6378136.6 + 100_000
idx(i, off) = 6(i - 1) + off
include(joinpath(Main.PROTOTYPE, "functions", "1_LOS_Metrics.jl"))
include(joinpath(Main.PROTOTYPE, "functions", "2_Laser_Forces_ver3.jl"))
include(joinpath(Main.PROTOTYPE, "functions", "3_Dynamics.jl"))
end

const RESIDUALS = Dict{String,Float64}()
const CHECK_COUNTS = Dict{String,Int}()

function compare(name, actual, expected, tolerance)
    a, b = actual isa Number ? [actual] : actual, expected isa Number ? [expected] : expected
    length(a) == length(b) || error("$name: length mismatch")
    residual = 0.0
    for (x, y) in zip(a, b)
        if isnan(x) || isnan(y)
            isnan(x) && isnan(y) || error("$name: undefined-value mismatch")
        else
            isfinite(x) && isfinite(y) || error("$name: nonfinite value")
            residual = max(residual, abs(x - y))
        end
    end
    RESIDUALS[name] = max(get(RESIDUALS, name, 0.0), residual)
    CHECK_COUNTS[name] = get(CHECK_COUNTS, name, 0) + length(a)
    residual <= tolerance || error("$name: residual $residual exceeds $tolerance")
    residual
end

filehash(path) = open(io -> bytes2hex(sha256(io)), path)
key(row) = (row.target_altitude_km, row.delta_i_deg, row.helpers)
casekey(row) = (row.scheduler, key(row)...)
endpointkey(row) = (casekey(row)..., row.orbits)
summary(root) = CSV.read(joinpath(root, "production_endpoint_summary.csv"), DataFrame)

function elements(r, v, mu; equatorial=false)
    raw = [Prototype.rv2coe(r[k], v[k], mu) for k in eachindex(r)]
    raan = equatorial ? fill(NaN, length(raw)) : [first(raw).Ω]
    if !equatorial
        for k in 2:length(raw)
            delta = raw[k].Ω - raw[k - 1].Ω
            push!(raan, raan[end] + atan(sin(delta), cos(delta)))
        end
    end
    (; a=getproperty.(raw, :a), e=getproperty.(raw, :e),
        i=rad2deg.(getproperty.(raw, :i)), raan=rad2deg.(raan))
end

function read_track(path)
    directory = dirname(path)
    info = TOML.parsefile(joinpath(directory, "run_info.toml"))
    settings = info["settings"]
    verification = TOML.parsefile(joinpath(directory, "verification.toml"))
    info["status"] == "verified" && verification["passed"] || error("Unverified run $path")
    digest = filehash(path)
    digest == info["checksums"]["trajectory.feather"] || error("Checksum mismatch: $path")
    tbl = Arrow.Table(path)
    n = settings["target_index"]
    r = [SVector{3,Float64}(ntuple(j -> tbl[Symbol("sc$(n)_pos_$j")][k], 3))
        for k in 1:length(tbl.time)]
    v = [SVector{3,Float64}(ntuple(j -> tbl[Symbol("sc$(n)_vel_$j")][k], 3))
        for k in 1:length(tbl.time)]
    time = collect(tbl.time)
    compare("exact_time_grid_s", time,
        collect(range(0.0, settings["duration_s"]; length=1001)), 0.0)
    compare("initial_position_m", r[1], info["initial_states"][n]["position_m"], 0.0)
    compare("initial_velocity_mps", v[1], info["initial_states"][n]["velocity_mps"], 0.0)
    diagnostics = Dict(string(name)[8:end] => collect(tbl[name])
        for name in propertynames(tbl) if startswith(string(name), "target_"))
    for axis in ("R", "T", "N")
        positive = diagnostics["dv_$(axis)_positive_mps"]
        negative = diagnostics["dv_$(axis)_negative_mps"]
        compare("RTN_branch_sum_mps", diagnostics["dv_$(axis)_mps"], positive + negative, 1e-8)
        diagnostics["dv_$(axis)_absolute_mps"] = positive - negative
    end
    if settings["scheduler"] == "no_laser"
        info["firing_intervals"] == 0 && info["actual_firing_s"] == 0 ||
            error("No-laser run contains firing")
        for name in propertynames(tbl)
            if occursin("laser", string(name)) || startswith(string(name), "target_")
                all(iszero, tbl[name]) || error("Nonzero baseline diagnostic: $path $name")
            end
        end
    end
    mu, j2, re, mass = (settings[x] for x in
        ("mu_m3_s2", "J2", "earth_equatorial_radius_m", "mass_kg"))
    energy = [mass * (dot(v[k], v[k]) / 2 - mu / norm(r[k]) +
        mu * j2 * re^2 / (2norm(r[k])^3) * (3(r[k][3] / norm(r[k]))^2 - 1))
        for k in eachindex(r)]
    work = collect(tbl[Symbol("sc$(n)_laser_work_J")])
    compare("target_J2_energy_minus_work_J", energy .- energy[1], work,
        20.0 + 1e-4maximum(abs, work))
    oe = elements(r, v, mu; equatorial=settings["target_inclination_deg"] == 0)
    (; directory, info, settings, digest, r, v, time, diagnostics, oe)
end

const OE_FIELDS = (("da_m", :a, 1e-5), ("de", :e, 1e-12),
    ("di_deg", :i, 2e-6), ("draan_deg", :raan, 1e-8))

function check_endpoints(track, rows)
    for row in eachrow(rows)
        row.run_id == track.info["run_id"] &&
            row.source_trajectory_sha256 == track.digest || error("Endpoint provenance mismatch")
        index = Int(row.row)
        index == (row.orbits == 80 ? 801 : 1001) || error("Unexpected endpoint index")
        compare("endpoint_time_s", row.time_s, track.time[index], 0.0)
        compare("endpoint_orbit_count", track.time[index] /
            (track.settings["duration_s"] / 100), row.orbits, 1e-12)
        for (column, field, tol) in OE_FIELDS
            values = getproperty(track.oe, field)
            compare("endpoint_$column", row[Symbol(column)], values[index] - values[1], tol)
        end
        for (column, values) in track.diagnostics
            compare("endpoint_$column", row[Symbol(column)], values[index], 1e-10)
        end
    end
end

function check_difference(on, off, rows)
    allowed = Set(("scheduler", "laser_power_W", "force_N"))
    filter(p -> first(p) ∉ allowed, on.settings) ==
        filter(p -> first(p) ∉ allowed, off.settings) || error("Paired settings differ")
    on.info["initial_states"] == off.info["initial_states"] || error("Paired initial states differ")
    on.info["source_hashes"] == off.info["source_hashes"] || error("Paired sources differ")
    on.time == off.time || error("Paired time grids differ")
    path = first(rows.source_difference_trajectory)
    metadata = TOML.parsefile(joinpath(dirname(path), "comparison_info.toml"))
    metadata["laser_trajectory_sha256"] == on.digest &&
        metadata["no_laser_trajectory_sha256"] == off.digest &&
        metadata["difference_trajectory_sha256"] == filehash(path) ||
        error("Difference provenance mismatch")
    table = Arrow.Table(path)
    compare("difference_time_s", table.time_s, on.time, 0.0)
    for (column, field, tol) in OE_FIELDS
        expected = getproperty(on.oe, field) - getproperty(off.oe, field)
        compare("difference_$column", table[Symbol(column)], expected, tol)
    end
    for (column, values) in on.diagnostics
        compare("difference_$column", table[Symbol(column)],
            values - off.diagnostics[column], 1e-10)
    end
    velocity = [begin
        radial = normalize(off.r[k])
        normal = normalize(cross(off.r[k], off.v[k]))
        transverse = cross(normal, radial)
        delta = on.v[k] - off.v[k]
        (dot(delta, radial), dot(delta, transverse), dot(delta, normal), norm(delta))
    end for k in eachindex(on.time)]
    for (component, label) in enumerate(("R", "T", "N", "norm"))
        field = Symbol("velocity_difference_$(label)_mps")
        compare(string(field), table[field], getindex.(velocity, component), 1e-10)
    end
    for row in eachrow(rows)
        row.run_id == on.info["run_id"] &&
            row.source_no_laser_run_id == off.info["run_id"] || error("Difference endpoint ID mismatch")
        for name in propertynames(table)
            name in propertynames(rows) || continue
            name in (:time_s, :row, :raan_defined) && continue
            tolerance = name == :di_deg ? 2e-6 : name == :da_m ? 1e-5 : 1e-10
            compare("difference_endpoint_$(name)", row[name], table[name][Int(row.row)], tolerance)
        end
    end
end

function check_figure_csvs(root, endpoints)
    for (window, folder) in ((80, "80orbits"), (100, "for_paper"))
        directory = joinpath(root, "images", folder)
        manifest = TOML.parsefile(joinpath(directory, "figure_manifest.toml"))
        for (file, digest) in manifest["figure_checksums"]
            filehash(joinpath(directory, file)) == digest || error("Figure checksum mismatch")
        end
        reference = Dict(endpointkey(row) => row for row in eachrow(endpoints))
        for file in readdir(directory)
            endswith(file, ".csv") && file != "continuous_access_thresholds.csv" || continue
            data = CSV.read(joinpath(directory, file), DataFrame)
            :scheduler in propertynames(data) || continue
            nrow(data) == count(==(window), endpoints.orbits) || error("Incomplete figure CSV")
            allunique(endpointkey.(eachrow(data))) || error("Duplicate figure cases")
            for row in eachrow(data)
                original = reference[endpointkey(row)]
                row.orbits == window || error("Incorrect figure orbit window")
                for field in propertynames(data)
                    isequal(row[field], original[field]) || error("Figure source mismatch: $file $field")
                end
            end
        end
    end
end

function audit()
    tables = summary.((OUTPUT3, OUTPUT4, OUTPUT5))
    nrow.(tables) == (324, 162, 324) || error("Incomplete campaign summaries")
    for table in tables
        allunique(endpointkey.(eachrow(table))) || error("Duplicate endpoint keys")
    end
    groups3 = Dict(casekey(first(eachrow(g))) => g for g in groupby(tables[1],
        [:scheduler, :target_altitude_km, :delta_i_deg, :helpers]))
    groups5 = Dict(casekey(first(eachrow(g))) => g for g in groupby(tables[3],
        [:scheduler, :target_altitude_km, :delta_i_deg, :helpers]))
    for (count, group) in enumerate(groupby(tables[2], [:target_altitude_km, :delta_i_deg, :helpers]))
        baseline = read_track(first(group.source_trajectory))
        check_endpoints(baseline, group)
        for policy in ("naive_next_entering", "positive_along_track")
            id = (policy, key(first(eachrow(group)))...)
            on = read_track(first(groups3[id].source_trajectory))
            check_endpoints(on, groups3[id])
            check_difference(on, baseline, groups5[id])
        end
        println("Audited geometry $count/81: ", key(first(eachrow(group))))
        flush(stdout)
        GC.gc()
    end
    for (root, table) in zip((OUTPUT3, OUTPUT4, OUTPUT5), tables)
        check_figure_csvs(root, table)
    end
    write_metrics("audit"; extra=Dict("trajectory_runs" => 243, "difference_pairs" => 162,
        "samples_per_run" => 1001, "endpoint_rows" => 810))
end

function parameters(info; native=false)
    s = info["settings"]
    Dict{Symbol,Any}(:mu => s["mu_m3_s2"], :c => s["speed_of_light_mps"],
        :N => 2, :masses => fill(s["mass_kg"], 2), :Pmatrix => zeros(2, 2),
        :cavity => Dict{Tuple{Int,Int},Dict{Symbol,Any}}(), :use_los => native,
        :max_range => native ? s["range_m"] : Inf, :min_range => 0.0,
        :use_J2 => true, :J2 => s["J2"], :Re => s["earth_equatorial_radius_m"],
        :use_drag => false, :target_ids => [2], :gve_schedule => "none")
end

function propagate(track; native=false, dtmax=10.0, power_scale=1.0)
    s = track.settings
    s["helper_count"] == 1 || error("The independent propagation test is bounded to one helper")
    p = parameters(track.info; native)
    u0 = vcat([vcat(x["position_m"], x["velocity_mps"]) for x in track.info["initial_states"]]...)
    for (i, state) in enumerate(track.info["initial_states"])
        altitude = i == 1 ? s["helper_altitude_km"] : s["target_altitude_km"]
        inclination = i == 1 ? 0.0 : s["target_inclination_deg"]
        anomaly = i == 1 ? 0.0 : s["target_true_anomaly_deg"]
        independent = collect(Prototype.state_from_OE(s["earth_equatorial_radius_m"] + altitude * 1000;
            i_deg=inclination, ν_deg=anomaly, μ=s["mu_m3_s2"]))
        compare("independent_initial_state", u0[6i-5:6i], independent, 1e-8)
    end
    append!(u0, zeros(4))
    cavity = Dict{Tuple{Int,Int},Dict{Symbol,Any}}((1, 2) =>
        Dict{Symbol,Any}(:B => s["magnification_B"], :Pin => power_scale * s["laser_power_W"]))
    firing = Arrow.Table(joinpath(track.directory, "laser_firing.feather"))
    boundaries = sort(unique(vcat(collect(firing.start_time_s), collect(firing.end_time_s))))
    function update!(p, t)
        active = any(firing.start_time_s[k] <= t < firing.end_time_s[k]
            for k in 1:length(firing.start_time_s))
        p[:cavity] = active ? cavity : empty(cavity)
    end
    native ? (p[:cavity] = cavity) : update!(p, 0.0)
    reference = native && s["scheduler"] == "positive_along_track" ? ForwardPrototype : Prototype
    function rhs!(du, u, p, t)
        reference.nbody_photon!(du, u, p, t)
        forces, _ = reference.laser_forces(u, p)
        force = get(forces, (1, 2), zero(SVector{3,Float64}))
        r, v = SVector{3}(u[7:9]), SVector{3}(u[10:12])
        radial, normal = normalize(r), normalize(cross(r, v))
        transverse = cross(normal, radial)
        du[13:15] .= SVector(dot(force, radial), dot(force, transverse), dot(force, normal)) / p[:masses][2]
        du[16] = dot(force, v)
    end
    events = filter(t -> 0 < t < last(track.time), boundaries)
    callback = native || isempty(events) ? nothing :
        PresetTimeCallback(events, integrator -> update!(integrator.p, integrator.t);
            save_positions=(false, false))
    solution = solve(ODEProblem(rhs!, u0, (0.0, last(track.time)), p), Vern9();
        reltol=1e-12, abstol=1e-12, dtmax, callback, saveat=track.time,
        save_everystep=false, dense=false)
    string(solution.retcode) == "Success" || error("Prototype integration failed: $(solution.retcode)")
    collect(solution.t) == track.time || error("Prototype time grid differs")
    r = [SVector{3}(u[7:9]) for u in solution.u]
    v = [SVector{3}(u[10:12]) for u in solution.u]
    oe = elements(r, v, s["mu_m3_s2"]; equatorial=s["target_inclination_deg"] == 0)
    (; r, v, oe, solution)
end

function check_force_conventions(track)
    s = track.settings
    state = vcat([vcat(x["position_m"], x["velocity_mps"]) for x in track.info["initial_states"]]...)
    p = parameters(track.info; native=true)
    p[:cavity][(1, 2)] = Dict{Symbol,Any}(:B => s["magnification_B"], :Pin => 10_000.0)
    ver2, _ = Prototype.laser_forces(state, p)
    ver3, _ = ForwardPrototype.laser_forces(state, p)
    expected = s["force_N"]
    compare("ver2_BP_over_c_force_N", norm(ver2[(1, 2)]), expected, 1e-15)
    compare("ver3_twice_BP_over_c_force_N", norm(ver3[(1, 2)]), 2expected, 1e-15)
    compare("ver2_recoil_force_N", ver2[(1, 2)], -ver2[(2, 1)], 0.0)
    compare("ver3_recoil_force_N", ver3[(1, 2)], -ver3[(2, 1)], 0.0)
    CSV.write(joinpath(RESULTS, "force_conventions.csv"),
        [(implementation="prototype_ver2", power_W=10_000, force_N=norm(ver2[(1, 2)])),
         (implementation="prototype_ver3", power_W=10_000, force_N=norm(ver3[(1, 2)])),
         (implementation="production_recorded", power_W=10_000, force_N=expected)])
end

function replay()
    rows = NamedTuple[]
    baseline = nothing
    production_baseline = nothing
    matched_force = nothing
    for policy in ("no_laser", "naive_next_entering", "positive_along_track")
        data = summary(policy == "no_laser" ? OUTPUT4 : OUTPUT3)
        selection = filter(row -> row.target_altitude_km == 1050 && row.delta_i_deg == 0.5 &&
            row.helpers == 1 && row.scheduler == policy && row.orbits == 100, data)
        nrow(selection) == 1 || error("Missing reference case")
        track = read_track(only(selection.source_trajectory))
        policy == "no_laser" && (production_baseline = track)
        policy == "naive_next_entering" && check_force_conventions(track)
        modes = policy == "no_laser" ? ("replay",) :
            policy == "naive_next_entering" ? ("replay", "native") :
            ("replay", "native", "native_matched_force", "native_matched_force_dt1")
        for mode in modes
            println("Prototype $policy $mode, 100 orbits"); flush(stdout)
            matched = startswith(mode, "native_matched_force")
            power_scale = matched ? 0.5 : 1.0
            dtmax = endswith(mode, "_dt1") ? 1.0 : 10.0
            result = propagate(track; native=mode != "replay", dtmax, power_scale)
            position_error = maximum(norm.(result.r - track.r))
            velocity_error = maximum(norm.(result.v - track.v))
            sma_error = maximum(abs.(result.oe.a - track.oe.a))
            if mode == "replay"
                compare("replay_position_norm_m", norm.(result.r - track.r), zeros(1001), 0.1)
                compare("replay_velocity_norm_mps", norm.(result.v - track.v), zeros(1001), 1e-4)
                compare("replay_semimajor_axis_m", result.oe.a, track.oe.a, 0.001)
                compare("replay_eccentricity", result.oe.e, track.oe.e, 1e-9)
                compare("replay_inclination_deg", result.oe.i, track.oe.i, 1e-6)
                compare("replay_RAAN_deg", result.oe.raan, track.oe.raan, 1e-6)
                for (j, axis) in enumerate(("R", "T", "N"))
                    compare("replay_laser_dv_$(axis)_mps", getindex.(result.solution.u, 12 + j),
                        track.diagnostics["dv_$(axis)_mps"], 1e-6)
                end
            end
            mode == "native_matched_force" && (matched_force = result)
            if mode == "native_matched_force_dt1"
                compare("matched_force_native_dt_convergence_sma_m", result.oe.a, matched_force.oe.a, 0.001)
            end
            policy == "no_laser" && (baseline = result)
            if mode == "replay"
                compare("replay_Output5_semimajor_axis_m", result.oe.a - baseline.oe.a,
                    track.oe.a - production_baseline.oe.a, 0.001)
            end
            for (orbits, index) in ((80, 801), (100, 1001))
                push!(rows, (; policy, mode, orbits, prototype_power_scale=power_scale, dtmax_s=dtmax,
                    prototype_da_m=result.oe.a[index] - result.oe.a[1],
                    production_da_m=track.oe.a[index] - track.oe.a[1],
                    endpoint_da_error_m=(result.oe.a[index] - result.oe.a[1]) -
                        (track.oe.a[index] - track.oe.a[1]),
                    prototype_on_minus_off_a_m=result.oe.a[index] - baseline.oe.a[index],
                    production_on_minus_off_a_m=track.oe.a[index] - production_baseline.oe.a[index],
                    max_position_error_m=position_error, max_velocity_error_mps=velocity_error,
                    max_sma_error_m=sma_error, source_trajectory=only(selection.source_trajectory)))
            end
            CSV.write(joinpath(RESULTS, "propagation_comparison.csv"), rows)
        end
    end
    write_metrics("replay"; extra=Dict("geometry" => "1 helper at 1000 km; target 1050 km, 0.5 deg, anomaly 1 deg",
        "duration_orbits" => 100, "reference_solver" => "Vern9, rtol=atol=1e-12, dtmax=10 s",
        "replay_scope" => "Original prototype dynamics/forces, matched constants and recorded firing times. Not an independent scheduling test.",
        "native_scope" => "Original ver2 range gating and ver3 positive-T gating at every RHS evaluation; native ver3 has twice production force at the same power.",
        "matched_force_scope" => "Native ver3 run with half input power ONLY in this test, compensating its 2BP/c formula to isolate continuous-vs-accepted-step gating. Source files unchanged. Repeated with dtmax 1 s.",
        "thresholds" => Dict("position_m" => 0.1, "velocity_mps" => 1e-4, "sma_m" => 0.001)))
end

function write_metrics(name; extra=Dict())
    paths = [joinpath(PROTOTYPE, "functions", file) for file in
        ("1_LOS_Metrics.jl", "2_Laser_Forces_ver2.jl", "2_Laser_Forces_ver3.jl",
            "3_Dynamics.jl", "5_OE_Converters.jl")]
    report = merge(Dict("passed" => true, "completed_utc" => string(now(UTC)),
        "prototype_root" => PROTOTYPE, "prototype_hashes" => Dict(path => filehash(path) for path in paths),
        "max_residuals" => copy(RESIDUALS), "comparison_counts" => copy(CHECK_COUNTS)), extra)
    open(joinpath(RESULTS, "$name.toml"), "w") do io
        TOML.print(io, report; sorted=true)
    end
end

function main()
    length(ARGS) == 1 && only(ARGS) in ("audit", "replay") ||
        error("Usage: compare_outputs.jl audit|replay")
    mkpath(RESULTS)
    name = only(ARGS)
    open(joinpath(RESULTS, "$name.toml"), "w") do io
        TOML.print(io, Dict("passed" => false, "status" => "running",
            "started_utc" => string(now(UTC))))
    end
    name == "audit" ? audit() : replay()
    println("$name passed")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
