if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "9_SpaceAGORA.jl-main_full integration")))
end

module ComparisonSimulation

using Arrow, CSV, DataFrames, Dates, DiffEqBase, LinearAlgebra, OrdinaryDiffEq
using Printf, SHA, SciMLBase, StaticArrays, TOML

include("scheduling_policies.jl")
using .ComparisonPolicies

const LIBRARY_ROOT = normpath(joinpath(@__DIR__, "..", "9_SpaceAGORA.jl-main_full integration"))
const GRAVITY_SOURCE = joinpath(LIBRARY_ROOT, "src", "models", "environment", "gravity",
    "central_acceleration.jl")
include(GRAVITY_SOURCE)

export Settings, propagate, initial_state, orbital_elements, target_table, run_comparison
export ComparisonPolicies, MU, EARTH_RADIUS_M, J2, target_period, duration, settings_dict

const MU = 3.98600436233e14
const EARTH_RADIUS_M = 6_378_136.6
const J2 = 0.00108263
const LIGHT_SPEED = 299_792_458.0
const POLICIES = (:naive_next_entering, :positive_along_track)

Base.@kwdef struct Settings
    target_altitude_km::Float64 = 1000.0
    helper_altitude_km::Float64 = 1050.0
    target_inclination_deg::Float64 = 0.5
    helper_inclination_deg::Float64 = 0.0
    target_anomaly_deg::Float64 = 1.0
    helpers::Int = 200
    orbits::Float64 = 80.0
    mass_kg::Float64 = 227.0
    range_m::Float64 = 200_000.0
    power_W::Float64 = 10_000.0
    magnification::Float64 = 100.0
    beta::Float64 = 1.0
    eta::Float64 = 1.0
    j2::Float64 = J2
    max_step_s::Float64 = 1.0
    reltol::Float64 = 1e-12
    abstol::Float64 = 1e-12
    samples_per_orbit::Int = 200
end

target_period(s::Settings) = 2pi * sqrt((EARTH_RADIUS_M + 1000s.target_altitude_km)^3 / MU)
duration(s::Settings) = s.orbits * target_period(s)
force_magnitude(s::Settings) = s.eta * s.beta * s.magnification * s.power_W / LIGHT_SPEED

function validate(s::Settings)
    all(name -> isfinite(getfield(s, name)), fieldnames(Settings)) ||
        throw(ArgumentError("All settings must be finite"))
    all(>(0), (s.target_altitude_km, s.helper_altitude_km, s.helpers, s.orbits,
        s.mass_kg, s.range_m, s.max_step_s, s.reltol, s.abstol, s.samples_per_orbit)) ||
        throw(ArgumentError("Altitudes, counts, mass, range, duration and solver settings must be positive"))
    all(>=(0), (s.power_W, s.magnification, s.beta, s.eta, s.j2)) ||
        throw(ArgumentError("Laser factors and J2 must be nonnegative"))
    all(i -> 0 <= i <= 180, (s.target_inclination_deg, s.helper_inclination_deg)) ||
        throw(ArgumentError("Inclinations must be in [0, 180] degrees"))
    nothing
end

function initial_state(s::Settings)
    validate(s)
    state = Matrix{Float64}(undef, 6, s.helpers + 1)
    for index in axes(state, 2)
        is_target = index == 1
        radius = EARTH_RADIUS_M + 1000 * (is_target ? s.target_altitude_km : s.helper_altitude_km)
        inclination = deg2rad(is_target ? s.target_inclination_deg : s.helper_inclination_deg)
        anomaly = deg2rad(is_target ? s.target_anomaly_deg : 360 * (index - 2) / s.helpers)
        speed = sqrt(MU / radius)
        state[:, index] .= (radius * cos(anomaly),
            radius * sin(anomaly) * cos(inclination),
            radius * sin(anomaly) * sin(inclination),
            -speed * sin(anomaly),
            speed * cos(anomaly) * cos(inclination),
            speed * cos(anomaly) * sin(inclination))
    end
    state
end

@inline position(u, index) = SVector{3,Float64}(u[1, index], u[2, index], u[3, index])
@inline velocity(u, index) = SVector{3,Float64}(u[4, index], u[5, index], u[6, index])

struct PropagatorParameters{P}
    settings::Settings
    policy::LaserPolicy
    planet::P
end

function dynamics!(du, u, parameters::PropagatorParameters, t)
    @inbounds for index in axes(u, 2)
        acceleration = _inverse_squared_j2_gravity_accel(position(u, index), parameters.planet)
        for axis in 1:3
            du[axis, index] = u[axis + 3, index]
            du[axis + 3, index] = acceleration[axis]
        end
    end
    policy = parameters.policy
    helper = policy.active_helper_idx
    if helper > 0
        acceleration = pair_force(policy, position(u, policy.target_idx), position(u, helper)) /
            parameters.settings.mass_kg
        @inbounds for axis in 1:3
            du[axis + 3, policy.target_idx] += acceleration[axis]
            du[axis + 3, helper] -= acceleration[axis]
        end
    end
    nothing
end

"""Osculating elements in the fixed Earth-equatorial inertial frame."""
function orbital_elements(position, velocity)
    all(isfinite, position) && all(isfinite, velocity) ||
        throw(ArgumentError("Nonfinite Cartesian state"))
    radius = norm(position)
    h = cross(position, velocity)
    radius > 0 && norm(h) > 0 || throw(ArgumentError("Degenerate orbit"))
    energy = dot(velocity, velocity) / 2 - MU / radius
    energy < 0 || throw(ArgumentError("Expected a bound orbit"))
    a = -MU / (2energy)
    e = norm(cross(velocity, h) / MU - position / radius)
    inclination = atan(hypot(h[1], h[2]), h[3])
    raan = hypot(h[1], h[2]) <= 1e-12 * norm(h) ? NaN : atan(h[1], -h[2])
    (a_km=a / 1000, eccentricity=e, inclination_deg=rad2deg(inclination),
        raan_deg=rad2deg(raan))
end

function propagate(s::Settings, schedule::Symbol; progress=false)
    state = initial_state(s)
    policy = LaserPolicy(schedule; helper_indices=collect(2:(s.helpers + 1)),
        range_m=s.range_m, force_N=force_magnitude(s))
    parameters = PropagatorParameters(s, policy, (μ=MU, J2=s.j2, Rp_e=EARTH_RADIUS_M))
    positions = Vector{SVector{3,Float64}}(undef, s.helpers + 1)
    velocities = similar(positions)
    history = (time_s=Float64[], active_helper_idx=Int[],
        distance_m=Union{Missing,Float64}[], along_track_projection=Union{Missing,Float64}[],
        in_range_helpers=Int[])
    last_progress = time()

    function reschedule!(integrator)
        all(isfinite, integrator.u) || error("Nonfinite state at t=$(integrator.t)")
        for index in eachindex(positions)
            positions[index] = position(integrator.u, index)
            velocities[index] = velocity(integrator.u, index)
            norm(positions[index]) > EARTH_RADIUS_M || error("Spacecraft reached the Earth")
        end
        helper = update_schedule!(policy, positions, velocities)
        distance = helper == 0 ? missing : norm(positions[1] - positions[helper])
        projection = helper == 0 ? missing :
            along_track_projection(positions[1], velocities[1], positions[helper])
        helper == 0 || distance <= s.range_m || error("Scheduler selected an out-of-range helper")
        schedule !== :positive_along_track || helper == 0 || projection > 0 ||
            error("Along-track policy selected a nonpositive projection")
        push!(history.time_s, Float64(integrator.t))
        push!(history.active_helper_idx, helper)
        push!(history.distance_m, distance)
        push!(history.along_track_projection, projection)
        push!(history.in_range_helpers, count(policy.in_range))
        # A changed force selection invalidates the integrator's cached derivative.
        DiffEqBase.u_modified!(integrator, true)
        if progress && (time() - last_progress >= 30 || integrator.t == duration(s))
            @printf("%s: %.1f%%, t=%.3f s, scheduler updates=%d\n",
                schedule, 100integrator.t / duration(s), integrator.t, length(history.time_s))
            flush(stdout)
            last_progress = time()
        end
        nothing
    end

    callback = DiffEqBase.DiscreteCallback((u, t, integrator) -> true, reschedule!;
        initialize=(cb, u, t, integrator) -> reschedule!(integrator),
        save_positions=(false, false))
    times = collect(range(0.0, duration(s); length=ceil(Int, s.orbits * s.samples_per_orbit) + 1))
    problem = ODEProblem(dynamics!, state, (0.0, duration(s)), parameters)
    start = time()
    solution = solve(problem, Tsit5(); callback, reltol=s.reltol, abstol=s.abstol,
        dtmax=s.max_step_s, saveat=times, save_everystep=false, dense=false,
        maxiters=10_000_000)
    elapsed_s = time() - start
    SciMLBase.successful_retcode(solution) || error("Solver failed: $(solution.retcode)")
    solution.t[end] == duration(s) || error("Simulation stopped before the requested duration")
    solution.t == times || error("Saved times differ from the requested sampling grid")
    all(u -> all(isfinite, u), solution.u) || error("Saved trajectory contains nonfinite states")
    length(history.time_s) == solution.stats.naccept + 1 ||
        error("Scheduler did not run exactly once per accepted step plus initialization")
    all(>(0), diff(history.time_s)) || error("Non-increasing scheduler times")
    maximum(diff(history.time_s)) <= s.max_step_s + 2eps(duration(s)) ||
        error("Maximum accepted step exceeded its configured limit")
    (; solution, history, policy, settings=s, elapsed_s)
end

function target_table(result)
    solution, settings = result.solution, result.settings
    elements = [orbital_elements(position(u, 1), velocity(u, 1)) for u in solution.u]
    all(oe -> all(isfinite, oe), elements) ||
        error("Target elements are nonfinite; RAAN is undefined for an equatorial target")
    table = DataFrame(elements)
    unwrapped = copy(table.raan_deg)
    for index in 2:length(unwrapped)
        difference = deg2rad(table.raan_deg[index] - table.raan_deg[index - 1])
        unwrapped[index] = unwrapped[index - 1] + rad2deg(atan(sin(difference), cos(difference)))
    end
    table.raan_deg = unwrapped
    insertcols!(table, 1, :time_s => solution.t, :time_days => solution.t ./ 86400,
        :target_orbits => solution.t ./ target_period(settings))
    for name in (:a_km, :eccentricity, :inclination_deg, :raan_deg)
        table[!, Symbol("delta_", name)] = table[!, name] .- table[1, name]
    end
    table
end

function settings_dict(s::Settings)
    settings = Dict{String,Any}(string(name) => getfield(s, name) for name in fieldnames(Settings))
    merge!(settings, Dict("mu_m3_s2" => MU, "earth_equatorial_radius_m" => EARTH_RADIUS_M,
        "speed_of_light_mps" => LIGHT_SPEED, "force_N" => force_magnitude(s),
        "initial_target_period_s" => target_period(s), "duration_s" => duration(s),
        "duration_days" => duration(s) / 86400, "target_index" => 1,
        "helper_indices" => collect(2:(s.helpers + 1)), "initial_eccentricity" => 0.0,
        "initial_RAAN_deg" => 0.0, "initial_argument_of_periapsis_deg" => 0.0,
        "helper_phase_rule_deg" => "360*(helper_index-2)/helpers",
        "frame" => "fixed Earth-equatorial inertial; J2 symmetry axis is +z",
        "epoch_label" => "2026-01-01T00:00:00 (no epoch-dependent forces)",
        "solver" => "Tsit5", "force_application" => "continuous equal/opposite accelerations in RHS",
        "force_range_gate" => "0 < separation <= range_m at every RHS evaluation",
        "scheduler_updates" => "initialization and every accepted solver step",
        "orbit_axis" => "elapsed time / initial Keplerian target period",
        "drag" => false, "solar_radiation_pressure" => false,
        "third_body_gravity" => false, "line_of_sight_occultation_gate" => false))
end

filehash(path) = open(io -> bytes2hex(sha256(io)), path)

function write_info(path, info)
    open(path * ".partial", "w") do io
        TOML.print(io, info; sorted=true)
    end
    mv(path * ".partial", path; force=true)
end

function save_case(s::Settings, schedule::Symbol, root)
    directory = joinpath(root, "results", string(schedule))
    ispath(directory) && error("Refusing to overwrite existing results: $directory")
    mkpath(directory)
    info_path = joinpath(directory, "run_info.toml")
    info = Dict{String,Any}("status" => "running", "policy" => string(schedule),
        "settings" => settings_dict(s), "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(), "started_utc" => string(now(UTC)),
        "source_sha256" => Dict(name => filehash(joinpath(@__DIR__, name))
            for name in ("scheduling_policies.jl", "simulation.jl")),
        "gravity_source" => relpath(GRAVITY_SOURCE, @__DIR__),
        "gravity_source_sha256" => filehash(GRAVITY_SOURCE),
        "project_sha256" => filehash(joinpath(LIBRARY_ROOT, "Project.toml")),
        "manifest_sha256" => filehash(joinpath(LIBRARY_ROOT, "Manifest.toml")))
    write_info(info_path, info)
    try
        println("Running $schedule: $(s.helpers) helpers, $(s.orbits) orbits, J2=$(s.j2)")
        flush(stdout)
        result = propagate(s, schedule; progress=true)
        table = target_table(result)
        CSV.write(joinpath(directory, "target_orbital_elements.csv"), table)
        CSV.write(joinpath(directory, "schedule.csv"), result.history)
        trajectory = DataFrame(time_s=result.solution.t,
            target_orbits=result.solution.t ./ target_period(s))
        for index in 1:(s.helpers + 1), (field, offset) in (("pos", 0), ("vel", 3)), axis in 1:3
            trajectory[!, Symbol("sc$(index)_$(field)_$axis")] =
                [u[offset + axis, index] for u in result.solution.u]
        end
        Arrow.write(joinpath(directory, "trajectory.feather"), trajectory; file=true, compress=:zstd)
        saved = Arrow.Table(joinpath(directory, "trajectory.feather"))
        saved.time_s == result.solution.t || error("Saved trajectory times failed verification")
        length(propertynames(saved)) == 2 + 6 * (s.helpers + 1) ||
            error("Saved trajectory has an incorrect column count")
        for name in propertynames(saved)
            all(isfinite, saved[name]) || error("Saved trajectory has nonfinite values in $name")
            first(saved[name]) == first(trajectory[!, name]) &&
                last(saved[name]) == last(trajectory[!, name]) ||
                error("Saved trajectory endpoints differ in $name")
        end
        selected_time = sum((result.history.active_helper_idx[index] > 0 ?
            result.history.time_s[index + 1] - result.history.time_s[index] : 0.0)
            for index in 1:(length(result.history.time_s) - 1))
        merge!(info, Dict("status" => "verified", "finished_utc" => string(now(UTC)),
            "solver_retcode" => string(result.solution.retcode),
            "accepted_steps" => result.solution.stats.naccept,
            "rejected_steps" => result.solution.stats.nreject,
            "scheduler_updates" => length(result.history.time_s),
            "link_activations" => result.policy.activations,
            "selected_link_time_s" => selected_time,
            "selected_link_duty_fraction" => selected_time / duration(s),
            "selected_link_time_note" => "Selection duration; RHS range gating can shorten actual firing",
            "trajectory_rows" => nrow(trajectory), "trajectory_columns" => ncol(trajectory),
            "wall_propagation_s" => result.elapsed_s,
            "initial_elements" => Dict(string(k) => table[1, k]
                for k in (:a_km, :eccentricity, :inclination_deg, :raan_deg)),
            "final_elements" => Dict(string(k) => table[end, k]
                for k in (:a_km, :eccentricity, :inclination_deg, :raan_deg)),
            "element_changes" => Dict(string(k) => table[end, k]
                for k in (:delta_a_km, :delta_eccentricity, :delta_inclination_deg, :delta_raan_deg)),
            "checksums" => Dict(name => filehash(joinpath(directory, name))
                for name in ("trajectory.feather", "target_orbital_elements.csv", "schedule.csv"))))
        write_info(info_path, info)
        println("VERIFIED $directory")
        flush(stdout)
        info
    catch err
        info["status"] = "failed"
        info["error"] = sprint(showerror, err)
        write_info(info_path, info)
        @error "Policy comparison failed; incomplete files retained" directory exception=(err, catch_backtrace())
        rethrow()
    end
end

function run_comparison(; settings=Settings(), root=@__DIR__)
    validate(settings)
    for policy in POLICIES
        directory = joinpath(root, "results", string(policy))
        ispath(directory) && error("Results already exist; choose a new output root: $directory")
    end
    [save_case(settings, policy, root) for policy in POLICIES]
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) || error("Run without arguments, or call run_comparison(settings=..., root=...) from Julia")
    ComparisonSimulation.run_comparison()
end
