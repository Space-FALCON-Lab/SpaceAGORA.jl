const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using Dates
using Printf
using LinearAlgebra
using StaticArrays
using SPICE
using CSV
using DataFrames
using Statistics

include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
using .SimulationModel

const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :SimulationEngine)
    include(joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
end
if !isdefined(@__MODULE__, :run_simulation)
    const run_simulation = SimulationEngine.run_simulation
end

const EM = SimulationModel.EnvironmentModels
const CB = SimulationModel.SimulationCallbacks

@inline function _pctl(x::Vector{Float64}, p::Float64)
    isempty(x) && return 0.0
    xs = sort(x)
    idx = clamp(Int(ceil(p * length(xs))), 1, length(xs))
    return xs[idx]
end

function parse_cli(args::Vector{String})
    opts = Dict{String, String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || error("Unsupported argument '$arg'. Use --key=value.")
        body = arg[3:end]
        if occursin("=", body)
            key, value = split(body, "=", limit=2)
            opts[key] = value
        else
            key = body
            if i == length(args) || startswith(args[i + 1], "--")
                opts[key] = "true"
            else
                i += 1
                opts[key] = args[i]
            end
        end
        i += 1
    end
    return opts
end

function _with_env(vars::Dict{String, String}, f::Function)
    old = Dict{String, Union{Nothing, String}}()
    for (k, v) in vars
        old[k] = haskey(ENV, k) ? ENV[k] : nothing
        ENV[k] = v
    end
    try
        return f()
    finally
        for (k, ov) in old
            if ov === nothing
                delete!(ENV, k)
            else
                ENV[k] = ov
            end
        end
    end
end

struct ScenarioDef
    name::String
    rp_alt_m::Float64
    ra_alt_m::Float64
    i_deg::Float64
    Ω_deg::Float64
    ω_deg::Float64
    ν_deg::Float64
    mission_time_s::Float64
end

function _scenario_defs(planet_name::String)
    key = lowercase(strip(planet_name))
    if key == "venus"
        return [
            ScenarioDef("drag_passage", 180e3, 430e3, 25.0, 20.0, 30.0, -20.0, 1200.0),
            ScenarioDef("entry", 160e3, 320e3, 22.0, 35.0, 10.0, -35.0, 700.0),
            ScenarioDef("orbit", 430e3, 430e3, 30.0, 15.0, 0.0, 0.0, 1800.0)
        ]
    elseif key == "earth"
        return [
            ScenarioDef("drag_passage", 120e3, 1000e3, 25.0, 20.0, 30.0, -20.0, 1200.0),
            ScenarioDef("entry", 120e3, 360e3, 22.0, 35.0, 10.0, -35.0, 700.0),
            ScenarioDef("orbit", 900e3, 900e3, 30.0, 15.0, 0.0, 0.0, 1800.0)
        ]
    elseif key == "titan"
        return [
            ScenarioDef("drag_passage", 320e3, 2200e3, 25.0, 20.0, 30.0, -20.0, 1200.0),
            ScenarioDef("entry", 280e3, 800e3, 22.0, 35.0, 10.0, -35.0, 700.0),
            ScenarioDef("orbit", 1500e3, 1500e3, 30.0, 15.0, 0.0, 0.0, 1800.0)
        ]
    end

    return [
        ScenarioDef("drag_passage", 110e3, 4500e3, 25.0, 20.0, 30.0, -20.0, 1800.0),
        ScenarioDef("entry", 80e3, 1500e3, 22.0, 35.0, 10.0, -45.0, 1200.0),
        ScenarioDef("orbit", 450e3, 450e3, 30.0, 15.0, 0.0, 0.0, 3600.0)
    ]
end

function _parse_scenarios(raw::String)
    txt = lowercase(strip(raw))
    if isempty(txt) || txt == "all"
        return Set(["drag_passage", "entry", "orbit"])
    end
    chosen = Set{String}()
    for token in split(txt, ",")
        v = lowercase(strip(token))
        v in ("drag_passage", "entry", "orbit") || error("Unsupported scenario '$v'.")
        push!(chosen, v)
    end
    return chosen
end

function _planet_from_name(name::String, spice_path::String)
    key = lowercase(strip(name))
    if key == "earth"
        return Earth("", spice_path)
    elseif key == "mars"
        return Mars("", spice_path)
    elseif key == "venus"
        return Venus("", spice_path)
    elseif key == "titan"
        return Titan("", spice_path)
    elseif key == "jupiter"
        return Jupiter("", spice_path)
    elseif key == "uranus"
        return Uranus("", spice_path)
    elseif key == "neptune"
        return Neptune("", spice_path)
    end
    throw(ArgumentError("Unsupported planet '$name'."))
end

function _make_args(planet_name::String, scenario::ScenarioDef, initial_time::InitialTime)::SimulationConfiguration
    spice_path = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
    planet = _planet_from_name(planet_name, spice_path)

    ic = InitialCondition(
        ra=planet.Rp_e + scenario.ra_alt_m,
        rp=planet.Rp_e + scenario.rp_alt_m,
        i=scenario.i_deg,
        ω=scenario.ω_deg,
        Ω=scenario.Ω_deg,
        ν=scenario.ν_deg
    )

    root = Link(root=true, m=391.0, ref_area=2.2 * 1.7)
    spacecraft = SpacecraftModel(
        joints=Joint[],
        links=Link[root],
        root=root,
        instant_actuation=true,
        prop_mass=0.0,
        inertia_tensor=root.inertia,
        n_reaction_wheels=0,
        n_thrusters=0,
        initial_condition=ic,
        id=100
    )

    density_model = Base.invokelatest(
        GRAMAtmosphereModel;
        planet_name=planet_name,
        initial_time=initial_time
    )

    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
            normalize=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=scenario.mission_time_s,
            orientation_sim=false,
            num_steps_to_save=1000
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=125.0,
            density_model=density_model,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=true
        ),
        dynamics_model=DynamicsModel([spacecraft], (InverseSquaredGravityModel(), AerodynamicCoefficientfM())),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=initial_time,
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            dt_max_orbit=1.0,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_atmosphere=0.2
        )
    )
end

function _run_case(planet_name::String, scenario::ScenarioDef, envvars::Dict{String, String}, initial_time::InitialTime)
    args = _make_args(planet_name, scenario, initial_time)
    elapsed_s = 0.0
    n_steps = 0
    retcode = ""
    sol = nothing
    _with_env(envvars, () -> begin
        SimulationModel.SimulationCallbacks._gram_runtime_stats_reset!()
        elapsed_s = @elapsed begin
            sol = run_simulation(args; isolate_state=false, return_solution=true)
            n_steps = length(sol.t)
            retcode = string(sol.retcode)
        end
    end)
    return (elapsed_s=elapsed_s, n_steps=n_steps, retcode=retcode, sol=sol, args=args)
end

function _extract_track(sol, planet)
    n = length(sol.t)
    t0 = Float64(sol.t[1])
    t = Vector{Float64}(undef, n)
    alt = Vector{Float64}(undef, n)
    lat = Vector{Float64}(undef, n)
    lon = Vector{Float64}(undef, n)
    @inbounds for k in 1:n
        t[k] = Float64(sol.t[k]) - t0
        pos = SVector{3, Float64}(sol.u[k].sc[1].pos)
        vel = SVector{3, Float64}(sol.u[k].sc[1].vel)
        pos_pp, _ = CB.r_intor_p!(pos, vel, planet)
        alt[k], lat[k], lon[k] = CB.rtolatlong(pos_pp, planet)
    end
    return (t=t, alt=alt, lat=lat, lon=lon)
end

function _evaluate_track_density(
    track,
    planet_name::String,
    scenario::ScenarioDef,
    envvars::Dict{String, String},
    initial_time::InitialTime
)
    args = _make_args(planet_name, scenario, initial_time)
    p = ODEParams(n_sats=1, args=args)
    density_model = args.environment_model.density_model
    EM.clear_gram_static_grid_cache!()
    EM.clear_gram_offline_surrogate_cache!()

    n = length(track.t)
    rho = Vector{Float64}(undef, n)
    temp = Vector{Float64}(undef, n)
    wind = Vector{SVector{3, Float64}}(undef, n)
    elapsed_s = 0.0

    _with_env(envvars, () -> begin
        elapsed_s = @elapsed begin
            @inbounds for k in 1:n
                rho[k], temp[k], wind[k] = Base.invokelatest(
                    getDensity,
                    density_model,
                    track.alt[k],
                    track.lat[k],
                    track.lon[k],
                    track.t[k],
                    true,
                    p
                )
            end
        end
    end)

    return (rho=rho, temp=temp, wind=wind, elapsed_s=elapsed_s)
end

function _summarize_errors(scenario::String, track, ref, eval)
    n = length(track.t)
    rho_abs = [abs(eval.rho[k] - ref.rho[k]) for k in 1:n]
    rho_rel = [rho_abs[k] / max(abs(ref.rho[k]), 1e-12) for k in 1:n]
    temp_abs = [abs(eval.temp[k] - ref.temp[k]) for k in 1:n]
    wind_ew_abs = [abs(eval.wind[k][1] - ref.wind[k][1]) for k in 1:n]
    wind_ns_abs = [abs(eval.wind[k][2] - ref.wind[k][2]) for k in 1:n]
    wind_up_abs = [abs(eval.wind[k][3] - ref.wind[k][3]) for k in 1:n]
    wind_abs = [norm(eval.wind[k] - ref.wind[k]) for k in 1:n]

    return (
        scenario=scenario,
        samples=n,
        t_span_s=track.t[end] - track.t[1],
        alt_min_km=minimum(track.alt) * 1e-3,
        alt_max_km=maximum(track.alt) * 1e-3,
        rho_abs_max=maximum(rho_abs),
        rho_abs_p95=_pctl(rho_abs, 0.95),
        rho_rel_max=maximum(rho_rel),
        rho_rel_p95=_pctl(rho_rel, 0.95),
        temp_abs_max=maximum(temp_abs),
        temp_abs_p95=_pctl(temp_abs, 0.95),
        wind_abs_max=maximum(wind_abs),
        wind_abs_p95=_pctl(wind_abs, 0.95),
        wind_ew_abs_max=maximum(wind_ew_abs),
        wind_ns_abs_max=maximum(wind_ns_abs),
        wind_up_abs_max=maximum(wind_up_abs),
        point_eval_s=ref.elapsed_s,
        surrogate_eval_s=eval.elapsed_s
    )
end

function _error_table(scenario::String, track, ref, eval)
    n = length(track.t)
    rho_abs = [abs(eval.rho[k] - ref.rho[k]) for k in 1:n]
    rho_rel = [rho_abs[k] / max(abs(ref.rho[k]), 1e-12) for k in 1:n]
    temp_abs = [abs(eval.temp[k] - ref.temp[k]) for k in 1:n]
    wind_abs = [norm(eval.wind[k] - ref.wind[k]) for k in 1:n]
    wind_ew_abs = [abs(eval.wind[k][1] - ref.wind[k][1]) for k in 1:n]
    wind_ns_abs = [abs(eval.wind[k][2] - ref.wind[k][2]) for k in 1:n]
    wind_up_abs = [abs(eval.wind[k][3] - ref.wind[k][3]) for k in 1:n]

    return DataFrame(
        scenario=fill(scenario, n),
        t_s=track.t,
        alt_m=track.alt,
        lat_deg=rad2deg.(track.lat),
        lon_deg=rad2deg.(track.lon),
        rho_ref=ref.rho,
        rho_eval=eval.rho,
        rho_abs_err=rho_abs,
        rho_rel_err=rho_rel,
        temp_ref=ref.temp,
        temp_eval=eval.temp,
        temp_abs_err=temp_abs,
        wind_ref_ew=[w[1] for w in ref.wind],
        wind_ref_ns=[w[2] for w in ref.wind],
        wind_ref_up=[w[3] for w in ref.wind],
        wind_eval_ew=[w[1] for w in eval.wind],
        wind_eval_ns=[w[2] for w in eval.wind],
        wind_eval_up=[w[3] for w in eval.wind],
        wind_abs_err=wind_abs,
        wind_ew_abs_err=wind_ew_abs,
        wind_ns_abs_err=wind_ns_abs,
        wind_up_abs_err=wind_up_abs
    )
end

function run_matrix()
    opts = parse_cli(copy(ARGS))

    planet_name = lowercase(get(opts, "planet", "mars"))
    selected = _parse_scenarios(get(opts, "scenarios", "all"))
    out_summary = get(opts, "out-summary", joinpath(REPO_ROOT, "output", "gram_real_sim_surrogate_matrix_summary.csv"))
    out_errors = get(opts, "out-errors", joinpath(REPO_ROOT, "output", "gram_real_sim_surrogate_matrix_errors.csv"))
    out_runtime = get(opts, "out-runtime", joinpath(REPO_ROOT, "output", "gram_real_sim_surrogate_matrix_runtime.csv"))

    surrogate_dir = strip(get(opts, "surrogate-dir", ""))

    wind_mode = planet_name == "earth" ? "nominal" : "perturbed"

    point_env = Dict(
        "SPACEAGORA_GRAM_TRACK_CACHE" => "off",
        "SPACEAGORA_GRAM_STATIC_GRID" => "off",
        "SPACEAGORA_GRAM_OFFLINE_SURROGATE" => "off",
        "SPACEAGORA_GRAM_WIND_MODE" => wind_mode,
        "SPACEAGORA_GRAM_PROFILE" => "1"
    )
    surrogate_env = Dict(
        "SPACEAGORA_GRAM_TRACK_CACHE" => "off",
        "SPACEAGORA_GRAM_STATIC_GRID" => "off",
        "SPACEAGORA_GRAM_OFFLINE_SURROGATE" => "on",
        "SPACEAGORA_GRAM_WIND_MODE" => wind_mode,
        "SPACEAGORA_GRAM_PROFILE" => "1"
    )
    if !isempty(surrogate_dir)
        surrogate_env["SPACEAGORA_GRAM_OFFLINE_SURROGATE_DIR"] = surrogate_dir
    end

    default_month = planet_name == "earth" ? 1 : 11
    initial_time = InitialTime(year=2001, month=default_month, day=6, hour=19, minute=0, second=32.0)

    rows_runtime = NamedTuple[]
    rows_summary = NamedTuple[]
    error_tables = DataFrame[]
    surrogate_dir_label = isempty(surrogate_dir) ? "<default per-planet path>" : surrogate_dir

    println("GRAM real simulation matrix: point_to_point vs offline_surrogate")
    println("planet=$planet_name")
    println("surrogate_dir=$surrogate_dir_label")

    for sc in _scenario_defs(planet_name)
        sc.name in selected || continue
        println("\nScenario: $(sc.name)")
        point_run = nothing
        surrogate_run = nothing

        try
            point_run = _run_case(planet_name, sc, point_env, initial_time)
            push!(rows_runtime, (
                scenario=sc.name,
                mode="point_to_point",
                elapsed_s=point_run.elapsed_s,
                n_steps=point_run.n_steps,
                retcode=point_run.retcode
            ))
        catch err
            @warn "Point-to-point run failed; preserving partial scenario results." scenario=sc.name planet=planet_name error=sprint(showerror, err)
            push!(rows_runtime, (
                scenario=sc.name,
                mode="point_to_point",
                elapsed_s=NaN,
                n_steps=0,
                retcode="ERROR"
            ))
        end

        try
            surrogate_run = _run_case(planet_name, sc, surrogate_env, initial_time)
            push!(rows_runtime, (
                scenario=sc.name,
                mode="offline_surrogate",
                elapsed_s=surrogate_run.elapsed_s,
                n_steps=surrogate_run.n_steps,
                retcode=surrogate_run.retcode
            ))
        catch err
            @warn "Offline surrogate run failed; preserving partial scenario results." scenario=sc.name planet=planet_name error=sprint(showerror, err)
            push!(rows_runtime, (
                scenario=sc.name,
                mode="offline_surrogate",
                elapsed_s=NaN,
                n_steps=0,
                retcode="ERROR"
            ))
        end

        if point_run === nothing || surrogate_run === nothing
            @warn "Skipping cross-mode error summary because one mode failed." scenario=sc.name planet=planet_name point_ok=(point_run !== nothing) surrogate_ok=(surrogate_run !== nothing)
            continue
        end

        try
            track = _extract_track(point_run.sol, point_run.args.environment_model.planet)
            point_eval = _evaluate_track_density(track, planet_name, sc, point_env, initial_time)
            surrogate_eval = _evaluate_track_density(track, planet_name, sc, surrogate_env, initial_time)

            err_summary = _summarize_errors(sc.name, track, point_eval, surrogate_eval)
            push!(rows_summary, merge(err_summary, (
                runtime_point_s=point_run.elapsed_s,
                runtime_surrogate_s=surrogate_run.elapsed_s,
                runtime_speedup=point_run.elapsed_s / max(surrogate_run.elapsed_s, 1e-9),
                runtime_point_steps=point_run.n_steps,
                runtime_surrogate_steps=surrogate_run.n_steps,
                runtime_point_retcode=point_run.retcode,
                runtime_surrogate_retcode=surrogate_run.retcode
            )))
            push!(error_tables, _error_table(sc.name, track, point_eval, surrogate_eval))
        catch err
            @warn "Cross-mode comparison failed after both runs completed." scenario=sc.name planet=planet_name error=sprint(showerror, err)
        end
    end

    runtime_df = DataFrame(rows_runtime)
    summary_df = DataFrame(rows_summary)
    errors_df = isempty(error_tables) ? DataFrame() : vcat(error_tables...)

    mkpath(dirname(out_runtime))
    mkpath(dirname(out_summary))
    mkpath(dirname(out_errors))
    CSV.write(out_runtime, runtime_df)
    CSV.write(out_summary, summary_df)
    CSV.write(out_errors, errors_df)

    println("\nRuntime:")
    show(runtime_df, allrows=true, allcols=true)
    println("\n\nSummary:")
    show(summary_df, allrows=true, allcols=true)
    println("\n\nSaved CSVs:")
    println("  runtime: $out_runtime")
    println("  summary: $out_summary")
    println("  errors : $out_errors")

    return (runtime=runtime_df, summary=summary_df, errors=errors_df)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_matrix()
end
