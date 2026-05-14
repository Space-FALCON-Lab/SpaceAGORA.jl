module AerobrakingPerturbationMC

using Arrow
using CSV
using DataFrames
using Dates
using Distributed
using LinearAlgebra
using Printf
using Statistics
using StaticArrays
using TOML

include(joinpath(@__DIR__, "..", "..", "..", "examples", "common.jl"))

import SpaceAGORA: getDensity, getDensityBatch!

const DEFAULT_PLANETS = [:mars, :earth, :venus]
const DEFAULT_PERIAPSIS_REGIMES = [:shallow, :nominal, :deep]
const DEFAULT_APOAPSIS_REGIMES = [:low, :medium, :high]
const DEFAULT_DYNAMICS_CASES = [
    :two_body,
    :j2,
    :harmonics_low,
    :srp,
    :third_body_sun,
    :gram_aero,
    :full_environment,
]
const DEFAULT_DENSITY_SCALES = Dict(:low => 0.75, :nominal => 1.0, :high => 1.25)

@inline _timestamp() = Dates.format(now(), "HH:MM:SS")

@inline function _sample_field(sample::NamedTuple, key::Symbol, default)
    return haskey(sample, key) ? getfield(sample, key) : default
end

function _sample_label(sample::NamedTuple)::String
    density = _sample_field(sample, :density_case, :none)
    density_part = density === :none ? "" : "/density=$(density)"
    return "$(sample.planet)/$(sample.periapsis_regime)/apo=$(sample.apoapsis_regime)/$(sample.dynamics_case)$(density_part)"
end

function _progress_prefix(sample::NamedTuple)::String
    idx = _sample_field(sample, :sample_index, 0)
    total = _sample_field(sample, :sample_count, 0)
    counter = idx > 0 && total > 0 ? "[$idx/$total] " : ""
    return "[aero-perturb $(_timestamp())] $(counter)worker=$(myid())"
end

function _sample_slug(sample::NamedTuple)::String
    idx = _sample_field(sample, :sample_index, 0)
    prefix = idx > 0 ? @sprintf("case_%03d", idx) : "case"
    density = _sample_field(sample, :density_case, :none)
    return join((
        prefix,
        String(sample.planet),
        String(sample.periapsis_regime),
        "apo_" * String(sample.apoapsis_regime),
        String(sample.dynamics_case),
        "density_" * String(density),
    ), "_")
end

function _log_sample(sample::NamedTuple, message::AbstractString)
    println("$(_progress_prefix(sample)) $message")
    flush(stdout)
    return nothing
end

function _arrow_safe_dataframe(df::DataFrame)::DataFrame
    out = copy(df)
    for name in names(out)
        col = out[!, name]
        if any(value -> value isa Symbol, col)
            out[!, name] = [
                value === missing ? missing : value isa Symbol ? String(value) : value
                for value in col
            ]
        end
    end
    return out
end

function _write_feather(path::String, df::DataFrame)::String
    Arrow.write(path, _arrow_safe_dataframe(df))
    return path
end

function _write_table_pair(base_path_without_ext::String, df::DataFrame)
    csv_path = base_path_without_ext * ".csv"
    feather_path = base_path_without_ext * ".feather"
    CSV.write(csv_path, df)
    _write_feather(feather_path, df)
    return (; csv_path, feather_path)
end

Base.@kwdef struct StudySpec
    planets::Vector{Symbol} = copy(DEFAULT_PLANETS)
    periapsis_regimes::Vector{Symbol} = copy(DEFAULT_PERIAPSIS_REGIMES)
    apoapsis_regimes::Vector{Symbol} = copy(DEFAULT_APOAPSIS_REGIMES)
    dynamics_cases::Vector{Symbol} = copy(DEFAULT_DYNAMICS_CASES)
    density_scales::Dict{Symbol, Float64} = copy(DEFAULT_DENSITY_SCALES)
    norbits::Int = 30
    procs::Int = parse(Int, get(ENV, "SPACEAGORA_AERO_PERTURB_PROCS", "0"))
    output_dir::String = joinpath(REPO_ROOT, "output", "aerobraking_perturbation_mc")
    results::Bool = true
    smoke::Bool = get(ENV, "SPACEAGORA_AERO_PERTURB_SMOKE", "0") == "1"
end

const GRAM_TRAJECTORY_DENSITY_ENV = Dict(
    "SPACEAGORA_GRAM_TRACK_CACHE" => "off",
    "SPACEAGORA_VACUUM_GRAM_CACHE" => "1",
    "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "24",
    "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => "900",
    "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "5000",
)

function configure_gram_trajectory_density!()
    applied = Dict{String, String}()
    for (name, default_value) in GRAM_TRAJECTORY_DENSITY_ENV
        override_name = "SPACEAGORA_AERO_PERTURB_" * replace(name, "SPACEAGORA_" => "")
        value = get(ENV, override_name, default_value)
        ENV[name] = value
        applied[name] = value
    end
    return applied
end

struct DensityScaleModel{M} <: SM.AbstractDensityModel
    base::M
    scale::Float64
end

function Base.getproperty(model::DensityScaleModel, name::Symbol)
    if name === :base || name === :scale
        return getfield(model, name)
    end
    base = getfield(model, :base)
    if hasproperty(base, name)
        return getproperty(base, name)
    end
    return getfield(model, name)
end

SM.SimulationCallbacks._is_gram_density_model(model::DensityScaleModel) =
    SM.SimulationCallbacks._is_gram_density_model(model.base)

SM.SimulationCallbacks._gram_track_trajectory_supported(model::DensityScaleModel) =
    SM.SimulationCallbacks._gram_track_trajectory_supported(model.base)

SM.SimulationCallbacks.density_model_threadsafe(model::DensityScaleModel)::Bool =
    SM.SimulationCallbacks.density_model_threadsafe(model.base)

function getDensity(model::DensityScaleModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p)
    rho, temp, wind_vec = Base.invokelatest(getDensity, model.base, h, lat, lon, el_time, wind, p)
    return model.scale * rho, temp, wind_vec
end

function getDensity(model::DensityScaleModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)
    rho, temp, wind_vec = Base.invokelatest(getDensity, model.base, h, lat, lon, el_time, wind)
    return model.scale * rho, temp, wind_vec
end

function getDensityBatch!(
    rhos,
    temps,
    winds,
    model::DensityScaleModel,
    hs,
    lats,
    lons,
    el_time,
    wind::Bool,
    p,
)
    Base.invokelatest(getDensityBatch!, rhos, temps, winds, model.base, hs, lats, lons, el_time, wind, p)
    @inbounds for i in eachindex(rhos)
        rhos[i] *= model.scale
    end
    return nothing
end

@inline _is_aero_case(case::Symbol)::Bool = case in (:gram_aero, :full_environment)

function _planet(planet_id::Symbol)
    if planet_id == :mars
        return Mars("", SPICE_PATH)
    elseif planet_id == :earth
        return Earth("", SPICE_PATH)
    elseif planet_id == :venus
        return Venus("", SPICE_PATH)
    end
    throw(ArgumentError("Unsupported planet '$planet_id'. Supported planets: mars, earth, venus."))
end

function _planet_name(planet_id::Symbol)::String
    planet_id == :mars && return "mars"
    planet_id == :earth && return "earth"
    planet_id == :venus && return "venus"
    return String(planet_id)
end

function _periapsis_altitude_m(planet_id::Symbol, regime::Symbol)::Float64
    table = if planet_id == :mars
        Dict(:shallow => 150e3, :nominal => 125e3, :deep => 105e3)
    elseif planet_id == :venus
        Dict(:shallow => 180e3, :nominal => 150e3, :deep => 125e3)
    else
        Dict(:shallow => 180e3, :nominal => 145e3, :deep => 120e3)
    end
    return table[regime]
end

function _apoapsis_altitude_m(planet_id::Symbol, regime::Symbol)::Float64
    table = if planet_id == :mars
        Dict(:low => 500e3, :medium => 2_000e3, :high => 4_500e3)
    elseif planet_id == :venus
        Dict(:low => 600e3, :medium => 5_000e3, :high => 20_000e3)
    else
        Dict(:low => 600e3, :medium => 5_000e3, :high => 30_000e3)
    end
    return table[regime]
end

function _harmonics_file(planet_id::Symbol)::String
    if planet_id == :mars
        return joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "Mars50c.csv")
    elseif planet_id == :earth
        return joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
    elseif planet_id == :venus
        return joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "MGNP180U.csv")
    end
    throw(ArgumentError("No harmonics file configured for '$planet_id'."))
end

function _density_model(planet_id::Symbol, dynamics_case::Symbol, density_scale::Float64)
    if !_is_aero_case(dynamics_case)
        return NoAtmosphereModel()
    end
    setup_gram_example!(@__MODULE__)
    base = Base.invokelatest(GRAMAtmosphereModel; planet_name=_planet_name(planet_id))
    return DensityScaleModel(base, density_scale)
end

function _spacecraft(planet, rp_alt_m::Float64, ra_alt_m::Float64)
    return make_three_body_spacecraft(
        bus_dims=(2.2, 2.6, 1.7),
        panel_dims=(0.01, 3.0, 1.5),
        bus_mass=350.0,
        panel_mass_each=10.0,
        panel_offset_y=2.6,
        ic=InitialCondition(
            ra=planet.Rp_e + ra_alt_m,
            rp=planet.Rp_e + rp_alt_m,
            i=93.0,
            ω=80.0,
            Ω=30.0,
            ν=180.0,
        ),
        reflection_coefficient=0.9,
        prop_mass=50.0,
        id=1,
    )
end

function _dynamic_effectors(planet_id::Symbol, planet, spacecraft, dynamics_case::Symbol)
    if dynamics_case == :two_body
        return (InverseSquaredGravityModel(),)
    elseif dynamics_case == :j2
        return (InverseSquaredJ2GravityModel(),)
    elseif dynamics_case == :harmonics_low
        return (GravitationalHarmonicsModel(4, 4, _harmonics_file(planet_id), planet),)
    elseif dynamics_case == :srp
        return (
            InverseSquaredGravityModel(),
            SolarRadiationPressureModel(spacecraft.root.reflection_coefficient, spacecraft.root.ref_area),
        )
    elseif dynamics_case == :third_body_sun
        return (
            InverseSquaredGravityModel(),
            NBodyGravityModel(body_names=("Sun",), primary_body_name=planet.name, planet=planet),
        )
    elseif dynamics_case == :gram_aero
        return (InverseSquaredGravityModel(), AerodynamicCoefficientfM())
    elseif dynamics_case == :full_environment
        return (
            NBodyGravityModel(body_names=("Sun",), primary_body_name=planet.name, planet=planet),
            GravitationalHarmonicsModel(4, 4, _harmonics_file(planet_id), planet),
            SolarRadiationPressureModel(spacecraft.root.reflection_coefficient, spacecraft.root.ref_area),
            AerodynamicCoefficientfM(),
        )
    end
    throw(ArgumentError("Unsupported dynamics case '$dynamics_case'."))
end

@inline function _period_s(planet, rp_alt_m::Float64, ra_alt_m::Float64)::Float64
    a = planet.Rp_e + 0.5 * (rp_alt_m + ra_alt_m)
    return 2π * sqrt(a^3 / planet.μ)
end

function _make_config(sample::NamedTuple, results_directory::String, results::Bool)
    planet = _planet(sample.planet)
    rp_alt_m = _periapsis_altitude_m(sample.planet, sample.periapsis_regime)
    ra_alt_m = _apoapsis_altitude_m(sample.planet, sample.apoapsis_regime)
    spacecraft = _spacecraft(planet, rp_alt_m, ra_alt_m)
    mission_time = sample.norbits * _period_s(planet, rp_alt_m, ra_alt_m)
    dynamic_effectors = _dynamic_effectors(sample.planet, planet, spacecraft, sample.dynamics_case)
    density_model = _density_model(sample.planet, sample.dynamics_case, sample.density_scale)

    return make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=mission_time,
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=dynamic_effectors,
        density_model=density_model,
        orientation_sim=false,
        keplerian=false,
        EI_km=max(220.0, rp_alt_m / 1e3 + 80.0),
        verbose=false,
        results=results,
        results_directory=results_directory,
    )
end

function make_samples(spec::StudySpec)
    samples = NamedTuple[]
    periapsis_regimes = spec.smoke ? spec.periapsis_regimes[1:min(end, 1)] : spec.periapsis_regimes
    apoapsis_regimes = spec.smoke ? spec.apoapsis_regimes[1:min(end, 1)] : spec.apoapsis_regimes
    dynamics_cases = spec.smoke ? [:two_body, :gram_aero] : spec.dynamics_cases
    density_items = sort(collect(spec.density_scales); by=x -> String(x.first))
    for planet in spec.planets
        for periapsis_regime in periapsis_regimes
            for apoapsis_regime in apoapsis_regimes
                for dynamics_case in dynamics_cases
                    density_iter = _is_aero_case(dynamics_case) ? density_items : [(:none => NaN)]
                    for density_pair in density_iter
                        push!(samples, (
                            planet=planet,
                            periapsis_regime=periapsis_regime,
                            apoapsis_regime=apoapsis_regime,
                            dynamics_case=dynamics_case,
                            density_case=Symbol(density_pair.first),
                            density_scale=Float64(density_pair.second),
                            norbits=spec.smoke ? 1 : spec.norbits,
                        ))
                    end
                end
            end
        end
    end
    return samples
end

function _state_pos_vel(u)
    sc = getproperty(u, :sc)[1]
    pos = SVector{3, Float64}(getproperty(sc, :pos))
    vel = SVector{3, Float64}(getproperty(sc, :vel))
    return pos, vel
end

function _elements_from_rv(pos::SVector{3, Float64}, vel::SVector{3, Float64}, μ::Float64)
    r = norm(pos)
    v2 = dot(vel, vel)
    hvec = cross(pos, vel)
    h = norm(hvec)
    nvec = cross(SVector{3, Float64}(0.0, 0.0, 1.0), hvec)
    n = norm(nvec)
    evec = cross(vel, hvec) / μ - pos / r
    e = norm(evec)
    energy = 0.5 * v2 - μ / r
    a = -μ / (2.0 * energy)
    inc = acos(clamp(hvec[3] / h, -1.0, 1.0))
    raan = n > 0.0 ? atan(nvec[2], nvec[1]) : 0.0
    argp = (n > 0.0 && e > 0.0) ? acos(clamp(dot(nvec, evec) / (n * e), -1.0, 1.0)) : 0.0
    if e > 0.0 && evec[3] < 0.0
        argp = 2π - argp
    end
    return (; a, e, inc, raan, argp, energy, rp=a * (1.0 - e), ra=a * (1.0 + e))
end

function _density_metrics(args, sol)
    density_model = args.environment_model.density_model
    if density_model isa NoAtmosphereModel
        return (peak_density=0.0, integrated_density=0.0, max_dynamic_pressure=0.0, integrated_dynamic_pressure=0.0, time_below_interface_s=0.0)
    end

    planet = args.environment_model.planet
    EI_m = args.environment_model.EI * 1e3
    peak_density = 0.0
    integrated_density = 0.0
    max_dynamic_pressure = 0.0
    integrated_dynamic_pressure = 0.0
    time_below_interface_s = 0.0

    for i in eachindex(sol.t)
        pos, vel = _state_pos_vel(sol.u[i])
        r = norm(pos)
        alt = r - planet.Rp_e
        lat = asind(clamp(pos[3] / r, -1.0, 1.0))
        lon = atan(pos[2], pos[1]) * 180.0 / π
        density_state = try
            getDensity(density_model, Float64(alt), Float64(lat), Float64(lon), Float64(sol.t[i]), true)
        catch
            return (peak_density=NaN, integrated_density=NaN, max_dynamic_pressure=NaN, integrated_dynamic_pressure=NaN, time_below_interface_s=time_below_interface_s)
        end
        rho, _, wind = density_state
        q = 0.5 * rho * norm(vel - wind)^2
        peak_density = max(peak_density, rho)
        max_dynamic_pressure = max(max_dynamic_pressure, q)
        if i > 1
            dt = sol.t[i] - sol.t[i - 1]
            integrated_density += rho * dt
            integrated_dynamic_pressure += q * dt
            if alt <= EI_m
                time_below_interface_s += dt
            end
        end
    end

    return (; peak_density, integrated_density, max_dynamic_pressure, integrated_dynamic_pressure, time_below_interface_s)
end

@inline function _state_mass_kg(u)::Float64
    sc = getproperty(u, :sc)[1]
    return Float64(getproperty(sc, :mass))
end

function _trajectory_dataframe_from_solution(sol)
    times = Float64[]
    pos1 = Float64[]
    pos2 = Float64[]
    pos3 = Float64[]
    vel1 = Float64[]
    vel2 = Float64[]
    vel3 = Float64[]
    mass = Float64[]
    for i in eachindex(sol.t)
        pos, vel = _state_pos_vel(sol.u[i])
        push!(times, Float64(sol.t[i]))
        push!(pos1, pos[1]); push!(pos2, pos[2]); push!(pos3, pos[3])
        push!(vel1, vel[1]); push!(vel2, vel[2]); push!(vel3, vel[3])
        push!(mass, _state_mass_kg(sol.u[i]))
    end
    return DataFrame(
        time=times,
        sc1_pos_1=pos1,
        sc1_pos_2=pos2,
        sc1_pos_3=pos3,
        sc1_vel_1=vel1,
        sc1_vel_2=vel2,
        sc1_vel_3=vel3,
        sc1_mass=mass,
    )
end

@inline function _row_float(row, name::Symbol, default::Float64=NaN)::Float64
    hasproperty(row, name) || return default
    value = getproperty(row, name)
    value === missing && return default
    return Float64(value)
end

function _trajectory_dataframe(args, sol, trajectory_csv_path::String)
    if isfile(trajectory_csv_path)
        return CSV.read(trajectory_csv_path, DataFrame)
    end
    return _trajectory_dataframe_from_solution(sol)
end

function _state_sample_from_row(row, spacecraft)
    pos = SVector{3, Float64}(
        _row_float(row, :sc1_pos_1),
        _row_float(row, :sc1_pos_2),
        _row_float(row, :sc1_pos_3),
    )
    vel = SVector{3, Float64}(
        _row_float(row, :sc1_vel_1),
        _row_float(row, :sc1_vel_2),
        _row_float(row, :sc1_vel_3),
    )
    fallback_mass = Float64(spacecraft.dry_mass + spacecraft.prop_mass)
    mass = _row_float(row, :sc1_mass, fallback_mass)
    return SM.StateSample(pos, vel, mass; spacecraft=spacecraft)
end

function _planet_frame_sample(args, state::SM.StateSample, t::Float64)
    planet = args.environment_model.planet
    ephemerides_model = args.environment_model.ephemerides_model
    et0 = SM.ephemerides_time_seconds(args.initial_time, ephemerides_model)
    et = et0 + t
    l_pi = SM.planet_frame_lpi(planet, et, ephemerides_model)
    pos_pp = SVector{3, Float64}(l_pi * state.pos_ii)
    vel_pp = SVector{3, Float64}(l_pi * (state.vel_ii - cross(planet.ω, state.pos_ii)))
    alt, lat, lon = SM.SimulationCallbacks.rtolatlong(pos_pp, planet)
    return SM.PlanetFrameSample(l_pi, pos_pp, vel_pp, alt, lat, lon)
end

function _solar_sample(args, t::Float64)
    planet = args.environment_model.planet
    ephemerides_model = args.environment_model.ephemerides_model
    et0 = SM.ephemerides_time_seconds(args.initial_time, ephemerides_model)
    et = et0 + t
    primary = SM.DynamicEffectors._spice_query_name(planet.name)
    sun_pos = SM.EphemeridesModels.spice_position_j2000_m("sun", et, primary)
    return SM.SolarEphemerisSample(SVector{3, Float64}(sun_pos))
end

function _third_body_sample(args, model::SM.NBodyGravityModel, t::Float64)
    ephemerides_model = args.environment_model.ephemerides_model
    et0 = SM.ephemerides_time_seconds(args.initial_time, ephemerides_model)
    et = et0 + t
    primary = SM.DynamicEffectors._spice_query_name(model.primary_body_name)
    positions = ntuple(length(model.body_names)) do k
        body = SM.DynamicEffectors._spice_query_name(model.body_names[k])
        SVector{3, Float64}(SM.EphemeridesModels.spice_position_j2000_m(body, et, primary))
    end
    return SM.ThirdBodyEphemerisSample(model.body_names, positions)
end

function _aero_force_from_row(row)::SVector{3, Float64}
    force = SVector{3, Float64}(0.0, 0.0, 0.0)
    for prefix in (:sc1_drag, :sc1_lift, :sc1_cross)
        force += SVector{3, Float64}(
            _row_float(row, Symbol(prefix, "_1"), 0.0),
            _row_float(row, Symbol(prefix, "_2"), 0.0),
            _row_float(row, Symbol(prefix, "_3"), 0.0),
        )
    end
    return force
end

function _central_force(args, state::SM.StateSample, env)::SVector{3, Float64}
    force, _ = SM.wrench(InverseSquaredGravityModel(), state, env, 0.0)
    return force
end

function _force_for_effector(args, state::SM.StateSample, row, effector, t::Float64)::SVector{3, Float64}
    planet_frame = nothing
    solar = nothing
    third_bodies = nothing
    if effector isa InverseSquaredJ2GravityModel || effector isa GravitationalHarmonicsModel
        planet_frame = _planet_frame_sample(args, state, t)
    elseif effector isa SolarRadiationPressureModel
        solar = _solar_sample(args, t)
    elseif effector isa NBodyGravityModel
        third_bodies = _third_body_sample(args, effector, t)
    elseif effector isa AerodynamicCoefficientfM
        return _aero_force_from_row(row)
    end
    env = SM.EnvironmentSample(args.environment_model.planet; planet_frame=planet_frame, solar=solar, third_bodies=third_bodies)
    force, _ = SM.wrench(effector, state, env, t)
    if effector isa InverseSquaredJ2GravityModel || effector isa GravitationalHarmonicsModel
        return force - _central_force(args, state, env)
    end
    return force
end

function _active_perturbation_force(args, sample::NamedTuple, state::SM.StateSample, row, t::Float64)::SVector{3, Float64}
    sample.dynamics_case == :two_body && return SVector{3, Float64}(0.0, 0.0, 0.0)

    total = SVector{3, Float64}(0.0, 0.0, 0.0)
    for effector in args.dynamics_model.dynamic_effectors
        if effector isa InverseSquaredGravityModel
            continue
        end
        total += _force_for_effector(args, state, row, effector, t)
    end
    return total
end

function _write_active_force_history!(args, sample::NamedTuple, sol, case_dir::String)
    trajectory_csv_path = joinpath(case_dir, "simulation_results.csv")
    trajectory_feather_path = joinpath(case_dir, "simulation_results.feather")
    trajectory_df = _trajectory_dataframe(args, sol, trajectory_csv_path)
    _write_feather(trajectory_feather_path, trajectory_df)
    spacecraft = args.dynamics_model.spacecraft[1]

    fx = Vector{Float64}(undef, nrow(trajectory_df))
    fy = similar(fx)
    fz = similar(fx)
    fmag = similar(fx)
    for (i, row) in enumerate(eachrow(trajectory_df))
        t = _row_float(row, :time)
        state = _state_sample_from_row(row, spacecraft)
        force = _active_perturbation_force(args, sample, state, row, t)
        fx[i] = force[1]
        fy[i] = force[2]
        fz[i] = force[3]
        fmag[i] = norm(force)
    end

    force_df = DataFrame(
        time=trajectory_df.time,
        active_perturbation_force_ii_1=fx,
        active_perturbation_force_ii_2=fy,
        active_perturbation_force_ii_3=fz,
        active_perturbation_force_mag=fmag,
    )
    force_paths = _write_table_pair(joinpath(case_dir, "active_force_history"), force_df)

    combined_df = copy(trajectory_df)
    combined_df[!, :active_perturbation_force_ii_1] = fx
    combined_df[!, :active_perturbation_force_ii_2] = fy
    combined_df[!, :active_perturbation_force_ii_3] = fz
    combined_df[!, :active_perturbation_force_mag] = fmag
    combined_paths = _write_table_pair(joinpath(case_dir, "trajectory_with_active_force"), combined_df)

    return (
        simulation_results_csv=isfile(trajectory_csv_path) ? trajectory_csv_path : "",
        simulation_results_feather=trajectory_feather_path,
        active_force_history_csv=force_paths.csv_path,
        active_force_history_feather=force_paths.feather_path,
        trajectory_with_active_force_csv=combined_paths.csv_path,
        trajectory_with_active_force_feather=combined_paths.feather_path,
    )
end

function run_sample(sample::NamedTuple; output_dir::String=joinpath(REPO_ROOT, "output"), results::Bool=true)
    total_elapsed = time()
    _log_sample(sample, "START $(_sample_label(sample))")
    case_dir = joinpath(output_dir, _sample_slug(sample))
    mkpath(case_dir)
    args = _make_config(sample, case_dir, results)
    planet = args.environment_model.planet
    rp_alt_m = _periapsis_altitude_m(sample.planet, sample.periapsis_regime)
    ra_alt_m = _apoapsis_altitude_m(sample.planet, sample.apoapsis_regime)
    _log_sample(
        sample,
        @sprintf(
            "CONFIG ready mission_time=%.1f s rp_alt=%.1f km ra_alt=%.1f km",
            args.mission_configuration.mission_time,
            rp_alt_m / 1e3,
            ra_alt_m / 1e3,
        ),
    )
    _log_sample(sample, "SOLVE start")
    elapsed_s = @elapsed sol = run_simulation(args; return_solution=true, isolate_state=false)
    retcode = hasproperty(sol, :retcode) ? String(Symbol(getproperty(sol, :retcode))) : "unknown"
    _log_sample(sample, @sprintf("SOLVE done elapsed=%.2f s retcode=%s saved=%d", elapsed_s, retcode, length(sol.t)))
    first_pos, first_vel = _state_pos_vel(first(sol.u))
    final_pos, final_vel = _state_pos_vel(last(sol.u))
    initial_el = _elements_from_rv(first_pos, first_vel, planet.μ)
    final_el = _elements_from_rv(final_pos, final_vel, planet.μ)
    density = _density_metrics(args, sol)
    _log_sample(sample, "FORCE history start")
    force_paths = _write_active_force_history!(args, sample, sol, case_dir)
    _log_sample(sample, "FORCE history written $(force_paths.active_force_history_csv)")
    _log_sample(
        sample,
        @sprintf(
            "DONE total=%.2f s final_ra_alt=%.3f km final_rp_alt=%.3f km",
            time() - total_elapsed,
            (final_el.ra - planet.Rp_e) / 1e3,
            (final_el.rp - planet.Rp_e) / 1e3,
        ),
    )

    return merge(sample, (
        requested_rp_alt_m=rp_alt_m,
        requested_ra_alt_m=ra_alt_m,
        mission_time_s=args.mission_configuration.mission_time,
        elapsed_s=elapsed_s,
        retcode=retcode,
        success=retcode in ("Success", "Terminated"),
        initial_a_m=initial_el.a,
        initial_e=initial_el.e,
        initial_rp_alt_m=initial_el.rp - planet.Rp_e,
        initial_ra_alt_m=initial_el.ra - planet.Rp_e,
        final_a_m=final_el.a,
        final_e=final_el.e,
        final_rp_alt_m=final_el.rp - planet.Rp_e,
        final_ra_alt_m=final_el.ra - planet.Rp_e,
        final_i_deg=rad2deg(final_el.inc),
        final_raan_deg=rad2deg(final_el.raan),
        final_argp_deg=rad2deg(final_el.argp),
        final_energy_j_kg=final_el.energy,
        delta_a_m=final_el.a - initial_el.a,
        delta_e=final_el.e - initial_el.e,
        delta_rp_alt_m=(final_el.rp - initial_el.rp),
        delta_ra_alt_m=(final_el.ra - initial_el.ra),
        samples_saved=length(sol.t),
        case_dir=case_dir,
    ), density, force_paths)
end

function _sample_failure_row(
    sample::NamedTuple;
    output_dir::String,
    stage::AbstractString,
    err=nothing,
    bt=nothing,
)
    case_dir = joinpath(output_dir, _sample_slug(sample))
    mkpath(case_dir)
    err_text = err === nothing ? "" : sprint(showerror, err, bt)
    failure_path = joinpath(case_dir, "failure.log")
    open(failure_path, "w") do io
        println(io, "stage=$(stage)")
        println(io, "worker=$(myid())")
        println(io, "sample=$sample")
        isempty(err_text) || println(io, err_text)
    end
    _log_sample(sample, "FAILED stage=$(stage) $(failure_path)")
    return merge(sample, (
        requested_rp_alt_m=NaN,
        requested_ra_alt_m=NaN,
        mission_time_s=NaN,
        elapsed_s=NaN,
        retcode="failed",
        success=false,
        failure_stage=String(stage),
        failure_error=err_text,
        failure_log=failure_path,
        initial_a_m=NaN,
        initial_e=NaN,
        initial_rp_alt_m=NaN,
        initial_ra_alt_m=NaN,
        final_a_m=NaN,
        final_e=NaN,
        final_rp_alt_m=NaN,
        final_ra_alt_m=NaN,
        final_i_deg=NaN,
        final_raan_deg=NaN,
        final_argp_deg=NaN,
        final_energy_j_kg=NaN,
        delta_a_m=NaN,
        delta_e=NaN,
        delta_rp_alt_m=NaN,
        delta_ra_alt_m=NaN,
        samples_saved=0,
        case_dir=case_dir,
        peak_density=NaN,
        integrated_density=NaN,
        max_dynamic_pressure=NaN,
        integrated_dynamic_pressure=NaN,
        time_below_interface_s=NaN,
        simulation_results_csv="",
        simulation_results_feather="",
        active_force_history_csv="",
        active_force_history_feather="",
        trajectory_with_active_force_csv="",
        trajectory_with_active_force_feather="",
    ))
end

function run_sample_guarded(sample::NamedTuple; output_dir::String=joinpath(REPO_ROOT, "output"), results::Bool=true)
    try
        row = run_sample(sample; output_dir=output_dir, results=results)
        return merge(row, (failure_stage="", failure_error="", failure_log=""))
    catch err
        return _sample_failure_row(sample; output_dir=output_dir, stage="sample_exception", err=err, bt=catch_backtrace())
    end
end

function _run_indexed_sample_remote(sample::NamedTuple, outdir::String, results::Bool)
    return run_sample_guarded(sample; output_dir=outdir, results=results)
end

function _run_samples_distributed(indexed_samples, outdir::String, results::Bool)
    available = collect(workers())
    isempty(available) && return map(sample -> run_sample_guarded(sample; output_dir=outdir, results=results), indexed_samples)

    rows = Vector{Any}(undef, length(indexed_samples))
    next_sample = Threads.Atomic{Int}(1)
    worker_tasks = Task[]

    for worker in available
        task = @async begin
            while true
                idx = Threads.atomic_add!(next_sample, 1)
                idx > length(indexed_samples) && break
                sample = indexed_samples[idx]
                try
                    rows[idx] = remotecall_fetch(_run_indexed_sample_remote, worker, sample, outdir, results)
                catch err
                    rows[idx] = _sample_failure_row(
                        sample;
                        output_dir=outdir,
                        stage="worker_$(worker)_terminated",
                        err=err,
                        bt=catch_backtrace(),
                    )
                    break
                end
            end
        end
        push!(worker_tasks, task)
    end

    foreach(wait, worker_tasks)
    for idx in eachindex(rows)
        if !isassigned(rows, idx)
            rows[idx] = _sample_failure_row(
                indexed_samples[idx];
                output_dir=outdir,
                stage="not_started_after_worker_exit",
            )
        end
    end
    return rows
end

function paired_deltas(results_df::DataFrame)
    baselines = results_df[
        (results_df.dynamics_case .== :two_body) .&
        (results_df.density_case .== :none),
        :,
    ]
    base_cols = [
        :planet,
        :periapsis_regime,
        :apoapsis_regime,
        :final_a_m,
        :final_e,
        :final_rp_alt_m,
        :final_ra_alt_m,
        :final_i_deg,
        :final_raan_deg,
        :final_argp_deg,
        :final_energy_j_kg,
        :max_dynamic_pressure,
        :integrated_dynamic_pressure,
        :time_below_interface_s,
    ]
    base = baselines[:, base_cols]
    rename!(base, Dict(c => Symbol("baseline_", c) for c in names(base) if !(Symbol(c) in (:planet, :periapsis_regime, :apoapsis_regime))))
    joined = leftjoin(results_df, base; on=[:planet, :periapsis_regime, :apoapsis_regime])
    for col in (:final_a_m, :final_e, :final_rp_alt_m, :final_ra_alt_m, :final_i_deg, :final_raan_deg, :final_argp_deg, :final_energy_j_kg, :max_dynamic_pressure, :integrated_dynamic_pressure, :time_below_interface_s)
        joined[!, Symbol("delta_vs_baseline_", col)] = joined[!, col] .- joined[!, Symbol("baseline_", col)]
    end
    return joined
end

function aggregate_deltas(delta_df::DataFrame)
    grouped = groupby(delta_df, [:planet, :periapsis_regime, :apoapsis_regime, :dynamics_case, :density_case])
    return combine(
        grouped,
        :delta_vs_baseline_final_ra_alt_m => mean => :mean_delta_final_ra_alt_m,
        :delta_vs_baseline_final_rp_alt_m => mean => :mean_delta_final_rp_alt_m,
        :delta_vs_baseline_final_a_m => mean => :mean_delta_final_a_m,
        :delta_vs_baseline_final_e => mean => :mean_delta_final_e,
        :delta_vs_baseline_max_dynamic_pressure => mean => :mean_delta_max_dynamic_pressure,
        :delta_vs_baseline_integrated_dynamic_pressure => mean => :mean_delta_integrated_dynamic_pressure,
        :elapsed_s => mean => :mean_elapsed_s,
        :success => (x -> count(identity, x) / length(x)) => :success_rate,
    )
end

function _manifest(spec::StudySpec, samples, outdir::String)
    gram_density_env = Dict(
        name => get(ENV, name, "")
        for name in sort(collect(keys(GRAM_TRAJECTORY_DENSITY_ENV)))
    )
    return Dict(
        "created_at" => string(now()),
        "repo_root" => REPO_ROOT,
        "julia_version" => string(VERSION),
        "workers" => nworkers(),
        "threads" => Threads.nthreads(),
        "output_dir" => outdir,
        "spec" => Dict(
            "planets" => string.(spec.planets),
            "periapsis_regimes" => string.(spec.periapsis_regimes),
            "apoapsis_regimes" => string.(spec.apoapsis_regimes),
            "dynamics_cases" => string.(spec.dynamics_cases),
            "density_scales" => Dict(String(k) => v for (k, v) in spec.density_scales),
            "norbits" => spec.norbits,
            "save_simulation_csv" => spec.results,
            "save_feather" => true,
            "active_force_history" => true,
            "gram_density_method" => "vacuum_predicted_gram_cache",
            "gram_density_env" => gram_density_env,
            "smoke" => spec.smoke,
            "sample_count" => length(samples),
        ),
    )
end

function run_study(spec::StudySpec=StudySpec())
    gram_density_env = configure_gram_trajectory_density!()
    samples = make_samples(spec)
    indexed_samples = [
        merge(sample, (sample_index=i, sample_count=length(samples)))
        for (i, sample) in enumerate(samples)
    ]
    stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    outdir = joinpath(spec.output_dir, stamp)
    mkpath(outdir)

    println("[aero-perturb] samples=$(length(samples)), norbits=$(spec.smoke ? 1 : spec.norbits), workers=$(nworkers()), output=$(outdir)")
    println("[aero-perturb] GRAM density method=vacuum_predicted_gram_cache env=$(gram_density_env)")
    println("[aero-perturb] progress logging is per simulation: START -> CONFIG -> SOLVE start -> SOLVE done -> DONE")
    flush(stdout)
    rows = if nworkers() > 1
        _run_samples_distributed(indexed_samples, outdir, spec.results)
    else
        map(sample -> run_sample_guarded(sample; output_dir=outdir, results=spec.results), indexed_samples)
    end

    results_df = DataFrame(rows)
    delta_df = paired_deltas(results_df)
    aggregate_df = aggregate_deltas(delta_df)

    _write_table_pair(joinpath(outdir, "results"), results_df)
    _write_table_pair(joinpath(outdir, "paired_deltas"), delta_df)
    _write_table_pair(joinpath(outdir, "aggregates"), aggregate_df)
    open(joinpath(outdir, "manifest.toml"), "w") do io
        TOML.print(io, _manifest(spec, samples, outdir))
    end

    println("[aero-perturb] wrote results/paired_deltas/aggregates as .csv and .feather, plus manifest.toml")
    return (; outdir, results=results_df, deltas=delta_df, aggregates=aggregate_df)
end

function _parse_symbols(raw::String)::Vector{Symbol}
    isempty(strip(raw)) && return Symbol[]
    return Symbol.(strip.(split(raw, ",")))
end

function spec_from_args(args::Vector{String})
    planets = copy(DEFAULT_PLANETS)
    dynamics_cases = copy(DEFAULT_DYNAMICS_CASES)
    procs = parse(Int, get(ENV, "SPACEAGORA_AERO_PERTURB_PROCS", "0"))
    norbits = 30
    output_dir = joinpath(REPO_ROOT, "output", "aerobraking_perturbation_mc")
    smoke = get(ENV, "SPACEAGORA_AERO_PERTURB_SMOKE", "0") == "1"
    results = true

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--planets"
            i += 1; planets = _parse_symbols(args[i])
        elseif arg == "--dynamics"
            i += 1; dynamics_cases = _parse_symbols(args[i])
        elseif arg == "--norbits"
            i += 1; norbits = parse(Int, args[i])
        elseif arg == "--procs"
            i += 1; procs = parse(Int, args[i])
        elseif arg == "--output-dir"
            i += 1; output_dir = abspath(args[i])
        elseif arg == "--smoke"
            smoke = true
        elseif arg == "--save-simulation-csv"
            results = true
        elseif arg == "--no-save-simulation-csv"
            results = false
        elseif arg in ("-h", "--help")
            println("""
            Usage:
              julia --project=. benchmarks/studies/aerobraking_perturbation_mc/main.jl [options]

            Options:
              --planets mars,earth,venus
              --dynamics two_body,j2,harmonics_low,srp,third_body_sun,gram_aero,full_environment
              --norbits N
              --procs N
              --output-dir PATH
              --smoke
              --save-simulation-csv
              --no-save-simulation-csv
            """)
            exit(0)
        else
            throw(ArgumentError("Unknown argument '$arg'. Use --help for options."))
        end
        i += 1
    end

    return StudySpec(
        planets=planets,
        dynamics_cases=dynamics_cases,
        norbits=norbits,
        procs=procs,
        output_dir=output_dir,
        results=results,
        smoke=smoke,
    )
end

end
