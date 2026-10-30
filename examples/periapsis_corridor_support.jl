module PeriapsisCorridorSupport

using CSV
using DataFrames
using LinearAlgebra
using Roots
using SHA
using StaticArrays
using TOML
using SpaceAGORA

export CorridorConfig, load_corridor_config, corridor_atmosphere_backend, CORRIDOR_ATMOSPHERE_BACKENDS, fitted_atmosphere_parameters,
    CorridorAtmosphereFit, corridor_density, corridor_scale_height, atmosphere_fit_relative_errors,
    corridor_coordinate, corridor_eccentricity_bounds, dimensionless_context,
    geodetic_altitude_latitude, geodetic_radius_m, geodetic_altitude_at_radius, corridor_radius_band,
    physical_apsis_events, write_csv, write_manifest, radial_velocity, drag_context_ratio

struct CorridorConfig
    path::String
    values::Dict{String, Any}
end

const CORRIDOR_ATMOSPHERE_BACKENDS = ("us76_piecewise", "gram")

"""Selected atmosphere backend; a config without `backend` is the US76 fit."""
corridor_atmosphere_backend(config::CorridorConfig) =
    String(get(config.values["atmosphere"], "backend", "us76_piecewise"))

function load_corridor_config(path::AbstractString=joinpath(@__DIR__, "periapsis_corridor_config.toml"))
    parsed = TOML.parsefile(path)
    mission = parsed["mission"]
    lo = Float64(mission["corridor_min_altitude_m"])
    mid = Float64(mission["reference_altitude_m"])
    hi = Float64(mission["corridor_max_altitude_m"])
    lo < mid < hi || throw(ArgumentError("periapsis corridor must satisfy min < reference < max"))
    isfinite(Float64(mission["horizon_s"])) && Float64(mission["horizon_s"]) > 0.0 ||
        throw(ArgumentError("the configured horizon must be finite and positive"))
    vehicle = parsed["vehicle"]
    for key in ("dry_mass_kg", "propellant_mass_kg", "ballistic_coefficient_kg_m2",
                "drag_coefficient", "thrust_n", "specific_impulse_s")
        value = Float64(vehicle[key])
        isfinite(value) && value > 0.0 || throw(ArgumentError("vehicle.$key must be finite and positive"))
    end
    integration = parsed["integration"]
    for key in ("reltol", "abstol", "dt_max_orbit_s", "dt_max_atmosphere_s")
        value = Float64(integration[key])
        isfinite(value) && value > 0.0 || throw(ArgumentError("integration.$key must be finite and positive"))
    end
    calibration = parsed["calibration"]
    for key in ("tolerance_m", "horizon_fraction_of_period")
        value = Float64(calibration[key])
        isfinite(value) && value > 0.0 || throw(ArgumentError("calibration.$key must be finite and positive"))
    end
    Int(calibration["max_iterations"]) >= 1 || throw(ArgumentError("calibration.max_iterations must be at least 1"))
    backend = String(get(parsed["atmosphere"], "backend", "us76_piecewise"))
    backend in CORRIDOR_ATMOSPHERE_BACKENDS || throw(ArgumentError(
        "atmosphere.backend must be one of $(join(CORRIDOR_ATMOSPHERE_BACKENDS, ", ")); got \"$backend\""))
    if backend == "gram"
        gram = get(parsed, "atmosphere_gram", nothing)
        gram isa AbstractDict && haskey(gram, "planet_name") && haskey(gram, "entry_interface_km") ||
            throw(ArgumentError("backend \"gram\" requires [atmosphere_gram] planet_name and entry_interface_km"))
        ei = Float64(gram["entry_interface_km"])
        isfinite(ei) && ei > 0.0 || throw(ArgumentError("atmosphere_gram.entry_interface_km must be finite and positive"))
    end
    comparison = parsed["comparison"]
    for key in ("planning_lead_deg", "max_burn_delta_v_mps", "max_burn_duration_s",
                "total_delta_v_budget_mps", "target_periapsis_altitude_m")
        value = Float64(comparison[key])
        isfinite(value) && value > 0.0 || throw(ArgumentError("comparison.$key must be finite and positive"))
    end
    return CorridorConfig(abspath(path), parsed)
end

"""
Piecewise exponential density fit. Layer `k` spans `node_altitudes_m[k]` to
`node_altitudes_m[k+1]`, has density `node_density_kg_m3[k]` at its lower node,
and scale height `scale_heights_m[k]`; density is continuous at every node.
`below`/`above` (`:extend` or `:zero`) set the behavior outside the node range.
"""
struct CorridorAtmosphereFit
    node_altitudes_m::Vector{Float64}
    node_density_kg_m3::Vector{Float64}
    scale_heights_m::Vector{Float64}
    below::Symbol
    above::Symbol
    temperature_k::Float64
    source::String
    table_altitudes_m::Vector{Float64}
    table_density_kg_m3::Vector{Float64}
end

@inline function _layer(fit::CorridorAtmosphereFit, h::Float64)
    return clamp(searchsortedlast(fit.node_altitudes_m, h), 1, length(fit.scale_heights_m))
end

"""Fitted density (kg/m3) at geodetic altitude `h_m` above the planet ellipsoid."""
function corridor_density(fit::CorridorAtmosphereFit, h_m::Real)
    h = Float64(h_m)
    h < fit.node_altitudes_m[1] && fit.below === :zero && return 0.0
    h > fit.node_altitudes_m[end] && fit.above === :zero && return 0.0
    k = _layer(fit, h)
    return fit.node_density_kg_m3[k] * exp(-(h - fit.node_altitudes_m[k]) / fit.scale_heights_m[k])
end

"""Scale height (m) of the layer that contains (or, outside the nodes, is extended to) `h_m`."""
corridor_scale_height(fit::CorridorAtmosphereFit, h_m::Real) = fit.scale_heights_m[_layer(fit, Float64(h_m))]

function fitted_atmosphere_parameters(config::CorridorConfig)
    atm = config.values["atmosphere"]
    table_h = Float64.(atm["table_altitudes_m"])
    table_rho = Float64.(atm["table_density_kg_m3"])
    nodes = Float64.(atm["node_altitudes_m"])
    length(table_h) == length(table_rho) && issorted(table_h; lt=<=) && all(>(0.0), table_rho) ||
        throw(ArgumentError("atmosphere table must have matching, increasing altitudes and positive densities"))
    length(nodes) >= 2 && issorted(nodes; lt=<=) ||
        throw(ArgumentError("atmosphere needs at least two increasing node altitudes"))
    node_rho = map(nodes) do h
        idx = findfirst(==(h), table_h)
        idx === nothing && throw(ArgumentError("node altitude $h m is not a table altitude"))
        table_rho[idx]
    end
    all(diff(node_rho) .< 0.0) || throw(ArgumentError("node densities must decrease with altitude"))
    scale_heights = diff(nodes) ./ log.(node_rho[1:end-1] ./ node_rho[2:end])
    behavior(key) = begin
        value = Symbol(atm[key])
        value in (:extend, :zero) || throw(ArgumentError("atmosphere.$key must be \"extend\" or \"zero\""))
        value
    end
    return CorridorAtmosphereFit(nodes, node_rho, scale_heights,
        behavior("below_bottom_node"), behavior("above_top_node"),
        Float64(atm["temperature_k"]), String(atm["source"]), table_h, table_rho)
end

"""Relative error `(fit - table) / table` at every table row."""
function atmosphere_fit_relative_errors(fit::CorridorAtmosphereFit)
    return [(corridor_density(fit, h) - rho) / rho
            for (h, rho) in zip(fit.table_altitudes_m, fit.table_density_kg_m3)]
end

function corridor_coordinate(altitude_m::Real, minimum_m::Real, maximum_m::Real)
    lo, hi = Float64(minimum_m), Float64(maximum_m)
    lo < hi || throw(ArgumentError("corridor minimum must be below its maximum"))
    return (Float64(altitude_m) - lo) / (hi - lo)
end

# Geodetic altitude is the repo's own: `rtolatlong` (Bowring closed form; the
# routine the engine feeds to GRAM and to the drag effector's planet-frame
# sample) and the ellipsoid position of `CartesianInitialCondition`'s geodetic
# constructor. Ellipsoid: flattening (Rp_e - Rp_p)/Rp_e from the planet model.
const _rtolatlong = SpaceAGORA.SimulationModel.SimulationCallbacks.rtolatlong

# Ellipsoid position of a geodetic (latitude, longitude, altitude); the same
# formula as the engine's geodetic initial-condition constructor, kept local
# because that helper is not available as a stable entry point.
function _geodetic_position(lat_rad::Float64, lon_rad::Float64, alt_m::Float64, planet)
    f = (planet.Rp_e - planet.Rp_p) / planet.Rp_e
    e2 = 1.0 - (1.0 - f)^2
    n = planet.Rp_e / sqrt(1.0 - e2 * sin(lat_rad)^2)
    return SVector{3, Float64}(
        (n + alt_m) * cos(lat_rad) * cos(lon_rad),
        (n + alt_m) * cos(lat_rad) * sin(lon_rad),
        ((1.0 - e2) * n + alt_m) * sin(lat_rad)
    )
end

"""`(geodetic altitude m, geodetic latitude rad)` of a planet-fixed position `r_pp_m`."""
function geodetic_altitude_latitude(r_pp_m, planet)
    llh = _rtolatlong(SVector{3, Float64}(r_pp_m), planet)
    return (Float64(llh[1]), Float64(llh[2]))
end

"""Geocentric radius (m) of the point at geodetic `altitude_m` and geodetic latitude `latitude_rad`."""
geodetic_radius_m(altitude_m::Real, latitude_rad::Real, planet) =
    norm(_geodetic_position(Float64(latitude_rad), 0.0, Float64(altitude_m), planet))

"""Geodetic altitude (m) of the point at geocentric `radius_m` along the planet-fixed direction of `r_pp_m`."""
function geodetic_altitude_at_radius(r_pp_m, radius_m::Real, planet)
    r = SVector{3, Float64}(r_pp_m)
    return geodetic_altitude_latitude(r * (Float64(radius_m) / norm(r)), planet)[1]
end

"""
Radius band `(rp_min, rp_max)` that the geodetic corridor `[minimum, maximum]`
occupies at geodetic latitude `latitude_rad`. The band is latitude dependent:
at 50 degrees the ellipsoid surface is about 12.6 km inside `Rp_e`.
"""
corridor_radius_band(latitude_rad::Real, minimum_altitude_m::Real, maximum_altitude_m::Real, planet) =
    (geodetic_radius_m(minimum_altitude_m, latitude_rad, planet),
     geodetic_radius_m(maximum_altitude_m, latitude_rad, planet))

"""Eccentricity band that keeps periapsis radius inside `[rp_min_radius_m, rp_max_radius_m]` at semi-major axis `a_m`."""
function corridor_eccentricity_bounds(a_m::Real, rp_min_radius_m::Real, rp_max_radius_m::Real)
    a, rp_min, rp_max = Float64(a_m), Float64(rp_min_radius_m), Float64(rp_max_radius_m)
    a > 0.0 && rp_min < rp_max || throw(ArgumentError("invalid semi-major axis or corridor"))
    a > rp_min || throw(ArgumentError("the periapsis corridor does not admit an elliptic orbit at this semi-major axis"))
    e_min = max(0.0, 1.0 - rp_max / a)
    e_max = min(prevfloat(1.0), 1.0 - rp_min / a)
    e_min < e_max || throw(ArgumentError("the periapsis corridor has no noncircular elliptic eccentricity band"))
    return (minimum_eccentricity=e_min, maximum_eccentricity=e_max)
end

function dimensionless_context(a_m::Real, e::Real, inclination_rad::Real,
                               j2::Real, j3::Real, gravity_radius_m::Real,
                               reference_radius_m::Real, rp_min_radius_m::Real,
                               rp_max_radius_m::Real)
    a, ecc = Float64(a_m), Float64(e)
    0.0 <= ecc < 1.0 || throw(ArgumentError("context requires an elliptic orbit"))
    p = a * (1.0 - ecc^2)
    rg = Float64(gravity_radius_m)
    pi_j2 = Float64(j2) * (rg / p)^2
    pi_j3_signed = Float64(j3) * (rg / p)^3
    # The radius band comes from the geodetic corridor at the event latitude
    # (`corridor_radius_band`); `reference_radius_m` is only the Rp_e/a scale.
    ebounds = corridor_eccentricity_bounds(a, rp_min_radius_m, rp_max_radius_m)
    return (p_m=p, periapsis_radius_m=a * (1.0 - ecc),
            periapsis_radius_above_reference_m=a * (1.0 - ecc) - Float64(reference_radius_m),
            reference_radius_over_a=Float64(reference_radius_m) / a,
            Pi_J2=pi_j2, Pi_J3=abs(pi_j3_signed),
            signed_Pi_J3=pi_j3_signed, J3_sign=sign(Float64(j3)),
            Pi_3B=0.0, Gamma_3B_J2=0.0,
            eccentricity_corridor=ebounds,
            corridor_eccentricity_width=ebounds.maximum_eccentricity - ebounds.minimum_eccentricity,
            inclination_rad=Float64(inclination_rad), drag_context=0.0,
            third_body_context=0.0, low_thrust_context=0.0)
end

@inline function drag_context_ratio(density_kg_m3::Real, relative_speed_mps::Real,
                                    ballistic_coefficient_kg_m2::Real,
                                    mu_m3_s2::Real, radius_m::Real)
    rho, speed, beta, mu, radius = Float64.((density_kg_m3, relative_speed_mps,
                                              ballistic_coefficient_kg_m2, mu_m3_s2, radius_m))
    rho >= 0.0 && speed >= 0.0 && beta > 0.0 && mu > 0.0 && radius > 0.0 ||
        throw(ArgumentError("drag context requires nonnegative density/speed and positive beta, mu, and radius"))
    drag_accel = 0.5 * rho * speed^2 / beta
    central_accel = mu / radius^2
    return (Pi_drag=drag_accel / central_accel,
            dynamic_pressure_pa=0.5 * rho * speed^2,
            drag_acceleration_mps2=drag_accel,
            central_acceleration_mps2=central_accel)
end

@inline function _state_rv(u)
    sc = hasproperty(u, :sc) ? u.sc[1] : u[1]
    mass = hasproperty(sc, :mass) ? Float64(sc.mass) : NaN
    return SVector{3, Float64}(sc.pos), SVector{3, Float64}(sc.vel), mass
end

@inline function radial_velocity(u)::Float64
    r, v, _ = _state_rv(u)
    return dot(r, v)
end

"""
Locate all physical-radius apsides by root-finding radial velocity (the
authority) on the dense solution. `geodetic_of(t, r_ii)` returns
`(geodetic altitude m, geodetic latitude rad)` at the root; only that
conversion defines the reported altitude.
"""
function physical_apsis_events(solution, geodetic_of; root_atol_s::Float64=1e-7)
    length(solution.t) >= 2 || throw(ArgumentError("solution needs at least two stored integration times"))
    events = NamedTuple[]
    f0 = radial_velocity(solution.u[1])
    for idx in 1:(length(solution.t) - 1)
        t0, t1 = Float64(solution.t[idx]), Float64(solution.t[idx + 1])
        f1 = radial_velocity(solution.u[idx + 1])
        # Periapsis: r·v changes negative to positive; apoapsis: positive to negative.
        if (f0 < 0.0 && f1 >= 0.0) || (f0 > 0.0 && f1 <= 0.0)
            dense_f0 = radial_velocity(solution(t0))
            dense_f1 = radial_velocity(solution(t1))
            r0, v0, _ = _state_rv(solution(t0))
            r1, v1, _ = _state_rv(solution(t1))
            tol0 = 1e-12 * max(norm(r0) * norm(v0), 1.0)
            tol1 = 1e-12 * max(norm(r1) * norm(v1), 1.0)
            root = if abs(dense_f0) <= tol0
                t0
            elseif abs(dense_f1) <= tol1
                t1
            elseif signbit(dense_f0) != signbit(dense_f1)
                find_zero(t -> radial_velocity(solution(t)), (t0, t1), Roots.Brent(); xatol=root_atol_s)
            elseif min(abs(f0), abs(f1)) <= 1e-9 * max(norm(r0) * norm(v0), 1.0)
                # The solver steps exactly onto the event time, so the saved
                # knot holds the root (|r.v| at roundoff level) while the dense
                # interpolant at that knot is off by interpolation error and can
                # carry the wrong sign. Accept the saved endpoint nearest zero.
                abs(f0) <= abs(f1) ? t0 : t1
            else
                # No dense sign change and no saved endpoint at the root:
                # do not manufacture a root.
                f0 = f1
                continue
            end
            r, v, mass = _state_rv(solution(root))
            kind = f0 < f1 ? "periapsis" : "apoapsis"
            if isempty(events) || root - events[end].time_s > 1e-3
                altitude_m, latitude_rad = geodetic_of(Float64(root), r)
                push!(events, (time_s=Float64(root), kind=kind,
                               physical_altitude_m=altitude_m, latitude_rad=latitude_rad,
                               radius_m=norm(r), position_ii_m=r, velocity_ii_mps=v,
                               mass_kg=mass,
                               radial_velocity_residual_m2_s=dot(r, v)))
            end
        end
        f0 = f1
    end
    return events
end

function write_csv(path::AbstractString, records)
    mkpath(dirname(path))
    df = DataFrame(records)
    CSV.write(path, df)
    return df
end

function write_manifest(path::AbstractString, config::CorridorConfig;
                        case_name::AbstractString, gravity_path::AbstractString,
                        solver_settings=Dict{String, Any}(), initial_calibration=NamedTuple[],
                        actual_initial_state=nothing, source_files=String[],
                        run_status::AbstractString="not_started",
                        atmosphere_provenance=Dict{String, Any}())
    gravity_file = abspath(gravity_path)
    manifest = Dict{String, Any}(
        "case" => String(case_name),
        "run_status" => String(run_status),
        "configuration_file" => config.path,
        "configuration_sha256" => bytes2hex(sha256(read(config.path))),
        "effective_configuration" => config.values,
        "gravity_file" => gravity_file,
        "gravity_file_sha256" => bytes2hex(sha256(read(gravity_file))),
        "solver_settings" => Dict{String, Any}(solver_settings),
        "atmosphere_backend" => corridor_atmosphere_backend(config),
        "atmosphere_fit_used_for_drag" => corridor_atmosphere_backend(config) == "us76_piecewise",
        "density_source" => corridor_atmosphere_backend(config) == "gram" ?
            "Earth-GRAM mean density via the sampled environment atmosphere (not the fit)" :
            "piecewise exponential fit to the 1976 U.S. Standard Atmosphere",
        "atmosphere_provenance" => Dict{String, Any}(atmosphere_provenance),
        "atmosphere_source" => String(config.values["atmosphere"]["source"]),
        "mean_conversion" => String(config.values["gravity"]["conversion_model"]),
        "altitude_reference" => String(config.values["mission"]["altitude_reference"]),
        "horizon_s" => Float64(config.values["mission"]["horizon_s"]),
        "corridor_min_altitude_m" => Float64(config.values["mission"]["corridor_min_altitude_m"]),
        "corridor_max_altitude_m" => Float64(config.values["mission"]["corridor_max_altitude_m"]),
        "beta_kg_m2" => Float64(config.values["vehicle"]["ballistic_coefficient_kg_m2"]),
        "ballistic_coefficient_is_assumed_not_measured" => true,
        "drag_coefficient" => Float64(config.values["vehicle"]["drag_coefficient"]),
        "constant_CdA_m2" => (Float64(config.values["vehicle"]["dry_mass_kg"]) +
            Float64(config.values["vehicle"]["propellant_mass_kg"])) /
            Float64(config.values["vehicle"]["ballistic_coefficient_kg_m2"]),
        "equivalent_reference_area_m2" => ((Float64(config.values["vehicle"]["dry_mass_kg"]) +
            Float64(config.values["vehicle"]["propellant_mass_kg"])) /
            Float64(config.values["vehicle"]["ballistic_coefficient_kg_m2"])) /
            Float64(config.values["vehicle"]["drag_coefficient"]),
        "drag_effector" => "constant-CdA ballistic drag proxy; attitude and placeholder geometry do not set drag area",
        "atmosphere_scale_heights_m" => fitted_atmosphere_parameters(config).scale_heights_m,
        "truth_includes_j3" => true,
        "claims_safety_from_mean_only" => false,
        "initial_calibration" => [Dict{String, Any}(string(key) => value for (key, value) in pairs(record))
                                  for record in initial_calibration],
        "source_sha256" => Dict{String, Any}(
            basename(abspath(source)) => bytes2hex(sha256(read(source)))
            for source in source_files if isfile(source)),
    )
    actual_initial_state === nothing || (manifest["actual_initial_state"] =
        Dict{String, Any}(string(key) => value for (key, value) in pairs(actual_initial_state)))
    project_root = normpath(joinpath(@__DIR__, ".."))
    manifest["git_revision"] = try
        strip(read(`git -C $project_root rev-parse HEAD`, String))
    catch
        "unavailable"
    end
    manifest["local_diff_summary"] = try
        strip(read(`git -C $project_root diff --stat`, String))
    catch
        "unavailable"
    end
    open(path, "w") do io
        TOML.print(io, manifest)
    end
    return manifest
end

end # module
