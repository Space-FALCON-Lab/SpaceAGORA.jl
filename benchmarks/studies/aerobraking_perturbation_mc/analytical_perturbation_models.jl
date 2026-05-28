module AnalyticalPerturbationModels

using Dates
using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "..", "..", "..", "examples", "common.jl"))

using .SimulationModel

const SM = SimulationModel

export analytical_series, analytical_labels, simulated_model_ratio_series

const PLANET_ORDER = ("mars", "earth", "venus", "titan")
const PLANET_K_HIGHER_HARMONICS = Dict(
    "mars" => 13.0e-5,
    "venus" => 1.2e-5,
    "earth" => 1.0e-5,
    "titan" => 1.0e-5,
)
const BASIC_HARMONICS_KAULA_TRUNCATION_FLOOR = 1.0e-9
const MODEL_KEYS = (:j2, :harmonics, :srp, :third_body, :drag, :full_environment)
const ANALYTICAL_SPICE_LOCK = ReentrantLock()

const MODEL_LABELS = Dict(
    :j2 => "J2",
    :harmonics => "Harmonics",
    :srp => "SRP",
    :third_body => "3rd Body",
    :drag => "Aero",
    :full_environment => "Full Env",
)

const _INITIAL_TIME = InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
const _EPHEMERIDES_MODEL = SM.SpiceEphemeridesModel()
const _MU_SYMBOL = Symbol("μ")
const _OMEGA_SYMBOL = Symbol("ω")

struct AnalyticalPlanetCache
    planet
    harmonics_model
    degree2_harmonics_model
    j2_model
    srp_model
    nbody_model
    et0_s::Float64
end

struct CaseContext
    planet_id::String
    planet
    mass_scale::Float64
    rp_m::Float64
    ra_m::Float64
    ra_hat::SVector{3, Float64}
    semi_major_m::Float64
    total_ref_area_m2::Float64
    nominal_mass_kg::Float64
    beta_kg_m2::Float64
    harmonics_model
    degree2_harmonics_model
    j2_model
    srp_model
    nbody_model
    et0_s::Float64
end

const _PLANET_CACHE = Dict{String, AnalyticalPlanetCache}()

@inline _planet_mu(planet)::Float64 = Float64(getproperty(planet, _MU_SYMBOL))
@inline _planet_omega(planet)::SVector{3, Float64} = SVector{3, Float64}(getproperty(planet, _OMEGA_SYMBOL))

function analytical_labels()
    return copy(MODEL_LABELS)
end

function _planet(planet_id::String)
    planet_id == "mars" && return Mars("", SPICE_PATH)
    planet_id == "earth" && return Earth("", SPICE_PATH)
    planet_id == "venus" && return Venus("", SPICE_PATH)
    planet_id == "titan" && return Titan("", SPICE_PATH)
    throw(ArgumentError("Unsupported planet '$planet_id'."))
end

function _planet_cache(planet_id::String)::AnalyticalPlanetCache
    cached = get(_PLANET_CACHE, planet_id, nothing)
    cached !== nothing && return cached

    planet = _planet(planet_id)
    degree, order = _harmonics_degree_order(planet_id)
    ref_area = _spacecraft_reference_area_m2()
    harmonics = GravitationalHarmonicsModel(
        degree,
        order,
        _harmonics_file(planet_id),
        planet;
        include_central=false,
    )
    degree2_harmonics = GravitationalHarmonicsModel(
        2,
        min(2, order),
        _harmonics_file(planet_id),
        planet;
        include_central=false,
    )
    nbody = NBodyGravityModel(
        body_names=_third_body_names(planet_id),
        primary_body_name=planet.name,
        planet=planet,
    )
    cached = AnalyticalPlanetCache(
        planet,
        harmonics,
        degree2_harmonics,
        InverseSquaredJ2GravityModel(),
        SolarRadiationPressureModel(0.9, ref_area),
        nbody,
        SM.ephemerides_time_seconds(_INITIAL_TIME, _EPHEMERIDES_MODEL),
    )
    _PLANET_CACHE[planet_id] = cached
    return cached
end

function _periapsis_altitude_m(planet_id::String, regime::String)::Float64
    table = if planet_id == "mars"
        Dict("shallow" => 150e3, "nominal" => 125e3, "deep" => 110e3)
    elseif planet_id == "venus"
        Dict("shallow" => 180e3, "nominal" => 150e3, "deep" => 135e3)
    elseif planet_id == "titan"
        Dict("shallow" => 900e3, "nominal" => 720e3, "deep" => 650e3)
    else
        Dict("shallow" => 180e3, "nominal" => 145e3, "deep" => 120e3)
    end
    return table[regime]
end

function _harmonics_file(planet_id::String)::String
    planet_id == "mars" && return joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "Mars50c.csv")
    planet_id == "earth" && return joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
    planet_id == "venus" && return joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "MGNP180U.csv")
    planet_id == "titan" && return joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "titan5.csv")
    throw(ArgumentError("No harmonics file configured for '$planet_id'."))
end

function _harmonics_degree_order(planet_id::String)::Tuple{Int, Int}
    planet_id == "titan" && return (5, 5)
    return (20, 20)
end

function _third_body_names(planet_id::String)
    planet_id == "earth" && return ("Sun", "Moon")
    planet_id == "titan" && return ("Saturn", "Sun")
    return ("Sun",)
end

function _spacecraft_reference_area_m2()::Float64
    bus_dims = (2.2, 2.6, 1.7)
    panel_dims = (0.01, 3.0, 1.5)
    return bus_dims[1] * bus_dims[3] + 2.0 * panel_dims[2] * panel_dims[3]
end

function _nominal_mass_kg(mass_scale::Float64)::Float64
    return (350.0 + 2.0 * 10.0 + 50.0) * mass_scale
end

function _drag_beta_from_history(tbl, fallback_mass::Float64, fallback_area::Float64)::Float64
    names_set = Set(Symbol.(propertynames(tbl)))
    required = (:sc1_drag_1, :sc1_drag_2, :sc1_drag_3, :dynamic_pressure_pa, :sc1_mass)
    if !all(name -> name in names_set, required)
        return fallback_mass / fallback_area
    end

    best_idx = 0
    best_q = -Inf
    for i in eachindex(tbl.dynamic_pressure_pa)
        q = Float64(tbl.dynamic_pressure_pa[i])
        if isfinite(q) && q > best_q
            best_q = q
            best_idx = i
        end
    end
    best_idx == 0 && return fallback_mass / fallback_area
    best_q > 0.0 || return fallback_mass / fallback_area

    drag = norm(SVector{3, Float64}(
        Float64(tbl.sc1_drag_1[best_idx]),
        Float64(tbl.sc1_drag_2[best_idx]),
        Float64(tbl.sc1_drag_3[best_idx]),
    ))
    mass = Float64(tbl.sc1_mass[best_idx])
    if isfinite(drag) && drag > 0.0 && isfinite(mass) && mass > 0.0
        cd_area = drag / best_q
        if isfinite(cd_area) && cd_area > 0.0
            return mass / cd_area
        end
    end
    return fallback_mass / fallback_area
end

function _apoapsis_direction(tbl)::SVector{3, Float64}
    names_set = Set(Symbol.(propertynames(tbl)))
    all(name -> name in names_set, (:sc1_pos_1, :sc1_pos_2, :sc1_pos_3)) || return SVector{3, Float64}(1.0, 0.0, 0.0)

    best_idx = 0
    best_r2 = -Inf
    for i in eachindex(tbl.sc1_pos_1)
        pos = SVector{3, Float64}(
            Float64(tbl.sc1_pos_1[i]),
            Float64(tbl.sc1_pos_2[i]),
            Float64(tbl.sc1_pos_3[i]),
        )
        r2 = dot(pos, pos)
        if isfinite(r2) && r2 > best_r2
            best_r2 = r2
            best_idx = i
        end
    end
    best_idx == 0 && return SVector{3, Float64}(1.0, 0.0, 0.0)

    ra_vec = SVector{3, Float64}(
        Float64(tbl.sc1_pos_1[best_idx]),
        Float64(tbl.sc1_pos_2[best_idx]),
        Float64(tbl.sc1_pos_3[best_idx]),
    )
    ra_norm = norm(ra_vec)
    ra_norm > 0.0 && isfinite(ra_norm) || return SVector{3, Float64}(1.0, 0.0, 0.0)
    return ra_vec / ra_norm
end

function _context(info, tbl)::CaseContext
    planet_id = String(info.planet)
    cached = _planet_cache(planet_id)
    planet = cached.planet
    rp_alt_m = _periapsis_altitude_m(planet_id, String(info.periapsis))
    ra_alt_m = Float64(info.apoapsis_alt_km) * 1e3
    rp_m = planet.Rp_e + rp_alt_m
    ra_m = planet.Rp_e + ra_alt_m
    ra_hat = _apoapsis_direction(tbl)
    a_m = 0.5 * (rp_m + ra_m)
    mass_scale = Float64(info.spacecraft_mass_scale)
    ref_area = _spacecraft_reference_area_m2()
    nominal_mass = _nominal_mass_kg(mass_scale)
    beta = _drag_beta_from_history(tbl, nominal_mass, ref_area)
    return CaseContext(
        planet_id,
        planet,
        mass_scale,
        rp_m,
        ra_m,
        ra_hat,
        a_m,
        ref_area,
        nominal_mass,
        beta,
        cached.harmonics_model,
        cached.degree2_harmonics_model,
        cached.j2_model,
        cached.srp_model,
        cached.nbody_model,
        cached.et0_s,
    )
end

function _state(tbl, idx::Int)
    pos = SVector{3, Float64}(
        Float64(tbl.sc1_pos_1[idx]),
        Float64(tbl.sc1_pos_2[idx]),
        Float64(tbl.sc1_pos_3[idx]),
    )
    vel = SVector{3, Float64}(
        Float64(tbl.sc1_vel_1[idx]),
        Float64(tbl.sc1_vel_2[idx]),
        Float64(tbl.sc1_vel_3[idx]),
    )
    mass = Float64(tbl.sc1_mass[idx])
    return pos, vel, mass
end

function _planet_frame(ctx::CaseContext, pos::SVector{3, Float64}, vel::SVector{3, Float64}, t::Float64)
    et = ctx.et0_s + t
    l_pi = SM.planet_frame_lpi(ctx.planet, et, _EPHEMERIDES_MODEL)
    pos_pp = SVector{3, Float64}(l_pi * pos)
    vel_pp = SVector{3, Float64}(l_pi * (vel - cross(_planet_omega(ctx.planet), pos)))
    r = norm(pos_pp)
    alt = r - ctx.planet.Rp_e
    lat = asin(clamp(pos_pp[3] / r, -1.0, 1.0))
    lon = atan(pos_pp[2], pos_pp[1])
    return SM.PlanetFrameSample(l_pi, pos_pp, vel_pp, alt, lat, lon)
end

function _third_body_positions(ctx::CaseContext, t::Float64)
    et = ctx.et0_s + t
    primary = SM.DynamicEffectors._spice_query_name(ctx.nbody_model.primary_body_name)
    positions = ntuple(length(ctx.nbody_model.body_names)) do k
        body = SM.DynamicEffectors._spice_query_name(ctx.nbody_model.body_names[k])
        SVector{3, Float64}(SM.EphemeridesModels.spice_position_j2000_m(body, et, primary))
    end
    return positions
end

function _sun_position(ctx::CaseContext, t::Float64)
    et = ctx.et0_s + t
    primary = SM.DynamicEffectors._spice_query_name(ctx.planet.name)
    return SVector{3, Float64}(SM.EphemeridesModels.spice_position_j2000_m("sun", et, primary))
end

@inline function _central_ratio_denominator(ctx::CaseContext, r::Float64)::Float64
    return _planet_mu(ctx.planet) / (r * r)
end

function _basic_j2(ctx::CaseContext)::Float64
    return 1.5 * ctx.planet.J2 * (ctx.planet.Rp_e / ctx.semi_major_m)^2
end

function _kaula_harmonics_truncation_degree(ctx::CaseContext)::Int
    max_degree = min(ctx.harmonics_model.L, size(ctx.harmonics_model.C, 1) - 1)
    max_degree < 3 && return max_degree
    K = PLANET_K_HIGHER_HARMONICS[ctx.planet_id]
    radius_ratio = ctx.planet.Rp_e / ctx.semi_major_m
    truncation_degree = 3
    for degree in 3:max_degree
        kaula_term = K / degree^2 * radius_ratio^degree
        kaula_term < BASIC_HARMONICS_KAULA_TRUNCATION_FLOOR && break
        truncation_degree = degree
    end
    return truncation_degree
end

function _basic_harmonics(ctx::CaseContext)::Float64
    total = 0.0
    radius_ratio = ctx.planet.Rp_e / ctx.semi_major_m
    for degree in 3:_kaula_harmonics_truncation_degree(ctx)
        total += abs(ctx.harmonics_model.C[degree + 1, 1]) * radius_ratio^degree
    end
    return total
end

function _basic_srp(ctx::CaseContext)::Float64
    sun_pos = _sun_position(ctx, 0.0)
    solar_pressure_1au_pa = 4.56e-6
    au_m = ctx.srp_model.AU_m
    sun_distance_m = norm(sun_pos)
    accel = solar_pressure_1au_pa * (au_m / sun_distance_m)^2 *
        ctx.srp_model.Cr * ctx.srp_model.A / ctx.nominal_mass_kg
    return accel / _central_ratio_denominator(ctx, ctx.semi_major_m)
end

function _basic_third_body(ctx::CaseContext)::Float64
    positions = _third_body_positions(ctx, 0.0)
    total = 0.0
    for (mu3, pos3) in zip(ctx.nbody_model.body_mus, positions)
        r3 = norm(pos3)
        r3 > 0.0 && isfinite(r3) || continue
        alpha = dot(pos3 / r3, ctx.ra_hat)
        total += mu3 * ctx.ra_m^2 / (r3^3 * _planet_mu(ctx.planet)) * sqrt(1.0 + 3.0 * alpha^2)
    end
    return total
end

function _saved_aero_force_ratio(ctx::CaseContext, tbl, idx::Int)::Float64
    names_set = Set(Symbol.(propertynames(tbl)))
    pos, _, mass = _state(tbl, idx)
    force = if all(name -> name in names_set, (:sc1_drag_1, :sc1_drag_2, :sc1_drag_3, :sc1_lift_1, :sc1_lift_2, :sc1_lift_3, :sc1_cross_1, :sc1_cross_2, :sc1_cross_3))
        SVector{3, Float64}(
            Float64(tbl.sc1_drag_1[idx]) + Float64(tbl.sc1_lift_1[idx]) + Float64(tbl.sc1_cross_1[idx]),
            Float64(tbl.sc1_drag_2[idx]) + Float64(tbl.sc1_lift_2[idx]) + Float64(tbl.sc1_cross_2[idx]),
            Float64(tbl.sc1_drag_3[idx]) + Float64(tbl.sc1_lift_3[idx]) + Float64(tbl.sc1_cross_3[idx]),
        )
    elseif :active_perturbation_force_mag in names_set
        return _saved_force_ratio(ctx, tbl, idx)
    else
        return NaN
    end
    return norm(force) / max(mass, eps(Float64)) / _central_ratio_denominator(ctx, norm(pos))
end

function _saved_aero_force_ratios(ctx::CaseContext, tbl)::Vector{Float64}
    :time in Symbol.(propertynames(tbl)) || return Float64[]
    return [_saved_aero_force_ratio(ctx, tbl, i) for i in eachindex(tbl.time)]
end

function _basic_drag(ctx::CaseContext, tbl)::Float64
    rho_p = 0.0
    if :density_kg_m3 in Symbol.(propertynames(tbl))
        finite = filter(x -> isfinite(x) && x >= 0.0, Float64.(tbl.density_kg_m3))
        !isempty(finite) && (rho_p = maximum(finite))
    end
    if rho_p <= 0.0
        finite_ratios = filter(isfinite, _saved_aero_force_ratios(ctx, tbl))
        !isempty(finite_ratios) && return maximum(finite_ratios)
    end
    return rho_p * ctx.rp_m / (2.0 * ctx.beta_kg_m2)
end

function _detailed_j2(ctx::CaseContext, pos::SVector{3, Float64}, vel::SVector{3, Float64}, t::Float64)::Float64
    frame = _planet_frame(ctx, pos, vel, t)
    r = norm(pos)
    sin_phi = sin(frame.lat_rad)
    return 1.5 * ctx.planet.J2 * (ctx.planet.Rp_e / r)^2 *
        sqrt(1.0 - 2.0 * sin_phi^2 + 5.0 * sin_phi^4)
end

function _detailed_j2_force(ctx::CaseContext, pos::SVector{3, Float64}, vel::SVector{3, Float64}, mass::Float64, t::Float64)::SVector{3, Float64}
    frame = _planet_frame(ctx, pos, vel, t)
    state = SM.StateSample(pos, vel, mass)
    env = SM.EnvironmentSample(ctx.planet; planet_frame=frame)
    j2_force, _ = SM.wrench(ctx.j2_model, state, env, t)
    central_force, _ = SM.wrench(InverseSquaredGravityModel(), state, env, t)
    return j2_force - central_force
end

function _c20_harmonics_force(ctx::CaseContext, model, pos::SVector{3, Float64}, vel::SVector{3, Float64}, mass::Float64, t::Float64)::SVector{3, Float64}
    model.L >= 2 || return SVector{3, Float64}(0.0, 0.0, 0.0)
    C20 = Float64(model.C[3, 1])
    isfinite(C20) && C20 != 0.0 || return SVector{3, Float64}(0.0, 0.0, 0.0)
    frame = _planet_frame(ctx, pos, vel, t)
    r = norm(frame.pos_pp)
    r > 0.0 || return SVector{3, Float64}(0.0, 0.0, 0.0)
    x, y, z = frame.pos_pp
    r2 = r * r
    z2 = z * z
    J2 = -sqrt(5.0) * C20
    common = 1.5 * _planet_mu(ctx.planet) * J2 * model.reference_radius_m^2 / (r2 * r2)
    accel_pp = SVector{3, Float64}(
        common * (x / r) * (5.0 * z2 / r2 - 1.0),
        common * (y / r) * (5.0 * z2 / r2 - 1.0),
        common * (z / r) * (5.0 * z2 / r2 - 3.0),
    )
    return mass * (frame.l_pi' * accel_pp)
end

function _detailed_j2_from_force(ctx::CaseContext, pos::SVector{3, Float64}, vel::SVector{3, Float64}, mass::Float64, t::Float64)::Float64
    return norm(_detailed_j2_force(ctx, pos, vel, mass, t)) /
        max(mass, eps(Float64)) / _central_ratio_denominator(ctx, norm(pos))
end

function _detailed_harmonics(ctx::CaseContext, pos::SVector{3, Float64}, vel::SVector{3, Float64}, mass::Float64, t::Float64)::Float64
    frame = _planet_frame(ctx, pos, vel, t)
    state = SM.StateSample(pos, vel, mass)
    env = SM.EnvironmentSample(ctx.planet; planet_frame=frame)
    force, _ = SM.wrench(ctx.harmonics_model, state, env, t)
    if ctx.harmonics_model.include_central
        central_force, _ = SM.wrench(InverseSquaredGravityModel(), state, env, t)
        force -= central_force
    end
    force -= _c20_harmonics_force(ctx, ctx.harmonics_model, pos, vel, mass, t)
    return norm(force) / max(mass, eps(Float64)) / _central_ratio_denominator(ctx, norm(pos))
end

function _detailed_srp(ctx::CaseContext, pos::SVector{3, Float64}, vel::SVector{3, Float64}, mass::Float64, t::Float64)::Float64
    sun_pos = _sun_position(ctx, t)
    state = SM.StateSample(pos, vel, mass)
    env = SM.EnvironmentSample(ctx.planet; solar=SM.SolarEphemerisSample(sun_pos))
    force, _ = SM.wrench(ctx.srp_model, state, env, t)
    return norm(force) / max(mass, eps(Float64)) / _central_ratio_denominator(ctx, norm(pos))
end

function _simulated_j2(ctx::CaseContext, pos::SVector{3, Float64}, vel::SVector{3, Float64}, mass::Float64, t::Float64)::Float64
    return _detailed_j2_from_force(ctx, pos, vel, mass, t)
end

function _simulated_harmonics(ctx::CaseContext, pos::SVector{3, Float64}, vel::SVector{3, Float64}, mass::Float64, t::Float64)::Float64
    return _detailed_harmonics(ctx, pos, vel, mass, t)
end

function _simulated_srp(ctx::CaseContext, pos::SVector{3, Float64}, vel::SVector{3, Float64}, mass::Float64, t::Float64)::Float64
    return _detailed_srp(ctx, pos, vel, mass, t)
end

function _detailed_third_body(ctx::CaseContext, pos::SVector{3, Float64}, t::Float64)::Float64
    r = norm(pos)
    positions = _third_body_positions(ctx, t)
    accel = SVector{3, Float64}(0.0, 0.0, 0.0)
    for (mu3, pos3) in zip(ctx.nbody_model.body_mus, positions)
        pos_sc_to_body = pos3 - pos
        accel += mu3 * (
            pos_sc_to_body / norm(pos_sc_to_body)^3 -
            pos3 / norm(pos3)^3
        )
    end
    return norm(accel) / _central_ratio_denominator(ctx, r)
end

function _instantaneous_semimajor_axis(ctx::CaseContext, pos::SVector{3, Float64}, vel::SVector{3, Float64})::Float64
    r = norm(pos)
    v2 = dot(vel, vel)
    mu = _planet_mu(ctx.planet)
    denom = 2.0 / r - v2 / mu
    abs(denom) <= eps(Float64) && return Inf
    return 1.0 / denom
end

function _effective_cda_over_m(tbl, idx::Int)::Float64
    names_set = Set(Symbol.(propertynames(tbl)))
    required = (:sc1_drag_1, :sc1_drag_2, :sc1_drag_3, :dynamic_pressure_pa, :sc1_mass)
    all(name -> name in names_set, required) || return NaN
    q = Float64(tbl.dynamic_pressure_pa[idx])
    mass = Float64(tbl.sc1_mass[idx])
    q > 0.0 && mass > 0.0 || return NaN
    drag = norm(SVector{3, Float64}(
        Float64(tbl.sc1_drag_1[idx]),
        Float64(tbl.sc1_drag_2[idx]),
        Float64(tbl.sc1_drag_3[idx]),
    ))
    drag > 0.0 || return NaN
    return drag / (q * mass)
end

function _detailed_drag(ctx::CaseContext, tbl, idx::Int, pos::SVector{3, Float64}, vel::SVector{3, Float64}, t::Float64)::Float64
    names_set = Set(Symbol.(propertynames(tbl)))
    :density_kg_m3 in names_set || return _saved_aero_force_ratio(ctx, tbl, idx)
    rho = Float64(tbl.density_kg_m3[idx])
    isfinite(rho) && rho > 0.0 || return _saved_aero_force_ratio(ctx, tbl, idx)
    cda_over_m = _effective_cda_over_m(tbl, idx)
    if !isfinite(cda_over_m) || cda_over_m <= 0.0
        cda_over_m = 1.0 / ctx.beta_kg_m2
    end
    frame = _planet_frame(ctx, pos, vel, t)
    r = norm(pos)
    a = _instantaneous_semimajor_axis(ctx, pos, vel)
    h = norm(cross(pos, vel))
    hvec = cross(pos, vel)
    cos_i = abs(h) <= eps(Float64) ? 0.0 : clamp(hvec[3] / h, -1.0, 1.0)
    mu = _planet_mu(ctx.planet)
    omega_atm = norm(_planet_omega(ctx.planet))
    cos_phi = cos(frame.lat_rad)
    bracket = mu * (2.0 / r - 1.0 / a) -
        2.0 * h * omega_atm * cos_i +
        omega_atm^2 * r^2 * cos_phi^2
    return cda_over_m * rho * r^2 / (2.0 * mu) * max(0.0, bracket)
end

function _saved_force_ratio(ctx::CaseContext, tbl, idx::Int)::Float64
    names_set = Set(Symbol.(propertynames(tbl)))
    :active_perturbation_force_mag in names_set || return NaN
    _, _, mass = _state(tbl, idx)
    r2 = Float64(tbl.sc1_pos_1[idx])^2 + Float64(tbl.sc1_pos_2[idx])^2 + Float64(tbl.sc1_pos_3[idx])^2
    denom = mass * _planet_mu(ctx.planet) / r2
    denom > 0.0 || return NaN
    ratio = Float64(tbl.active_perturbation_force_mag[idx]) / denom
    return isfinite(ratio) ? ratio : NaN
end

function _simulated_drag(ctx::CaseContext, tbl, idx::Int)::Float64
    names_set = Set(Symbol.(propertynames(tbl)))
    required = (:sc1_drag_1, :sc1_drag_2, :sc1_drag_3, :sc1_lift_1, :sc1_lift_2, :sc1_lift_3, :sc1_cross_1, :sc1_cross_2, :sc1_cross_3)
    all(name -> name in names_set, required) || return _saved_force_ratio(ctx, tbl, idx)
    pos, _, mass = _state(tbl, idx)
    force = SVector{3, Float64}(
        Float64(tbl.sc1_drag_1[idx]) + Float64(tbl.sc1_lift_1[idx]) + Float64(tbl.sc1_cross_1[idx]),
        Float64(tbl.sc1_drag_2[idx]) + Float64(tbl.sc1_lift_2[idx]) + Float64(tbl.sc1_cross_2[idx]),
        Float64(tbl.sc1_drag_3[idx]) + Float64(tbl.sc1_lift_3[idx]) + Float64(tbl.sc1_cross_3[idx]),
    )
    return norm(force) / max(mass, eps(Float64)) / _central_ratio_denominator(ctx, norm(pos))
end

function _empty_series(n::Int)
    return Dict(key => fill(NaN, n) for key in MODEL_KEYS)
end

function _analytical_series_unlocked(info, tbl, indices::AbstractVector{<:Integer})
    ctx = _context(info, tbl)
    n = length(indices)
    basic = _empty_series(n)
    detailed = _empty_series(n)

    basic_values = Dict(
        :j2 => _basic_j2(ctx),
        :harmonics => _basic_harmonics(ctx),
        :srp => _basic_srp(ctx),
        :third_body => _basic_third_body(ctx),
        :drag => _basic_drag(ctx, tbl),
    )
    basic_values[:full_environment] =
        basic_values[:j2] + basic_values[:harmonics] + basic_values[:srp] +
        basic_values[:third_body] + basic_values[:drag]

    for key in keys(basic_values)
        basic[key] .= basic_values[key]
    end

    for (out_idx, idx_raw) in enumerate(indices)
        idx = Int(idx_raw)
        t = Float64(tbl.time[idx])
        pos, vel, mass = _state(tbl, idx)
        detailed[:j2][out_idx] = _detailed_j2(ctx, pos, vel, t)
        detailed[:harmonics][out_idx] = _detailed_harmonics(ctx, pos, vel, mass, t)
        detailed[:srp][out_idx] = _detailed_srp(ctx, pos, vel, mass, t)
        detailed[:third_body][out_idx] = _detailed_third_body(ctx, pos, t)
        detailed[:drag][out_idx] = _detailed_drag(ctx, tbl, idx, pos, vel, t)
        detailed[:full_environment][out_idx] = sum(v for v in (
            detailed[:j2][out_idx],
            detailed[:harmonics][out_idx],
            detailed[:srp][out_idx],
            detailed[:third_body][out_idx],
            detailed[:drag][out_idx],
        ) if isfinite(v))
    end

    return (; basic, detailed)
end

function analytical_series(info, tbl, indices::AbstractVector{<:Integer})
    lock(ANALYTICAL_SPICE_LOCK) do
        _analytical_series_unlocked(info, tbl, indices)
    end
end

function _simulated_model_ratio_series_unlocked(info, tbl, model_key::Symbol, indices::AbstractVector{<:Integer})
    ctx = _context(info, tbl)
    out = fill(NaN, length(indices))
    for (out_idx, idx_raw) in enumerate(indices)
        idx = Int(idx_raw)
        t = Float64(tbl.time[idx])
        pos, vel, mass = _state(tbl, idx)
        out[out_idx] = if model_key == :j2
            _simulated_j2(ctx, pos, vel, mass, t)
        elseif model_key == :harmonics
            _simulated_harmonics(ctx, pos, vel, mass, t)
        elseif model_key == :srp
            _simulated_srp(ctx, pos, vel, mass, t)
        elseif model_key == :third_body
            _saved_force_ratio(ctx, tbl, idx)
        elseif model_key == :drag
            _simulated_drag(ctx, tbl, idx)
        else
            _saved_force_ratio(ctx, tbl, idx)
        end
    end
    return out
end

function simulated_model_ratio_series(info, tbl, model_key::Symbol, indices::AbstractVector{<:Integer})
    lock(ANALYTICAL_SPICE_LOCK) do
        _simulated_model_ratio_series_unlocked(info, tbl, model_key, indices)
    end
end

end
