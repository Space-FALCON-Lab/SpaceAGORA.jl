@inline function srp_cannonball_accel(
    pos_ii::SVector{3, Float64},
    pos_primary_sun::SVector{3, Float64},
    planet_radius_m::Float64,
    p_srp_unscaled::Float64,
    reflection_coefficient::Float64,
    reference_area_m2::Float64,
    mass_kg::Float64;
    AU_m::Float64=149_597_870_700.0
)::SVector{3, Float64}
    if !(isfinite(mass_kg) && mass_kg > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(reference_area_m2) && reference_area_m2 > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(reflection_coefficient) && reflection_coefficient >= 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    r_sun_to_sc = pos_ii - pos_primary_sun
    r_sun_to_sc_mag = norm(r_sun_to_sc)
    if !(isfinite(r_sun_to_sc_mag) && r_sun_to_sc_mag > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    eclipse_ratio = eclipse_area_calc(pos_ii, pos_primary_sun, planet_radius_m)
    P_srp = p_srp_unscaled * (AU_m / r_sun_to_sc_mag)^2
    return _cannonball_radiation_accel(
        r_sun_to_sc,
        P_srp * eclipse_ratio,
        reflection_coefficient,
        reference_area_m2,
        mass_kg,
    )
end

@inline function _cannonball_radiation_accel(
    direction_ii::SVector{3, Float64},
    pressure_n_m2::Float64,
    reflection_coefficient::Float64,
    reference_area_m2::Float64,
    mass_kg::Float64,
)::SVector{3, Float64}
    if !(isfinite(mass_kg) && mass_kg > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(reference_area_m2) && reference_area_m2 > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(reflection_coefficient) && reflection_coefficient >= 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(pressure_n_m2) && pressure_n_m2 > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    direction_mag = norm(direction_ii)
    if !(isfinite(direction_mag) && direction_mag > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    scale = reflection_coefficient * reference_area_m2 * pressure_n_m2 / mass_kg / direction_mag
    return scale * direction_ii
end

@inline function _lambert_phase_function(alpha_rad::Float64)::Float64
    if !isfinite(alpha_rad)
        return 0.0
    end
    α = clamp(alpha_rad, 0.0, π)
    return max(0.0, (sin(α) + (π - α) * cos(α)) / π)
end

@inline function planetary_albedo_accel(
    pos_ii::SVector{3, Float64},
    pos_primary_sun::SVector{3, Float64},
    planet_radius_m::Float64,
    p_srp_unscaled::Float64,
    reflection_coefficient::Float64,
    reference_area_m2::Float64,
    mass_kg::Float64,
    planet_albedo::Float64;
    AU_m::Float64=149_597_870_700.0
)::SVector{3, Float64}
    if !(isfinite(planet_radius_m) && planet_radius_m > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(planet_albedo) && planet_albedo > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    spacecraft_range_m = norm(pos_ii)
    sun_range_m = norm(pos_primary_sun)
    if !(isfinite(spacecraft_range_m) && spacecraft_range_m > planet_radius_m)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(sun_range_m) && sun_range_m > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    cos_alpha = clamp(dot(pos_primary_sun, pos_ii) / (sun_range_m * spacecraft_range_m), -1.0, 1.0)
    phase = _lambert_phase_function(acos(cos_alpha))
    if phase <= 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    P_solar_at_planet = p_srp_unscaled * (AU_m / sun_range_m)^2
    pressure = P_solar_at_planet * planet_albedo * phase * (planet_radius_m / spacecraft_range_m)^2
    return _cannonball_radiation_accel(
        pos_ii,
        pressure,
        reflection_coefficient,
        reference_area_m2,
        mass_kg,
    )
end

@inline function planetary_ir_accel(
    pos_ii::SVector{3, Float64},
    planet_radius_m::Float64,
    reflection_coefficient::Float64,
    reference_area_m2::Float64,
    mass_kg::Float64,
    planet_ir_flux_w_m2::Float64,
)::SVector{3, Float64}
    if !(isfinite(planet_radius_m) && planet_radius_m > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(planet_ir_flux_w_m2) && planet_ir_flux_w_m2 > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    spacecraft_range_m = norm(pos_ii)
    if !(isfinite(spacecraft_range_m) && spacecraft_range_m > planet_radius_m)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    pressure = (planet_ir_flux_w_m2 / 299_792_458.0) * (planet_radius_m / spacecraft_range_m)^2
    return _cannonball_radiation_accel(
        pos_ii,
        pressure,
        reflection_coefficient,
        reference_area_m2,
        mass_kg,
    )
end

function srp(
    planet,
    p_srp_unscaled::Float64,
    reflection_coefficient::Float64,
    reference_area_m2::Float64,
    mass_kg::Float64,
    pos_ii::AbstractVector,
    et::Float64;
    AU_m::Float64=149_597_870_700.0
)
    pos_ii_sv = SVector{3, Float64}(pos_ii)
    primary_body_name = _spice_query_name(planet.name)
    pos_primary_sun_j2000_m = spice_position_j2000_m("sun", et, primary_body_name)
    pos_primary_sun = SVector{3, Float64}(pos_primary_sun_j2000_m)
    return srp_cannonball_accel(
        pos_ii_sv,
        pos_primary_sun,
        planet.Rp_e,
        p_srp_unscaled,
        reflection_coefficient,
        reference_area_m2,
        mass_kg;
        AU_m=AU_m
    )
end

function calcForceTorque(model::SolarRadiationPressureModel, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    if model.A == 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    planet = param.args.environment_model.planet
    pos_ii = SVector{3, Float64}(x[1], x[2], x[3])
    mass = Float64(x[7])
    et = param.shared_buffers.et_start[] + param.shared_buffers.current_time[]
    primary_body_name = _spice_query_name(planet.name)
    spice_rhs_memo_enabled = param.shared_buffers.spice_rhs_memo_enabled[]
    spice_rhs_memo = param.shared_buffers.spice_rhs_memo
    pos_primary_sun = if model.direct || model.albedo
        cache_entry = param.shared_buffers.srp_sun_ephemeris_cache[]
        pos_primary_sun_j2000_m = if cache_entry isa SRPSunEphemerisCache
            cached = _srp_sun_position_from_cache_j2000_m(cache_entry, et)
            cached === nothing ? _srp_sun_position_from_spice_j2000_m(et, primary_body_name, spice_rhs_memo_enabled, spice_rhs_memo, param.shared_buffers.spice_runtime_counters.srp_spkpos_runtime_calls) : cached
        else
            _srp_sun_position_from_spice_j2000_m(et, primary_body_name, spice_rhs_memo_enabled, spice_rhs_memo, param.shared_buffers.spice_runtime_counters.srp_spkpos_runtime_calls)
        end
        SVector{3, Float64}(pos_primary_sun_j2000_m)
    else
        SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    accel_ii = _srp_total_acceleration_ii(model, planet, pos_ii, pos_primary_sun, mass)
    force_ii = mass * SVector{3, Float64}(accel_ii)
    return force_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline environment_requirements(model::SolarRadiationPressureModel) = EffectorEnvironmentRequirements(solar=(model.direct || model.albedo))
@inline gravity_backbone_kick_structure(::SolarRadiationPressureModel) = :velocity_kick_explicit

@noinline function _srp_force_ii(model::SolarRadiationPressureModel, planet, pos_ii::SVector{3, Float64}, sun_pos_ii::SVector{3, Float64}, mass::Float64)::SVector{3, Float64}
    return mass * _srp_total_acceleration_ii(model, planet, pos_ii, sun_pos_ii, mass)
end

@inline function _srp_total_acceleration_ii(
    model::SolarRadiationPressureModel,
    planet,
    pos_ii::SVector{3, Float64},
    pos_primary_sun::SVector{3, Float64},
    mass_kg::Float64,
)::SVector{3, Float64}
    if model.A == 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    accel_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    if model.direct
        accel_ii .+= srp_cannonball_accel(
            pos_ii,
            pos_primary_sun,
            planet.Rp_e,
            4.56e-6,
            model.Cr,
            model.A,
            mass_kg;
            AU_m=model.AU_m
        )
    end
    if model.albedo
        accel_ii .+= planetary_albedo_accel(
            pos_ii,
            pos_primary_sun,
            planet.Rp_e,
            4.56e-6,
            model.Cr,
            model.A,
            mass_kg,
            model.planet_albedo;
            AU_m=model.AU_m
        )
    end
    if model.ir
        accel_ii .+= planetary_ir_accel(
            pos_ii,
            planet.Rp_e,
            model.Cr,
            model.A,
            mass_kg,
            model.planet_ir_flux_w_m2,
        )
    end
    return SVector{3, Float64}(accel_ii)
end

@inline function wrench(
    model::SolarRadiationPressureModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    if model.A == 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    pos_primary_sun = if model.direct || model.albedo
        solar = env.solar
        solar === nothing && throw(ArgumentError("SolarRadiationPressureModel wrench requires env.solar."))
        solar.sun_pos_ii
    else
        SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    return _srp_force_ii(model, env.planet, x.pos_ii, pos_primary_sun, x.mass_kg), SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline function gravity_backbone_kick_acceleration_ii(
    model::SolarRadiationPressureModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::SVector{3, Float64}
    pos_primary_sun = if model.direct || model.albedo
        solar = env.solar
        solar === nothing && throw(ArgumentError("SolarRadiationPressureModel gravity_backbone kick requires env.solar."))
        solar.sun_pos_ii
    else
        SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    return _srp_total_acceleration_ii(model, env.planet, x.pos_ii, pos_primary_sun, x.mass_kg)
end

"""
    FacetSolarRadiationPressureModel(; AU_m=149_597_870_700.0, p_srp_1au=4.56e-6)

Per-facet solar radiation pressure on the flat plates attached to each link with
[`add_facet!`](@ref) (`link.SRP_facets`). For an illuminated facet of area `A`,
unit normal `n̂`, specular coefficient `δ` and diffuse coefficient `ρ`,

    F = -P A cosθ [(1 - δ) ŝ + 2 (ρ/3 + δ cosθ) n̂] · ν,    cosθ = n̂ · ŝ > 0,

where `ŝ` points from the spacecraft to the Sun, `P` is `p_srp_1au` scaled to the
spacecraft's distance from the Sun, and `ν` is the shadow fraction from
`eclipse_area_calc` (Montenbruck and Gill, *Satellite Orbits*, Sec. 3.4). This is
the facet model of the former `srp!`, now evaluated through the `wrench` hook.

Each facet's force acts at its centre of pressure (`facet.cp`, link frame), so the
model also returns the torque about the root body's centre of mass in the root
body frame. Facets are oriented by the propagated attitude; without one
(`orientation_sim = false`) the body axes are taken as J2000 and no torque is
returned. A spacecraft with no facets feels no force from this model.
"""
struct FacetSolarRadiationPressureModel <: AbstractForceTorqueModel
    AU_m::Float64
    p_srp_1au::Float64
    function FacetSolarRadiationPressureModel(AU_m::Real, p_srp_1au::Real)
        AU_f, p_f = Float64(AU_m), Float64(p_srp_1au)
        (isfinite(AU_f) && AU_f > 0.0) || throw(ArgumentError("FacetSolarRadiationPressureModel.AU_m must be > 0 m, got $AU_f."))
        (isfinite(p_f) && p_f >= 0.0) || throw(ArgumentError("FacetSolarRadiationPressureModel.p_srp_1au must be >= 0 N/m^2, got $p_f."))
        return new(AU_f, p_f)
    end
end
FacetSolarRadiationPressureModel(; AU_m::Real=149_597_870_700.0, p_srp_1au::Real=4.56e-6) =
    FacetSolarRadiationPressureModel(AU_m, p_srp_1au)

@inline environment_requirements(::FacetSolarRadiationPressureModel) = EffectorEnvironmentRequirements(solar=true)

@inline function wrench(
    model::FacetSolarRadiationPressureModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    zero3 = SVector{3, Float64}(0.0, 0.0, 0.0)
    spacecraft = x.spacecraft
    spacecraft === nothing && return zero3, zero3
    solar = env.solar
    solar === nothing && throw(ArgumentError("FacetSolarRadiationPressureModel wrench requires env.solar."))
    sun_pos_ii = solar.sun_pos_ii
    r_sc_sun = sun_pos_ii - x.pos_ii
    d_sun = norm(r_sc_sun)
    (isfinite(d_sun) && d_sun > 0.0) || return zero3, zero3
    shadow = eclipse_area_calc(x.pos_ii, sun_pos_ii, env.planet.Rp_e)
    shadow == 0.0 && return zero3, zero3
    s_hat = r_sc_sun / d_sun
    P = model.p_srp_1au * (model.AU_m / d_sun)^2 * shadow
    has_attitude = x.q_ib !== nothing
    R_ib = has_attitude ? rot(x.q_ib) : SMatrix{3, 3, Float64}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

    force_ii = zero3
    torque_body = zero3
    @inbounds for link in spacecraft.links
        isempty(link.SRP_facets) && continue
        # Link frame -> root body frame (the aero wrench's convention).
        R_root_link = link.root ? SMatrix{3, 3, Float64}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0) :
            rot(SVector{4, Float64}(link.q...))'
        for facet in link.SRP_facets
            n_link = rot(SVector{4, Float64}(facet.attitude))' * SVector{3, Float64}(facet.normal_vector)
            n_ii = normalize(R_ib' * (R_root_link * n_link))
            cosθ = dot(n_ii, s_hat)
            cosθ > 0.0 || continue
            F_ii = -P * facet.area * cosθ * ((1.0 - facet.δ) * s_hat + 2.0 * (facet.ρ / 3.0 + facet.δ * cosθ) * n_ii)
            force_ii += F_ii
            if has_attitude
                cp = SVector{3, Float64}(facet.cp)
                lever_root = link.root ? cp : SVector{3, Float64}(link.r) + R_root_link * cp
                torque_body += cross(lever_root, R_ib * F_ii)
            end
        end
    end
    return force_ii, torque_body
end

# Per-link kernel at a live link pose (articulated spacecraft, opt-in; see `link_wrench`): the facets of ONE
# link, with the link's own COM position (eclipse and distance to the Sun) and attitude. Force at the link
# COM and torque about it, both inertial. Link `r`/`q` are never read.
@inline link_wrench_capable(::FacetSolarRadiationPressureModel) = true

function link_wrench(
    model::FacetSolarRadiationPressureModel,
    link,
    xl::LinkStateSample,
    env::EnvironmentSample,
    t::Float64,
    p,
    sat_idx::Int,
)::NTuple{5, SVector{3, Float64}}
    zero3 = SVector{3, Float64}(0.0, 0.0, 0.0)
    zero5 = (zero3, zero3, zero3, zero3, zero3)
    isempty(link.SRP_facets) && return zero5
    solar = env.solar
    solar === nothing && throw(ArgumentError("FacetSolarRadiationPressureModel link wrench requires env.solar."))
    sun_pos_ii = solar.sun_pos_ii
    r_sc_sun = sun_pos_ii - xl.pos_ii
    d_sun = norm(r_sc_sun)
    (isfinite(d_sun) && d_sun > 0.0) || return zero5
    shadow = eclipse_area_calc(xl.pos_ii, sun_pos_ii, env.planet.Rp_e)
    shadow == 0.0 && return zero5
    s_hat = r_sc_sun / d_sun
    P = model.p_srp_1au * (model.AU_m / d_sun)^2 * shadow
    R_link = rot(xl.q_ib)'                                 # link frame -> inertial
    force_ii = zero3
    torque_ii = zero3
    @inbounds for facet in link.SRP_facets
        n_link = rot(SVector{4, Float64}(facet.attitude))' * SVector{3, Float64}(facet.normal_vector)
        n_ii = normalize(R_link * n_link)
        cosθ = dot(n_ii, s_hat)
        cosθ > 0.0 || continue
        F_ii = -P * facet.area * cosθ * ((1.0 - facet.δ) * s_hat + 2.0 * (facet.ρ / 3.0 + facet.δ * cosθ) * n_ii)
        force_ii += F_ii
        torque_ii += cross(R_link * SVector{3, Float64}(facet.cp), F_ii)
    end
    return force_ii, torque_ii, zero3, zero3, zero3
end

@inline function _srp_sun_position_from_spice_j2000_m(
    et::Float64,
    primary_body_name::String,
    memo_enabled::Bool,
    memo::SpiceRhsMemo,
    counter::Base.Threads.Atomic{Int64}
)::SVector{3, Float64}
    return memo_enabled ?
        _srp_sun_position_from_spice_memoized_j2000_m(et, primary_body_name, memo, counter) :
        _srp_sun_position_from_spice_direct_j2000_m(et, primary_body_name, counter)
end

@inline function _srp_sun_position_from_spice_direct_j2000_m(
    et::Float64,
    primary_body_name::String,
    counter::Base.Threads.Atomic{Int64}
)::SVector{3, Float64}
    Base.Threads.atomic_add!(counter, 1)
    return spice_position_j2000_m("sun", et, primary_body_name)
end

@inline function _srp_sun_position_from_spice_memoized_j2000_m(
    et::Float64,
    primary_body_name::String,
    memo::SpiceRhsMemo,
    counter::Base.Threads.Atomic{Int64}
)::SVector{3, Float64}
    sun_query_name = "sun"
    return lock(memo.lock) do
        _spice_rhs_memo_reset_if_stale!(memo, et, primary_body_name)
        if haskey(memo.body_positions_j2000_m, sun_query_name)
            return memo.body_positions_j2000_m[sun_query_name]
        end
        Base.Threads.atomic_add!(counter, 1)
        position = spice_position_j2000_m(sun_query_name, et, primary_body_name)
        memo.body_positions_j2000_m[sun_query_name] = position
        return position
    end
end

@noinline function _srp_sun_position_from_cache_j2000_m(cache::SRPSunEphemerisCache, et::Float64)::Union{Nothing, SVector{3, Float64}}
    ets = cache.ets
    n_samples = length(ets)
    n_samples >= 2 || return nothing
    if et < ets[1] || et > ets[n_samples]
        return nothing
    end

    idx = searchsortedlast(ets, et)
    if idx <= 0
        return nothing
    elseif idx >= n_samples
        return cache.positions_j2000_m[n_samples]
    end

    et0 = ets[idx]
    et1 = ets[idx + 1]
    if et1 <= et0
        return cache.positions_j2000_m[idx]
    end
    tau = (et - et0) / (et1 - et0)
    p1 = cache.positions_j2000_m[idx]
    p2 = cache.positions_j2000_m[idx + 1]
    if idx <= 1 || (idx + 2) > n_samples
        return _interp_vec3_linear(p1, p2, tau)
    end
    p0 = cache.positions_j2000_m[idx - 1]
    p3 = cache.positions_j2000_m[idx + 2]
    return _interp_vec3_catmull_rom(p0, p1, p2, p3, tau)
end
