using StaticArrays
using LinearAlgebra

@inline function _extract_sample_pos_vel(x)
    if hasproperty(x, :pos_ii) && hasproperty(x, :vel_ii)
        return x.pos_ii, x.vel_ii
    end
    if hasproperty(x, :pos) && hasproperty(x, :vel)
        return SVector{3, Float64}(x[1], x[2], x[3]), SVector{3, Float64}(x[4], x[5], x[6])
    end
    return SVector{3, Float64}(x[1], x[2], x[3]), SVector{3, Float64}(x[4], x[5], x[6])
end

@inline function _extract_sample_mass_kg(x)::Float64
    if hasproperty(x, :mass_kg)
        return Float64(x.mass_kg)
    end
    if hasproperty(x, :mass)
        return Float64(x[7])
    end
    return length(x) >= 7 ? Float64(x[7]) : NaN
end

@inline function build_state_sample(sc_view, spacecraft, orientation_sim::Bool)::StateSample
    pos_ii, vel_ii = _extract_sample_pos_vel(sc_view)
    q_ib = (orientation_sim && hasproperty(sc_view, :q)) ? SVector{4, Float64}(sc_view.q) : nothing
    ω_body = (orientation_sim && hasproperty(sc_view, :ω)) ? SVector{3, Float64}(sc_view.ω) : nothing
    return StateSample(
        pos_ii,
        vel_ii,
        _extract_sample_mass_kg(sc_view);
        q_ib=q_ib,
        ω_body=ω_body,
        spacecraft=spacecraft,
    )
end

@inline function _planet_lpi_at_engine(p, t::Float64)::SMatrix{3, 3, Float64, 9}
    return SimulationModel.SimulationCallbacks._planet_lpi_at(p, t)
end

@inline function sample_planet_frame(x, p, sat_idx::Int, t::Float64)::PlanetFrameSample
    pos_ii, vel_ii = _extract_sample_pos_vel(x)
    planet = p.args.environment_model.planet
    l_pi = _planet_lpi_at_engine(p, t)
    pos_pp, vel_pp = SimulationModel.SimulationCallbacks._planet_relative_state(pos_ii, vel_ii, planet, l_pi)
    alt, lat, lon = SimulationModel.SimulationCallbacks.rtolatlong(pos_pp, planet)
    return PlanetFrameSample(l_pi, pos_pp, vel_pp, alt, lat, lon)
end

# Variant that accepts a pre-computed l_pi so the caller avoids acquiring harmonics_lpi_lock.
# Used in the parallel atmosphere pre-sample phase where l_pi is computed once before @batch.
@inline function sample_planet_frame_with_lpi(x, planet, l_pi::SMatrix{3, 3, Float64, 9})::PlanetFrameSample
    pos_ii, vel_ii = _extract_sample_pos_vel(x)
    pos_pp, vel_pp = SimulationModel.SimulationCallbacks._planet_relative_state(pos_ii, vel_ii, planet, l_pi)
    alt, lat, lon = SimulationModel.SimulationCallbacks.rtolatlong(pos_pp, planet)
    return PlanetFrameSample(l_pi, pos_pp, vel_pp, alt, lat, lon)
end

@inline function _sample_atmosphere_from_planet_frame(
    x,
    planet_frame::PlanetFrameSample,
    p,
    sat_idx::Int,
    t::Float64;
    write_buffers::Bool=true,
)::AtmosphereSample
    callbacks = SimulationModel.SimulationCallbacks
    static_cfg = p.shared_buffers.density_static_config[]
    cache_cfg, stats_enabled, target_include_j2 = if static_cfg === nothing
        # Fallback for callers that construct ODEParams without the standard
        # execution.jl init sequence (e.g. some test harnesses): correct but
        # re-parses ENV every call, same as before this cache was added.
        (
            callbacks._gram_track_cache_config(),
            callbacks._gram_runtime_stats_enabled(),
            callbacks._gram_track_cache_target_use_j2() &&
                callbacks._uses_j2_gravity_effector(p.args.dynamics_model.dynamic_effectors),
        )
    else
        (static_cfg.cache_cfg, static_cfg.stats_enabled, static_cfg.target_include_j2)
    end
    density_model = callbacks._density_model_for_sat(p, sat_idx)
    caches = p.shared_buffers.gram_density_cache
    pos_ii, vel_ii = _extract_sample_pos_vel(x)
    current_mass_kg = _extract_sample_mass_kg(x)
    rho, T, wind_vec = callbacks._density_state_from_kinematics!(
        p,
        sat_idx,
        pos_ii,
        vel_ii,
        current_mass_kg,
        planet_frame.alt_m,
        planet_frame.lat_rad,
        planet_frame.lon_rad,
        t,
        density_model,
        cache_cfg,
        stats_enabled,
        target_include_j2,
        caches,
    )
    if write_buffers
        callbacks._write_density_buffers!(p, sat_idx, rho, T, wind_vec, t)
    end
    return AtmosphereSample(rho, T, wind_vec)
end

@inline function sample_atmosphere(x, p, sat_idx::Int, t::Float64; write_buffers::Bool=true)::AtmosphereSample
    planet_frame = sample_planet_frame(x, p, sat_idx, t)
    return _sample_atmosphere_from_planet_frame(
        x,
        planet_frame,
        p,
        sat_idx,
        t;
        write_buffers=write_buffers,
    )
end

@inline function _buffered_atmosphere_valid(p, sat_idx::Int, t::Float64)::Bool
    times = p.shared_buffers.density_sample_t
    return sat_idx <= length(times) && times[sat_idx] == t
end

@inline function sample_buffered_atmosphere(x, p, sat_idx::Int, t::Float64)::AtmosphereSample
    if !_buffered_atmosphere_valid(p, sat_idx, t)
        return sample_atmosphere(x, p, sat_idx, t; write_buffers=true)
    end
    rho = sat_idx <= length(p.shared_buffers.densities) ? p.shared_buffers.densities[sat_idx] : 0.0
    T = sat_idx <= length(p.shared_buffers.temperatures) ? p.shared_buffers.temperatures[sat_idx] : p.args.environment_model.planet.T_ref
    wind_vec = sat_idx <= length(p.shared_buffers.winds) ? p.shared_buffers.winds[sat_idx] : SVector{3, Float64}(0.0, 0.0, 0.0)
    return AtmosphereSample(rho, T, wind_vec)
end

@inline function sample_buffered_planet_frame(p, sat_idx::Int)::PlanetFrameSample
    return PlanetFrameSample(
        p.shared_buffers.rhs_flat_planet_lpi[],
        p.shared_buffers.rhs_flat_planet_pos_pp[][sat_idx],
        p.shared_buffers.rhs_flat_planet_vel_pp[][sat_idx],
        p.shared_buffers.rhs_flat_planet_alt_m[][sat_idx],
        p.shared_buffers.rhs_flat_planet_lat_rad[][sat_idx],
        p.shared_buffers.rhs_flat_planet_lon_rad[][sat_idx],
    )
end

@inline function sample_solar_ephemeris(x, p, sat_idx::Int, t::Float64)::SolarEphemerisSample
    planet = p.args.environment_model.planet
    et = p.shared_buffers.et_start[] + t
    primary_body_name = SimulationModel.DynamicEffectors._spice_query_name(planet.name)
    spice_rhs_memo_enabled = p.shared_buffers.spice_rhs_memo_enabled[]
    spice_rhs_memo = p.shared_buffers.spice_rhs_memo
    perturbation_effectors = SimulationModel.DynamicEffectors.PerturbationEffectors
    cache_entry = p.shared_buffers.srp_sun_ephemeris_cache[]
    pos_primary_sun_j2000_m = if cache_entry isa SimulationModel.SRPSunEphemerisCache
        cached = SimulationModel.DynamicEffectors._srp_sun_position_from_cache_j2000_m(cache_entry, et)
        cached === nothing ?
            perturbation_effectors._srp_sun_position_from_spice_j2000_m(
                et,
                primary_body_name,
                spice_rhs_memo_enabled,
                spice_rhs_memo,
                p.shared_buffers.spice_runtime_counters.srp_spkpos_runtime_calls,
            ) :
            cached
    else
        perturbation_effectors._srp_sun_position_from_spice_j2000_m(
            et,
            primary_body_name,
            spice_rhs_memo_enabled,
            spice_rhs_memo,
            p.shared_buffers.spice_runtime_counters.srp_spkpos_runtime_calls,
        )
    end
    return SolarEphemerisSample(SVector{3, Float64}(pos_primary_sun_j2000_m))
end

@inline function sample_third_body_ephemerides(model::SimulationModel.NBodyGravityModel, x, p, sat_idx::Int, t::Float64)
    et = p.shared_buffers.et_start[] + t
    primary_body_name = SimulationModel.DynamicEffectors._spice_query_name(model.primary_body_name)
    spice_rhs_memo_enabled = p.shared_buffers.spice_rhs_memo_enabled[]
    spice_rhs_memo = p.shared_buffers.spice_rhs_memo
    cache_entry = p.shared_buffers.nbody_ephemeris_cache[]
    perturbation_effectors = SimulationModel.DynamicEffectors.PerturbationEffectors
    positions_ii = ntuple(length(model.body_names)) do k
        body_name_spice = SimulationModel.DynamicEffectors._spice_query_name(model.body_names[k])
        pos_primary_body_j2000_m = if cache_entry isa SimulationModel.NBodyEphemerisCache
            cached = SimulationModel.DynamicEffectors._nbody_body_position_from_cache_j2000_m(
                cache_entry,
                et,
                body_name_spice,
                primary_body_name,
            )
            cached === nothing ?
                perturbation_effectors._nbody_body_position_from_spice_j2000_m(
                    body_name_spice,
                    et,
                    primary_body_name,
                    spice_rhs_memo_enabled,
                    spice_rhs_memo,
                    p.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls,
                ) :
                cached
        else
            perturbation_effectors._nbody_body_position_from_spice_j2000_m(
                body_name_spice,
                et,
                primary_body_name,
                spice_rhs_memo_enabled,
                spice_rhs_memo,
                p.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls,
            )
        end
        return SVector{3, Float64}(pos_primary_body_j2000_m)
    end
    return ThirdBodyEphemerisSample(model.body_names, positions_ii)
end

@inline function sample_environment(
    req::EffectorEnvironmentRequirements,
    model,
    x,
    p,
    sat_idx::Int,
    t::Float64;
    write_buffers::Bool=false,
)::EnvironmentSample
    planet = p.args.environment_model.planet
    need_planet_frame = req.planet_frame || req.atmosphere
    planet_frame = need_planet_frame ? sample_planet_frame(x, p, sat_idx, t) : nothing
    atmosphere = req.atmosphere ? _sample_atmosphere_from_planet_frame(x, planet_frame, p, sat_idx, t; write_buffers=write_buffers) : nothing
    solar = req.solar ? sample_solar_ephemeris(x, p, sat_idx, t) : nothing
    third_bodies = isempty(req.third_body_names) ? nothing : sample_third_body_ephemerides(model, x, p, sat_idx, t)
    return EnvironmentSample(
        planet;
        planet_frame=req.planet_frame ? planet_frame : nothing,
        atmosphere=atmosphere,
        solar=solar,
        third_bodies=third_bodies,
    )
end

@inline function _sample_reusable_planet_frame(req::EffectorEnvironmentRequirements, x, p, sat_idx::Int, t::Float64)
    if req.planet_frame || req.atmosphere
        return p.shared_buffers.rhs_planet_frame_prefilled[] ?
            sample_buffered_planet_frame(p, sat_idx) :
            sample_planet_frame(x, p, sat_idx, t)
    end
    return nothing
end

@inline function _sample_reusable_atmosphere(req::EffectorEnvironmentRequirements, x, planet_frame, p, sat_idx::Int, t::Float64)
    req.atmosphere || return nothing
    if p.shared_buffers.rhs_atmosphere_prefilled[]
        return sample_buffered_atmosphere(x, p, sat_idx, t)
    end
    return _sample_atmosphere_from_planet_frame(x, planet_frame, p, sat_idx, t; write_buffers=false)
end

@inline function _sample_reusable_solar(req::EffectorEnvironmentRequirements, x, p, sat_idx::Int, t::Float64)
    req.solar || return nothing
    return (p.shared_buffers.rhs_solar_prefilled[] && p.shared_buffers.rhs_flat_solar_t[] == t) ?
        SolarEphemerisSample(p.shared_buffers.rhs_flat_solar_pos_ii[]) :
        sample_solar_ephemeris(x, p, sat_idx, t)
end

# Variant of sample_environment that reads reusable flat-RHS component buffers when
# available. Used by wrench-based effectors so repeated effectors for the same satellite
# do not rebuild the same planet-frame, atmosphere, or solar samples.
@inline function sample_environment_with_reusable_buffers(
    req::EffectorEnvironmentRequirements,
    model,
    x,
    p,
    sat_idx::Int,
    t::Float64,
)::EnvironmentSample
    planet = p.args.environment_model.planet
    sampled_planet_frame = _sample_reusable_planet_frame(req, x, p, sat_idx, t)
    atmosphere = _sample_reusable_atmosphere(req, x, sampled_planet_frame, p, sat_idx, t)
    solar = _sample_reusable_solar(req, x, p, sat_idx, t)
    third_bodies = isempty(req.third_body_names) ? nothing : sample_third_body_ephemerides(model, x, p, sat_idx, t)
    return EnvironmentSample(
        planet;
        planet_frame=req.planet_frame ? sampled_planet_frame : nothing,
        atmosphere=atmosphere,
        solar=solar,
        third_bodies=third_bodies,
    )
end

@inline function sample_environment_with_buffered_atm(
    req::EffectorEnvironmentRequirements,
    model,
    x,
    p,
    sat_idx::Int,
    t::Float64,
)::EnvironmentSample
    return sample_environment_with_reusable_buffers(req, model, x, p, sat_idx, t)
end

@inline function _wrench_method_available(effector::SimulationModel.DynamicEffectors.GravitationalHarmonicsModel)::Bool
    # The legacy RHS path reuses per-satellite harmonics scratch buffers; the
    # generic wrench hook allocates a scratch workspace per call.
    return false
end

@inline function _wrench_method_available(effector)::Bool
    return hasmethod(
        SimulationModel.wrench,
        Tuple{typeof(effector), SimulationModel.StateSample, SimulationModel.EnvironmentSample, Float64},
    )
end
