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

# The same answer as SimulationCallbacks._uses_j2_gravity_effector, by peeling
# the effector tuple instead of iterating it: iterating a heterogeneous tuple
# yields an abstract element and boxed it once per spacecraft per sample.
@inline _effectors_use_j2_gravity(::Tuple{})::Bool = false
@inline _effectors_use_j2_gravity(effs::Tuple)::Bool =
    first(effs) isa SimulationModel.InverseSquaredJ2GravityModel || _effectors_use_j2_gravity(Base.tail(effs))

@inline function _sample_atmosphere_from_planet_frame(
    x,
    planet_frame::PlanetFrameSample,
    p,
    sat_idx::Int,
    t::Float64;
    write_buffers::Bool=true,
)::AtmosphereSample
    callbacks = SimulationModel.SimulationCallbacks
    cb_env = callbacks._callback_env_config(p)
    density_model = callbacks._density_model_for_sat(p, sat_idx)
    # Freeze-per-step mode (see CallbackEnvConfig.density_freeze_per_step): reuse
    # the once-per-accepted-step sample from shared_buffers for every caller and
    # call site, not just the sample_buffered_atmosphere path -- wrench-based
    # effectors (the `AerodynamicCoefficientfM` path taken outside flat-mode
    # atmosphere prefill) call straight into this function with
    # write_buffers=false, bypassing that path entirely.
    # A fixed grid still varies with position and must check its spatial domain.
    if cb_env.density_freeze_per_step && !(density_model isa SimulationModel.EnvironmentModels.GRAMGridAtmosphereModel)
        times = p.shared_buffers.density_sample_t
        if sat_idx <= length(times) && isfinite(times[sat_idx])
            rho = sat_idx <= length(p.shared_buffers.densities) ? p.shared_buffers.densities[sat_idx] : 0.0
            T = sat_idx <= length(p.shared_buffers.temperatures) ? p.shared_buffers.temperatures[sat_idx] : p.args.environment_model.planet.T_ref
            wind_vec = sat_idx <= length(p.shared_buffers.winds) ? p.shared_buffers.winds[sat_idx] : SVector{3, Float64}(0.0, 0.0, 0.0)
            return AtmosphereSample(rho, T, wind_vec)
        end
    end
    cache_cfg = cb_env.gram_track_cache
    stats_enabled = cb_env.gram_runtime_stats_enabled
    target_include_j2 = cb_env.gram_track_cache_target_use_j2 &&
        _effectors_use_j2_gravity(p.args.dynamics_model.dynamic_effectors)
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


# ── Uniform-light atmosphere batch (flat-route pre-sample) ───────────────────
#
# The flat route's atmosphere pre-sample used to resolve, per spacecraft, the
# callback env snapshot, the spacecraft's density model (a read from an
# abstract-eltype vector, so everything downstream dispatched dynamically) and
# the GRAM cache settings, and then call the scalar density model once per
# spacecraft. When every spacecraft shares one density model and nothing about
# the query depends on per-spacecraft cache state, the whole constellation is
# answered by one `getDensityBatch!` call over the planet-frame altitude,
# latitude and longitude buffers the same pre-sample has just filled.
#
# Returns that shared model, or `nothing` when any condition fails and the
# per-spacecraft path must run unchanged:
# - every spacecraft is active (the per-spacecraft path skips inactive ones and
#   leaves their buffers alone; a batch over 1:N would not);
# - every spacecraft resolves to the same density model object, exactly as
#   `_density_model_for_sat` would resolve it;
# - the model is not native GRAM (`density_model_work_is_heavy`) and not a GRAM
#   grid snapshot -- their locking and threading are not this path's to change;
# - neither freeze-per-step nor the vacuum GRAM look-ahead cache nor the GRAM
#   track cache is on, since each of those makes the value depend on
#   per-spacecraft cache state rather than on the current position;
# - GRAM runtime statistics are off (the per-spacecraft path counts calls);
# - every output buffer is at least N long, the same guard
#   `_write_density_buffers!` applies per spacecraft.
#
# Bit-identity rests on `getDensityBatch!`'s method for the model evaluating the
# same expressions as the scalar `getDensity` the per-spacecraft path calls; the
# built-in analytic models' batch methods are the scalar expressions verbatim,
# and every other model falls back to calling the scalar method per element.
# test/unit/dynamics/aero_batch_parity_tests.jl asserts it bit for bit.
function _uniform_light_density_model(p, num_sats::Int)
    num_sats >= 1 || return nothing
    callbacks = SimulationModel.SimulationCallbacks
    cb_env = callbacks._callback_env_config(p)
    cb_env.density_freeze_per_step && return nothing
    cb_env.vacuum_gram_cache_enabled && return nothing
    cb_env.gram_runtime_stats_enabled && return nothing
    sb = p.shared_buffers
    (length(sb.densities) >= num_sats && length(sb.temperatures) >= num_sats &&
     length(sb.winds) >= num_sats && length(sb.density_sample_t) >= num_sats) || return nothing
    length(p.is_active) >= num_sats || return nothing
    @inbounds for i in 1:num_sats
        p.is_active[i] || return nothing
    end
    model = callbacks._density_batch_model_for_callback(
        sb.density_models, p.args.environment_model.density_model, num_sats,
    )
    model === nothing && return nothing
    callbacks.density_model_work_is_heavy(model) && return nothing
    model isa SimulationModel.EnvironmentModels.GRAMGridAtmosphereModel && return nothing
    callbacks._gram_track_cache_enabled(cb_env.gram_track_cache, model) && return nothing
    return model
end

# One density query for the whole constellation, written where the
# per-spacecraft path would have written it.
function _fill_uniform_light_atmosphere!(
    p,
    t::Float64,
    num_sats::Int,
    model,
    alts::AbstractVector{Float64},
    lats::AbstractVector{Float64},
    lons::AbstractVector{Float64},
)::Nothing
    sb = p.shared_buffers
    SimulationModel.getDensityBatch!(
        view(sb.densities, 1:num_sats),
        view(sb.temperatures, 1:num_sats),
        view(sb.winds, 1:num_sats),
        model,
        view(alts, 1:num_sats),
        view(lats, 1:num_sats),
        view(lons, 1:num_sats),
        t,
        true,
        p,
    )
    SimulationModel.SimulationCallbacks._write_density_time_buffers!(p, num_sats, t)
    return nothing
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
    sat_idx <= length(times) || return false
    # Freeze-per-step mode: the density DiscreteCallback (get_density_callback)
    # already resamples once per accepted step; trust that value for every RHS
    # stage evaluation within the step instead of requiring an exact time match
    # (see CallbackEnvConfig.density_freeze_per_step for the rationale). Guard
    # against the pre-first-callback-firing NaN default.
    if SimulationModel.SimulationCallbacks._callback_env_config(p).density_freeze_per_step
        return isfinite(times[sat_idx])
    end
    return times[sat_idx] == t
end

@inline function sample_buffered_atmosphere(x, p, sat_idx::Int, t::Float64)::AtmosphereSample
    # Time alone does not identify a grid query: solver stages can share a time
    # and have different positions, including outside grid bounds. Concurrent
    # effectors may call this path for one satellite, so leave writes to the
    # callback/prefill owners when bypassing the buffer for a read-only grid.
    density_model = SimulationModel.SimulationCallbacks._density_model_for_sat(p, sat_idx)
    if density_model isa SimulationModel.EnvironmentModels.GRAMGridAtmosphereModel
        return sample_atmosphere(x, p, sat_idx, t; write_buffers=false)
    end
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
    # `map` over the model's own body-name tuple, not `ntuple` over its length:
    # the length is part of `NBodyGravityModel`'s type, so mapping the tuple
    # gives the compiler a concrete result type, while `ntuple(f, n::Int)` with
    # a runtime `n` does not — it boxed the closure and built the tuple
    # dynamically once per spacecraft per derivative evaluation, which is where
    # most of this path's allocation went (docs/architecture/third_body_cost.md).
    # Body order and values are unchanged; `positions_ii[k]` still belongs to
    # `model.body_names[k]`.
    positions_ii = map(model.body_names) do body_name
        body_name_spice = SimulationModel.DynamicEffectors._spice_query_name(body_name)
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

@inline function _wrench_method_available(::SimulationModel.DynamicEffectors.GravitationalHarmonicsModel)::Bool
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
