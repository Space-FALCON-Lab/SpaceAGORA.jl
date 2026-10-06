using ..DynamicEffectors.AerodynamicEffectors: _thermal_incidence_mode, _aero_link_angles

@inline function _heat_rate_buffer_for_sat!(p, sat_idx::Int)
    links = p.args.dynamics_model.spacecraft[sat_idx].links
    n_links = length(links)
    heat_rates = p.shared_buffers.heat_rates[sat_idx]
    if length(heat_rates) != n_links
        resize!(heat_rates, n_links)
    end
    fill!(heat_rates, 0.0)
    return heat_rates
end

function _compute_stage_heat_rates!(
    p,
    x,
    sat_idx::Int,
    t::Float64;
    use_buffered_density::Bool=false,
)
    links = p.args.dynamics_model.spacecraft[sat_idx].links
    isempty(links) && return _heat_rate_buffer_for_sat!(p, sat_idx)
    # Vacuum short-circuit. Without an atmosphere there is no aerothermal
    # heating, so this call can only ever fall out of the `rho <= 0` guard
    # below with an all-zero buffer -- but not before paying for a full
    # sample_planet_frame (rtolatlong: atan/asin/sqrt per satellite) plus a
    # density sample, once per satellite per RHS stage. Profiling a
    # 1024-satellite L20 vacuum constellation put that at ~41% of the entire
    # solve. The type check is static, so this compiles away for every
    # configuration that does have an atmosphere.
    if p.args.environment_model.density_model isa NoAtmosphereModel
        return _heat_rate_buffer_for_sat!(p, sat_idx)
    end
    heat_rates = _heat_rate_buffer_for_sat!(p, sat_idx)
    engine = _simulation_engine_module()
    planet_frame = engine.sample_planet_frame(x, p, sat_idx, t)
    atmosphere = use_buffered_density ?
        engine.sample_buffered_atmosphere(x, p, sat_idx, t) :
        engine.sample_atmosphere(x, p, sat_idx, t; write_buffers=false)

    rho = atmosphere.rho_kg_m3
    T = atmosphere.temperature_k
    wind = atmosphere.wind_pp
    if !isfinite(rho) || !isfinite(T) || rho <= 0.0 || T <= 0.0
        return heat_rates
    end

    planet = p.args.environment_model.planet
    thermal_model = p.args.environment_model.thermal_model
    uD, uN, uE = latlongtoNED((planet_frame.alt_m, planet_frame.lat_rad, planet_frame.lon_rad))
    wE, wN, wU = wind
    wind_pp = wN * uN + wE * uE - wU * uD
    vel_pp_rw = planet_frame.vel_pp - wind_pp
    v = norm(vel_pp_rw)
    sound_velocity = sqrt(planet.γ * planet.R * T)
    if !isfinite(v) || !isfinite(sound_velocity) || v <= 0.0 || sound_velocity <= 0.0
        return heat_rates
    end

    mach = v / sound_velocity
    S = sqrt(planet.γ * 0.5) * mach
    orientation_sim = p.args.mission_configuration.orientation_sim
    incidence_mode = _thermal_incidence_mode(p.args.dynamics_model.dynamic_effectors, orientation_sim)
    spacecraft = p.args.dynamics_model.spacecraft[sat_idx]
    q_root = incidence_mode !== nothing && orientation_sim ?
        engine.build_state_sample(x, spacecraft, true).q_ib : nothing
    vel_pi = orientation_sim ? planet_frame.l_pi' * vel_pp_rw : SVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for j in eachindex(links)
        # Typed wrenches intentionally do not mutate Link angles. Derive heating
        # from the same current geometry, even if thermal sampling runs first.
        # No recognized aero model: preserve custom/legacy stored-angle inputs.
        alpha = if incidence_mode === nothing
            links[j].α
        else
            angle, _, _ = _aero_link_angles(spacecraft, links[j], 1,
                orientation_sim, vel_pi, 0.0, incidence_mode, q_root)
            angle
        end
        if !isfinite(alpha)
            continue
        end
        qdot = getHeatRate(thermal_model, S, T, rho, v, alpha)
        heat_rates[j] = (isfinite(qdot) && qdot > 0.0) ? qdot : 0.0
    end
    return heat_rates
end

# Per-satellite cost class for the thermal callback's thread decision.
#
# `_compute_stage_heat_rates!` is called here with `use_buffered_density=true`
# (it reads shared_buffers rather than evaluating a density model), so its
# per-satellite cost is one `sample_planet_frame` (a planet-relative
# position/velocity transform through `rtolatlong`) plus one `getHeatRate`
# call per thermal link on that spacecraft (a Maxwellian free-molecular
# heating model -- a handful of `erf`/`exp`/`sqrt` evaluations). The
# `sample_planet_frame` part is fixed per satellite; the per-link loop is what
# grows, and it is what a fanned-out dispatch actually has more of to hand a
# worker.
#
# Measured driving _compute_stage_heat_rates! directly for an 8-spacecraft
# shape (mirrors the density-guard measurement shape), median of 9 timed
# batches with GC disabled during each batch, at thread widths 1 (serial), 2,
# 4 and 8, varying links-per-spacecraft. `_thread_worker_count` caps the
# worker count at num_sats, so width 8 -- one worker per satellite -- is what
# an 8-satellite auto dispatch actually reaches under any budget >= 8
# (the shape this callback was seen regressing on used 24 threads); width 2
# is the achieved width only under a much narrower budget.
#
# us/call, serial vs. width 8, by links-per-spacecraft:
#   1: 45.71 / 48.21   2: 48.84 / 49.67    4: 46.72 / 50.38
#   8: 49.44 / 49.46   16: 51.74 / 50.97   32: 61.62 / 51.57
#   64: 78.11 / 54.90  128: 109.83 / 59.51  256: 170.64 / 72.19
#   512: 298.93 / 93.09
#
# At and below 8 links, width 8 is a dead heat with serial or slightly worse
# (dispatch overhead not amortized). 16 links is the first point where
# threaded is measurably faster (-1.5%), and every larger count measured
# widens that margin monotonically (-16% at 32, -69% at 512) -- consistent
# with a real crossover rather than noise. Width 2 does not reach break-even
# until roughly 64 links, so a narrow-budget run stays conservative near the
# threshold; that is judged an acceptable trade against leaving every
# ordinary few-link vehicle threaded for no reason.  See
# test/unit/parallel/thermal_callback_light_work_tests.jl for the guard this
# backs.
const THERMAL_CALLBACK_HEAVY_LINK_THRESHOLD = 16  # measured, see comment above

"""
    _thermal_callback_work_is_heavy(p, num_sats) -> Bool

True when at least one of the first `num_sats` spacecraft has enough thermal
links that the per-satellite body of `_compute_stage_heat_rates!` is worth a
threaded dispatch. `p === nothing` (the no-run-state overload of
`_thermal_callback_thread_decision`) answers `true` so that call path is
unaffected -- it exists only for direct unit testing, not for a real run.
"""
@inline function _thermal_callback_work_is_heavy(p, num_sats::Int)::Bool
    p === nothing && return true
    spacecraft = p.args.dynamics_model.spacecraft
    limit = min(num_sats, length(spacecraft))
    limit <= 0 && return false
    @inbounds for i in 1:limit
        length(spacecraft[i].links) >= THERMAL_CALLBACK_HEAVY_LINK_THRESHOLD && return true
    end
    return false
end

function get_thermal_callback(num_sats::Int, args::SimulationConfiguration)
    function update_thermal_sat!(i::Int, p, u, t::Float64)
        _compute_stage_heat_rates!(p, u.sc[i], i, t; use_buffered_density=true)
        return nothing
    end

    condition(u, t, integrator) = true

    function affect!(integrator)
        p = integrator.p
        u = integrator.u
        decision = _thermal_callback_thread_decision(
            p, num_sats;
            heavy_work=_thermal_callback_work_is_heavy(p, num_sats)
        )
        use_threads = decision.use_threads
        started_ns = time_ns()
        if use_threads
            ParallelPolicy.threaded_foreach_persistent(:thermal_callback, num_sats, decision.allotment) do i
                @inbounds update_thermal_sat!(i, p, u, Float64(integrator.t))
            end
        else
            @inbounds for i in 1:num_sats
                update_thermal_sat!(i, p, u, Float64(integrator.t))
            end
        end
        if decision.policy_applied
            ParallelPolicy.record_policy_observation!(
                :thermal_callback;
                mode=decision.mode,
                num_items=num_sats,
                use_threads=use_threads,
                elapsed_ns=(time_ns() - started_ns),
                env=_policy_env_config(p),
                ctx=ParallelPolicy.policy_context_hint(p)
            )
        end
    end

    # Housekeeping, not an event: the affect refreshes the thermal samples in
    # p and never touches u, so the before/after saves of the DiscreteCallback
    # default would only append two more copies of every accepted step.
    return DiscreteCallback(condition, affect!; initialize=(cb, u, t, integrator) -> affect!(integrator),
        save_positions=(false, false))
end
