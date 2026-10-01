# Opt-in GRAM perturbed density in the dynamics.
#
# By default SpaceAGORA flies GRAM's UNPERTURBED mean density: every native path
# returns DynamicsStateC.density, and GRAM's seeded perturbation (DensityStateC
# .perturbedDensity) is computed and discarded. That perturbation is a correlated
# random walk advanced by EVERY set_position!/update! on an instance, so reading it
# inside the right-hand side is ill-posed: an adaptive integrator evaluates the
# same state repeatedly (stages, error estimates, rejected steps), each evaluation
# advances the walk, and the RHS stops being a function of (u, t).
#
# The modes below make it well-posed by keeping the walk OUT of the RHS:
#
#   step (design A)  After each ACCEPTED step, advance a separate walk instance
#                    to the accepted state (trajectory order) and hold
#                    r = perturbedDensity / density until the next accepted step,
#                    as GNC commands are held between updates.
#   pass (design B)  At the first accepted state at or below the entry interface,
#                    predict the pass drag-free (two-body + J2 RK4, the vacuum
#                    GRAM cache's predictor) until it climbs back above EI, sweep
#                    the walk in order along the prediction at a fixed time
#                    spacing, and use the piecewise-linear r(t) through those
#                    knots for the rest of the pass (last knot held if the real
#                    pass outlasts the prediction).
#   naive_rhs        NEGATIVE CONTROL ONLY: the RHS multiplies by the ratio read
#                    from the mean instance's own last update, i.e. the ill-posed
#                    in-RHS walk. For measuring that ill-posedness, never for
#                    results.
#
# In every mode the factor applies only at or below EI, where the mean density
# itself comes from GRAM; above EI the polynomial fit is used unchanged (r = 1).
# The RHS sees mean x r. The walk instance is a second native GRAM instance per
# spacecraft built from the mean model's own constructor recipe -- same seed,
# same scales -- and is touched only by this file's callback.
#
# Env (read once when the callback set is built):
#   SPACEAGORA_GRAM_DENSITY_PERTURBATION        off (default) | step | pass | naive_rhs
#   SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_DT_S   B knot spacing, s (default 1.0)
#   SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_MAX_S  B prediction cap, s (default 7200)
#   SPACEAGORA_GRAM_DENSITY_PERTURBATION_LOG    diagnostics CSV path (default: none)
#
# With the mode off nothing is installed: the SharedBuffers slot stays `nothing`
# and every density path returns exactly what it did before.

@inline function _gram_density_perturbation_mode()::Symbol
    raw = lowercase(strip(get(ENV, "SPACEAGORA_GRAM_DENSITY_PERTURBATION", "off")))
    raw in ("", "off", "0", "false", "no", "none") && return :off
    raw in ("step", "a", "per_step", "accepted_step") && return :step
    raw in ("pass", "b", "per_pass", "lookahead") && return :pass
    raw in ("naive_rhs", "naive") && return :naive_rhs
    throw(ArgumentError(
        "Unsupported SPACEAGORA_GRAM_DENSITY_PERTURBATION='$raw'. Use one of: off, step, pass, naive_rhs."
    ))
end

@inline function _gram_density_perturbation_pass_dt_s()::Float64
    v = _parse_float_env("SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_DT_S", 1.0)
    (isfinite(v) && v > 0.0) || throw(ArgumentError("SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_DT_S must be > 0, got $v"))
    return v
end

@inline function _gram_density_perturbation_pass_max_s()::Float64
    v = _parse_float_env("SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_MAX_S", 7200.0)
    (isfinite(v) && v > 0.0) || throw(ArgumentError("SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_MAX_S must be > 0, got $v"))
    return v
end

function _new_gram_density_perturbation_state(mode::Symbol, num_sats::Int, ei_m::Float64,
                                              pass_dt_s::Float64, pass_max_s::Float64,
                                              log_path::String)::GramDensityPerturbationState
    return GramDensityPerturbationState(
        mode, ei_m, pass_dt_s, pass_max_s, log_path,
        Any[nothing for _ in 1:num_sats],
        zeros(Int, num_sats),
        ones(Float64, num_sats),
        fill(false, num_sats),
        zeros(Int, num_sats),
        fill(false, num_sats),
        fill(NaN, num_sats),
        [Float64[] for _ in 1:num_sats],
        [Float64[] for _ in 1:num_sats],
        Int8[], Int[], Int[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[],
    )
end

@inline function _gram_perturbation_log!(st::GramDensityPerturbationState, kind::Integer, sat::Int,
                                         t::Float64, alt::Float64, r::Float64, sigma::Float64,
                                         mean::Float64, aux::Float64)::Nothing
    push!(st.log_kind, Int8(kind)); push!(st.log_sat, sat); push!(st.log_pass, st.pass_count[sat])
    push!(st.log_t, t); push!(st.log_alt, alt); push!(st.log_r, r); push!(st.log_sigma, sigma)
    push!(st.log_mean, mean); push!(st.log_aux, aux)
    return nothing
end

# r from one walk sample. GRAM clamps perturbedDensity to [0.1, 10] x mean, so r
# is finite and positive whenever the mean is.
@inline function _gram_ratio(perturbed::Float64, mean::Float64)::Float64
    return mean > 0.0 && isfinite(perturbed) ? perturbed / mean : 1.0
end

# ── The factor the RHS applies ───────────────────────────────────────────────

@inline function _gram_pass_factor(st::GramDensityPerturbationState, sat::Int, t::Float64)::Float64
    st.pass_active[sat] || return 1.0
    knots = st.pass_r[sat]
    n = length(knots)
    n == 0 && return 1.0
    n == 1 && return @inbounds knots[1]
    x = (t - st.pass_t0[sat]) / st.pass_dt_s
    x <= 0.0 && return @inbounds knots[1]
    x >= n - 1 && return @inbounds knots[n]
    i = floor(Int, x)
    f = x - i
    @inbounds return (1.0 - f) * knots[i + 1] + f * knots[i + 2]
end

"""
    _apply_gram_density_perturbation(p, sat_idx, t, alt, rho, T, wind)

The density every RHS-side path hands the dynamics: `rho` itself when no
perturbation mode is installed (the default -- the same value, not a product
with 1.0), otherwise `rho * r` at or below the entry interface.
"""
@inline function _apply_gram_density_perturbation(p, sat_idx::Int, t::Float64, alt::Float64,
                                                  rho::Float64, T::Float64,
                                                  wind::SVector{3, Float64})::Tuple{Float64, Float64, SVector{3, Float64}}
    st = p.shared_buffers.gram_density_perturbation[]
    st === nothing && return rho, T, wind
    alt > st.ei_m && return rho, T, wind
    sat_idx <= length(st.held_r) || return rho, T, wind
    if st.mode === :step
        return rho * (@inbounds st.held_r[sat_idx]), T, wind
    elseif st.mode === :pass
        return rho * _gram_pass_factor(st, sat_idx, t), T, wind
    elseif st.mode === :naive_rhs
        model = _density_model_for_sat(p, sat_idx)
        model isa EnvironmentModels.GRAMAtmosphereModel || return rho, T, wind
        perturbed, mean, _, _ = EnvironmentModels._gram_last_density_state(model)
        return rho * _gram_ratio(perturbed, mean), T, wind
    end
    return rho, T, wind
end

# ── Design B: drag-free pass prediction and the in-order sweep ──────────────

function _build_gram_pass!(st::GramDensityPerturbationState, p, sat::Int,
                           pos::SVector{3, Float64}, vel::SVector{3, Float64}, t0::Float64)::Nothing
    planet = p.args.environment_model.planet
    walk = st.walk_models[sat]
    dt = st.pass_dt_s
    knots = st.pass_r[sat]; alts = st.pass_alt[sat]
    empty!(knots); empty!(alts)
    n_max = max(2, ceil(Int, st.pass_max_s / dt) + 1)
    pos_k = pos; vel_k = vel
    for k in 1:n_max
        t_k = t0 + (k - 1) * dt
        l_pi = _planet_lpi_at(p, t_k)
        alt_k, lat_k, lon_k = rtolatlong(SVector{3, Float64}(l_pi * pos_k), planet)
        perturbed, mean, sigma, relstep = EnvironmentModels._gram_walk_sample(walk, alt_k, lat_k, lon_k, t_k)
        st.walk_calls[sat] += 1
        r = _gram_ratio(perturbed, mean)
        push!(knots, r); push!(alts, alt_k)
        isempty(st.log_path) || _gram_perturbation_log!(st, 2, sat, t_k, alt_k, r, sigma, mean, relstep)
        # Stop at the first knot back above EI (it closes the interpolant), or
        # if the prediction reaches the surface.
        (k > 1 && alt_k > st.ei_m) && break
        alt_k < 0.0 && break
        k < n_max && ((pos_k, vel_k) = _vacuum_rk4_step(pos_k, vel_k, planet, dt))
    end
    st.pass_t0[sat] = t0
    st.pass_active[sat] = true
    return nothing
end

@inline function _gram_pass_predicted_alt(st::GramDensityPerturbationState, sat::Int, t::Float64)::Float64
    alts = st.pass_alt[sat]
    n = length(alts)
    n == 0 && return NaN
    x = (t - st.pass_t0[sat]) / st.pass_dt_s
    x <= 0.0 && return alts[1]
    x >= n - 1 && return NaN   # outside the prediction
    i = floor(Int, x); f = x - i
    return (1.0 - f) * alts[i + 1] + f * alts[i + 2]
end

# ── The per-accepted-step callback ───────────────────────────────────────────

function _write_gram_perturbation_log(st::GramDensityPerturbationState)::Nothing
    isempty(st.log_path) && return nothing
    mkpath(dirname(abspath(st.log_path)))
    open(st.log_path, "w") do io
        println(io, "kind,sat,pass,t_s,alt_m,r,sigma_frac,mean_density_kgm3,aux")
        @inbounds for k in eachindex(st.log_kind)
            println(io, st.log_kind[k], ",", st.log_sat[k], ",", st.log_pass[k], ",",
                    repr(st.log_t[k]), ",", repr(st.log_alt[k]), ",", repr(st.log_r[k]), ",",
                    repr(st.log_sigma[k]), ",", repr(st.log_mean[k]), ",", repr(st.log_aux[k]))
        end
    end
    open(st.log_path * ".summary.toml", "w") do io
        println(io, "mode = \"", st.mode, "\"")
        println(io, "ei_m = ", repr(st.ei_m))
        println(io, "pass_dt_s = ", repr(st.pass_dt_s))
        println(io, "walk_calls = ", st.walk_calls)
        println(io, "pass_count = ", st.pass_count)
        println(io, "log_rows = ", length(st.log_kind))
    end
    return nothing
end

"""
    get_gram_density_perturbation_callback(num_sats, args) -> Union{Nothing, DiscreteCallback}

`nothing` unless `SPACEAGORA_GRAM_DENSITY_PERTURBATION` selects a mode. The
callback must run before any callback that stages density into the shared
buffers (get_callbacks places it so), and it invalidates each spacecraft's
staged sample whenever that spacecraft's factor may have changed, so no RHS
evaluation after an update reuses a density computed with the old factor.
"""
function get_gram_density_perturbation_callback(num_sats::Int, args::SimulationConfiguration)
    mode = _gram_density_perturbation_mode()
    mode === :off && return nothing
    pass_dt = _gram_density_perturbation_pass_dt_s()
    pass_max = _gram_density_perturbation_pass_max_s()
    log_path = String(strip(get(ENV, "SPACEAGORA_GRAM_DENSITY_PERTURBATION_LOG", "")))
    ei_m = Float64(args.environment_model.EI) * 1e3

    function invalidate!(p, i::Int)
        times = p.shared_buffers.density_sample_t
        i <= length(times) && (times[i] = NaN)
        return nothing
    end

    function update_sat!(st::GramDensityPerturbationState, p, u, t::Float64, i::Int)
        kin = _stage_environment_kinematics(u.sc[i], p, t)
        in_atm = kin.alt <= st.ei_m
        entered = in_atm && !st.in_atm_prev[i]
        entered && (st.pass_count[i] += 1)
        if st.mode === :step
            r_new = 1.0
            if in_atm
                perturbed, mean, sigma, relstep = EnvironmentModels._gram_walk_sample(st.walk_models[i], kin.alt, kin.lat, kin.lon, t)
                st.walk_calls[i] += 1
                r_new = _gram_ratio(perturbed, mean)
                isempty(st.log_path) || _gram_perturbation_log!(st, 1, i, t, kin.alt, r_new, sigma, mean, relstep)
            end
            if r_new != st.held_r[i]
                st.held_r[i] = r_new
                invalidate!(p, i)
            end
        elseif st.mode === :pass
            if in_atm && !st.pass_active[i]
                _build_gram_pass!(st, p, i, kin.pos_ii, kin.vel_ii, t)
                invalidate!(p, i)
            elseif !in_atm && st.pass_active[i]
                st.pass_active[i] = false
                invalidate!(p, i)
            end
            if in_atm && st.pass_active[i] && !isempty(st.log_path)
                pred = _gram_pass_predicted_alt(st, i, t)
                _gram_perturbation_log!(st, 3, i, t, kin.alt, _gram_pass_factor(st, i, t), NaN, NaN, pred - kin.alt)
            end
        elseif st.mode === :naive_rhs
            if in_atm && !isempty(st.log_path)
                model = _density_model_for_sat(p, i)
                if model isa EnvironmentModels.GRAMAtmosphereModel
                    perturbed, mean, sigma, relstep = EnvironmentModels._gram_last_density_state(model)
                    _gram_perturbation_log!(st, 4, i, t, kin.alt, _gram_ratio(perturbed, mean), sigma, mean, relstep)
                end
            end
        end
        st.in_atm_prev[i] = in_atm
        return nothing
    end

    state_ref = Ref{Union{Nothing, GramDensityPerturbationState}}(nothing)

    function initialize(cb, u, t, integrator)
        p = integrator.p
        # A fresh state on every initialization: a solver-policy re-solve from
        # the start then replays a fresh walk instead of continuing a used one.
        st = _new_gram_density_perturbation_state(mode, num_sats, ei_m, pass_dt, pass_max, log_path)
        if mode === :step || mode === :pass
            for i in 1:num_sats
                model = _density_model_for_sat(p, i)
                model isa EnvironmentModels.GRAMAtmosphereModel || throw(ArgumentError(
                    "SPACEAGORA_GRAM_DENSITY_PERTURBATION=$(mode) needs a native GRAMAtmosphereModel " *
                    "for spacecraft $i; got $(typeof(model))."
                ))
                st.walk_models[i] = EnvironmentModels._gram_walk_clone(model)
            end
        end
        state_ref[] = st
        p.shared_buffers.gram_density_perturbation[] = st
        for i in 1:num_sats
            p.is_active[i] && update_sat!(st, p, u, Float64(t), i)
        end
        return nothing
    end

    function affect!(integrator)
        st = state_ref[]
        st === nothing && return nothing
        p = integrator.p
        u = integrator.u
        t = Float64(integrator.t)
        for i in 1:num_sats
            p.is_active[i] && update_sat!(st, p, u, t, i)
        end
        return nothing
    end

    function finalize(cb, u, t, integrator)
        st = state_ref[]
        st === nothing || _write_gram_perturbation_log(st)
        return nothing
    end

    condition(u, t, integrator) = true
    return DiscreteCallback(condition, affect!; initialize=initialize, finalize=finalize,
                            save_positions=(false, false))
end
