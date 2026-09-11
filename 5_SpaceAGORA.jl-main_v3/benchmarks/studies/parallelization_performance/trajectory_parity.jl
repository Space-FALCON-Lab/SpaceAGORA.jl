Base.@kwdef struct PPCParityThresholds
    pos_rms_rel::Float64 = 1e-8
    pos_max_rel::Float64 = 1e-6
    vel_rms_rel::Float64 = 1e-8
    vel_max_rel::Float64 = 1e-6
    q_angle_max_rad::Float64 = 1e-8
    omega_max_rel::Float64 = 1e-7
    mass_max_rel::Float64 = 1e-10
    event_time_abs_s::Float64 = 1e-6
end

@inline function _ppc_vec3(state, field::Symbol)
    hasproperty(state, field) || return nothing
    value = getproperty(state, field)
    try
        return SVector{3, Float64}(value)
    catch
        return nothing
    end
end

@inline function _ppc_vec4(state, field::Symbol)
    hasproperty(state, field) || return nothing
    value = getproperty(state, field)
    try
        return SVector{4, Float64}(value)
    catch
        return nothing
    end
end

@inline function _ppc_mass(state)
    for field in (:mass, :m)
        if hasproperty(state, field)
            value = try
                Float64(getproperty(state, field))
            catch
                NaN
            end
            isfinite(value) && return value
        end
    end
    return missing
end

function _ppc_state_at(sol, t::Float64)
    isempty(sol.t) && return nothing
    idx = searchsortedfirst(sol.t, t)
    if idx <= 1
        return sol.u[1]
    elseif idx > length(sol.t)
        return sol.u[end]
    end
    before = idx - 1
    after = idx
    return abs(sol.t[before] - t) <= abs(sol.t[after] - t) ? sol.u[before] : sol.u[after]
end

@inline function _ppc_rel_error(num::Float64, den::Float64)::Float64
    scale = max(abs(den), eps(Float64))
    return abs(num) / scale
end

@inline function _ppc_quat_angle(q_ref::SVector{4, Float64}, q_cmp::SVector{4, Float64})::Float64
    dot_abs = abs(dot(q_ref / norm(q_ref), q_cmp / norm(q_cmp)))
    return 2.0 * acos(clamp(dot_abs, -1.0, 1.0))
end

function _ppc_zero_crossings(signal::Vector{Float64}, times::Vector{Float64}; direction::Symbol)::Vector{Float64}
    crossings = Float64[]
    for i in 2:min(length(signal), length(times))
        y0 = signal[i - 1]
        y1 = signal[i]
        if direction == :up && y0 < 0.0 && y1 >= 0.0
            push!(crossings, times[i - 1] + (0.0 - y0) * (times[i] - times[i - 1]) / (y1 - y0))
        elseif direction == :down && y0 > 0.0 && y1 <= 0.0
            push!(crossings, times[i - 1] + (0.0 - y0) * (times[i] - times[i - 1]) / (y1 - y0))
        end
    end
    return crossings
end

function ppc_event_times(sol, args)::NamedTuple
    if isempty(sol.t) || isempty(sol.u)
        return (periapsis=Float64[], interface=Float64[])
    end
    radial = Float64[]
    interface = Float64[]
    planet_radius = args.environment_model.planet.Rp_e
    interface_alt = args.environment_model.EI * 1e3
    for u in sol.u
        sc = u.sc[1]
        pos = _ppc_vec3(sc, :pos)
        vel = _ppc_vec3(sc, :vel)
        if pos === nothing || vel === nothing
            push!(radial, NaN)
            push!(interface, NaN)
        else
            push!(radial, dot(pos, vel))
            push!(interface, norm(pos) - planet_radius - interface_alt)
        end
    end
    times = Float64.(sol.t)
    return (
        periapsis=_ppc_zero_crossings(radial, times; direction=:up),
        interface=_ppc_zero_crossings(interface, times; direction=:down)
    )
end

function ppc_compare_trajectories(sol_ref, sol_cmp, args; sample_count::Int=128)::NamedTuple
    isempty(sol_ref.t) && error("Reference trajectory has no saved times.")
    isempty(sol_cmp.t) && error("Candidate trajectory has no saved times.")
    t0 = max(Float64(first(sol_ref.t)), Float64(first(sol_cmp.t)))
    tf = min(Float64(last(sol_ref.t)), Float64(last(sol_cmp.t)))
    tf >= t0 || error("No overlapping trajectory time interval.")
    samples = sample_count <= 1 ? [t0] : collect(range(t0, tf; length=sample_count))
    n_sats = length(args.dynamics_model.spacecraft)
    pos_rel = Float64[]
    vel_rel = Float64[]
    q_angles = Float64[]
    omega_rel = Float64[]
    mass_rel = Float64[]

    for t in samples
        u_ref = _ppc_state_at(sol_ref, t)
        u_cmp = _ppc_state_at(sol_cmp, t)
        u_ref === nothing && continue
        u_cmp === nothing && continue
        for sat_idx in 1:n_sats
            sat_idx <= length(u_ref.sc) || continue
            sat_idx <= length(u_cmp.sc) || continue
            s_ref = u_ref.sc[sat_idx]
            s_cmp = u_cmp.sc[sat_idx]
            p_ref = _ppc_vec3(s_ref, :pos)
            p_cmp = _ppc_vec3(s_cmp, :pos)
            v_ref = _ppc_vec3(s_ref, :vel)
            v_cmp = _ppc_vec3(s_cmp, :vel)
            if p_ref !== nothing && p_cmp !== nothing
                push!(pos_rel, _ppc_rel_error(norm(p_cmp - p_ref), norm(p_ref)))
            end
            if v_ref !== nothing && v_cmp !== nothing
                push!(vel_rel, _ppc_rel_error(norm(v_cmp - v_ref), norm(v_ref)))
            end
            q_ref = _ppc_vec4(s_ref, :q)
            q_cmp = _ppc_vec4(s_cmp, :q)
            if q_ref !== nothing && q_cmp !== nothing
                push!(q_angles, _ppc_quat_angle(q_ref, q_cmp))
            end
            w_ref = _ppc_vec3(s_ref, :ω)
            w_cmp = _ppc_vec3(s_cmp, :ω)
            if w_ref !== nothing && w_cmp !== nothing
                push!(omega_rel, _ppc_rel_error(norm(w_cmp - w_ref), norm(w_ref)))
            end
            m_ref = _ppc_mass(s_ref)
            m_cmp = _ppc_mass(s_cmp)
            if !(m_ref isa Missing) && !(m_cmp isa Missing)
                push!(mass_rel, _ppc_rel_error(Float64(m_cmp) - Float64(m_ref), Float64(m_ref)))
            end
        end
    end

    ref_events = ppc_event_times(sol_ref, args)
    cmp_events = ppc_event_times(sol_cmp, args)
    event_count_equal = length(ref_events.periapsis) == length(cmp_events.periapsis) &&
        length(ref_events.interface) == length(cmp_events.interface)
    event_deltas = Float64[]
    for (a, b) in zip(ref_events.periapsis, cmp_events.periapsis)
        push!(event_deltas, abs(a - b))
    end
    for (a, b) in zip(ref_events.interface, cmp_events.interface)
        push!(event_deltas, abs(a - b))
    end

    thresholds = PPCParityThresholds()
    pos_rms = isempty(pos_rel) ? 0.0 : sqrt(mean(abs2, pos_rel))
    vel_rms = isempty(vel_rel) ? 0.0 : sqrt(mean(abs2, vel_rel))
    pos_max = isempty(pos_rel) ? 0.0 : maximum(pos_rel)
    vel_max = isempty(vel_rel) ? 0.0 : maximum(vel_rel)
    q_max = isempty(q_angles) ? 0.0 : maximum(q_angles)
    omega_max = isempty(omega_rel) ? 0.0 : maximum(omega_rel)
    mass_max = isempty(mass_rel) ? 0.0 : maximum(mass_rel)
    event_time_max = isempty(event_deltas) ? 0.0 : maximum(event_deltas)
    pass = pos_rms <= thresholds.pos_rms_rel &&
        pos_max <= thresholds.pos_max_rel &&
        vel_rms <= thresholds.vel_rms_rel &&
        vel_max <= thresholds.vel_max_rel &&
        q_max <= thresholds.q_angle_max_rad &&
        omega_max <= thresholds.omega_max_rel &&
        mass_max <= thresholds.mass_max_rel &&
        event_count_equal &&
        event_time_max <= thresholds.event_time_abs_s

    return (
        pass=pass,
        samples=length(samples),
        pos_rel_rms=pos_rms,
        pos_rel_p90=isempty(pos_rel) ? 0.0 : quantile(pos_rel, 0.9),
        pos_rel_max=pos_max,
        vel_rel_rms=vel_rms,
        vel_rel_p90=isempty(vel_rel) ? 0.0 : quantile(vel_rel, 0.9),
        vel_rel_max=vel_max,
        q_angle_max_rad=q_max,
        omega_rel_max=omega_max,
        mass_rel_max=mass_max,
        event_count_equal=event_count_equal,
        event_time_abs_max_s=event_time_max,
        ref_periapsis_count=length(ref_events.periapsis),
        cmp_periapsis_count=length(cmp_events.periapsis),
        ref_interface_count=length(ref_events.interface),
        cmp_interface_count=length(cmp_events.interface)
    )
end
