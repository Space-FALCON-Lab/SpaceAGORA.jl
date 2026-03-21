@kwdef struct TimeTerminalCondition <: AbstractTerminalCondition
    time_limit::Float64
end

@kwdef struct OrbitTerminalCondition <: AbstractTerminalCondition
    num_orbits::Int
end

@inline function terminal_condition_value(
    tc::TimeTerminalCondition,
    u,
    time::Float64,
    integrator;
    idx::Int=1
)
    return tc.time_limit - time
end

@inline function terminal_condition_value(
    tc::OrbitTerminalCondition,
    u,
    time::Float64,
    integrator;
    idx::Int=1
)
    return Float64(tc.num_orbits - integrator.p.orbit_counter[idx])
end

@inline function check_terminal_condition(tc::AbstractTerminalCondition, u, time::Float64, integrator; idx::Int=1)::Bool
    return terminal_condition_value(tc, u, time, integrator; idx=idx) <= 0.0
end

@inline function estimated_terminal_time(tc::TimeTerminalCondition, u0, args)::Float64
    return tc.time_limit
end

function estimated_terminal_time(tc::OrbitTerminalCondition, u0, args)::Float64
    planet = args.environment_model.planet
    max_period_s = 0.0

    @inbounds for i in eachindex(u0.sc)
        pos = SimulationModel.SVector{3, Float64}(u0.sc[i].pos)
        vel = SimulationModel.SVector{3, Float64}(u0.sc[i].vel)
        oe = try
            SimulationModel.rvtoorbitalelement(pos, vel, planet)
        catch
            nothing
        end
        oe === nothing && continue

        a = Float64(oe[1])
        e = Float64(oe[2])
        if !isfinite(a) || !isfinite(e) || a <= 0.0 || e < 0.0 || e >= 1.0
            continue
        end

        n = sqrt(planet.μ / a^3)
        if !isfinite(n) || n <= 0.0
            continue
        end

        period_s = 2pi / n
        if isfinite(period_s) && period_s > max_period_s
            max_period_s = period_s
        end
    end

    if !(isfinite(max_period_s) && max_period_s > 0.0)
        max_period_s = 24.0 * 3600.0
    end

    return max(1.0, Float64(tc.num_orbits) * max_period_s * 1.25)
end
