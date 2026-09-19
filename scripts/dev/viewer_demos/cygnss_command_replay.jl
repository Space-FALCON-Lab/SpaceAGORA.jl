"""Recorded wheel-command reconstruction. This is a development tool, not an actuator model.

Columns are signed axial wheel torques in N m, with a caller-supplied channel gain.
A piecewise-linear interpolation uses actual packet times; its exact integral
advances wheel momentum from one measured initial condition. The body receives
`-A*tau_w - omega × H_w`. Do not also install ReactionWheelMomentumModel: that
would apply the same wheel reaction twice. The engine owns `-omega × I*omega`.
Rod and environmental torques are separate. A control demand is not necessarily
an executed command; the loader requires the caller to name the channel explicitly.
"""
module CygnssCommandReplay
using SpaceAGORA, StaticArrays, LinearAlgebra, Arrow, DataFrames
const SM = SpaceAGORA.SimulationModel
export CommandedWheelReplay, command_values, load_wheel_commands

struct CommandedWheelReplay <: SM.AbstractForceTorqueModel
    times::Vector{Float64}
    torque::Matrix{Float64}
    momentum::Matrix{Float64}
    axes::Matrix{Float64}
    offset::Float64
    spacecraft_index::Int
end

function CommandedWheelReplay(times, torque, axes, initial_momentum;
                              time_offset_s=0.0, spacecraft_index=1, gain=1.0)
    t = Float64.(times); u = Float64.(torque); A = Float64.(axes)
    n = length(t); nw = length(initial_momentum)
    n >= 2 && all(isfinite, t) && all(>(0), diff(t)) ||
        throw(ArgumentError("command times must be finite and strictly increasing"))
    size(u) == (n, nw) && size(A) == (3, nw) ||
        throw(ArgumentError("commands must be N by wheels and axes 3 by wheels"))
    all(isfinite, u) && all(isfinite, A) && all(isfinite, initial_momentum) &&
        isfinite(gain) && isfinite(time_offset_s) && spacecraft_index >= 1 ||
        throw(ArgumentError("non-finite command inputs or invalid spacecraft index"))
    u .*= gain
    h = Matrix{Float64}(undef, n, nw); h[1, :] .= initial_momentum
    for i in 2:n
        h[i, :] .= h[i-1, :] .+ (t[i]-t[i-1])/2 .* (u[i-1, :] .+ u[i, :])
    end
    return CommandedWheelReplay(t, u, h, A, Float64(time_offset_s), spacecraft_index)
end

"Body momentum and its rate at run time t; no extrapolation past the supplied record."
function command_values(m::CommandedWheelReplay, t::Real)
    tt = Float64(t) + m.offset
    first(m.times) <= tt <= last(m.times) || throw(DomainError(tt, "outside command record"))
    k = min(searchsortedlast(m.times, tt), length(m.times)-1)
    dt = tt - m.times[k]; span = m.times[k+1] - m.times[k]
    u0 = @view m.torque[k, :]; u1 = @view m.torque[k+1, :]
    du = (u1-u0)/span
    u = u0 + dt*du
    h = m.momentum[k, :] + dt*u0 + (dt^2/2)*du
    return (momentum=SVector{3,Float64}(m.axes*h), rate=SVector{3,Float64}(m.axes*u), axial_momentum=h)
end

SM.environment_requirements(::CommandedWheelReplay) = SM.EffectorEnvironmentRequirements()
function SM.wrench(m::CommandedWheelReplay, x::SM.StateSample, env::SM.EnvironmentSample, t::Float64)
    z = SVector{3,Float64}(0,0,0)
    x.ω_body === nothing && return (z,z)
    v = command_values(m,t)
    return z, SM.wheel_reaction_torque(v.momentum, v.rate, x.ω_body)
end
function SM.wrench_caching!(m::CommandedWheelReplay, x::SM.StateSample, env::SM.EnvironmentSample,
                            t::Float64, p, sat_idx::Int)
    z = SVector{3,Float64}(0,0,0)
    sat_idx == m.spacecraft_index || return (z,z)
    return SM.wrench(m,x,env,t)
end

"Read an explicit command channel on its own original packet-time grid."
function load_wheel_commands(path::AbstractString; channel::Symbol)
    channel in (:rwCmd, :tqDmdCtrl) || throw(ArgumentError("channel must be rwCmd or tqDmdCtrl"))
    df = DataFrame(Arrow.Table(path))
    cols = ["$(channel)_$k" for k in 0:2]
    all(c -> c in names(df), ["t_rel"; cols]) || throw(ArgumentError("missing $channel columns in $path"))
    t = Float64.(df.t_rel); u = hcat((Float64.(df[!,c]) for c in cols)...)
    length(t) >= 2 && all(isfinite,t) && all(>(0),diff(t)) && all(isfinite,u) ||
        throw(ArgumentError("invalid command time series in $path"))
    return (times=t, torque=u, channel=channel)
end
end
