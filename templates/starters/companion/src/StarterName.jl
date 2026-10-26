"""
    StarterName

A SpaceAGORA companion that serves the external-propagation boundary
(`SpaceAGORA.ExternalPropagationServices`): spacecraft that an outside integrator
owns for a whole run. `CircularOrbitPropagator` is an example owner to replace
with the companion's own engine. The companion requires the exact contract
version of the facade before it defines any method, and again at load.
"""
module StarterName

import SpaceAGORA
using LinearAlgebra: cross, dot, norm
using StaticArrays: SVector

isdefined(SpaceAGORA, :ExternalPropagationServices) || throw(ArgumentError(
    "This SpaceAGORA revision has no ExternalPropagationServices facade; use the commit in SpaceAGORASource.toml."))
const EP = SpaceAGORA.ExternalPropagationServices
const CORE_CONTRACT_VERSION = v"0.0.0"
EP.require_version(CORE_CONTRACT_VERSION)

# Recheck when a precompiled companion loads against another core checkout.
__init__() = EP.require_version(CORE_CONTRACT_VERSION)

export CircularOrbitPropagator

"""
    CircularOrbitPropagator(spacecraft; step_s=0.5, mu=3.986004418e14)

Owns the spacecraft at the run indices `spacecraft` and moves each on the
circular orbit through its initial position and velocity about a point mass of
gravitational parameter `mu` (m^3/s^2), in whole steps of `step_s` seconds. The
attitude stays at its initial value. An example to replace with the companion's
own engine.
"""
struct CircularOrbitPropagator <: EP.AbstractExternalPropagator
    spacecraft::Vector{Int}
    step_s::Float64
    mu::Float64
end
CircularOrbitPropagator(spacecraft; step_s=0.5, mu=3.986004418e14) =
    CircularOrbitPropagator(collect(Int, spacecraft), Float64(step_s), Float64(mu))

# One owned spacecraft's circle: radius, rate, in-plane unit vectors, attitude.
struct Circle
    radius::Float64
    rate::Float64
    e1::SVector{3, Float64}
    e2::SVector{3, Float64}
    q::SVector{4, Float64}
end

# The run's own state: built from u0 in external_prepare, never shared between runs.
mutable struct CircleRuntime
    step_s::Float64
    steps::Int
    circles::Vector{Circle}
end

function circle(initial, mu)
    radius = norm(initial.pos)
    e1 = initial.pos / radius
    normal = cross(initial.pos, initial.vel)
    e2 = cross(normal, initial.pos) / norm(cross(normal, initial.pos))
    return Circle(radius, sqrt(mu / radius^3), e1, e2, initial.q)
end

position(c::Circle, t) = c.radius * (cos(c.rate * t) * c.e1 + sin(c.rate * t) * c.e2)
velocity(c::Circle, t) = c.radius * c.rate * (-sin(c.rate * t) * c.e1 + cos(c.rate * t) * c.e2)

EP.external_spacecraft(owner::CircularOrbitPropagator) = copy(owner.spacecraft)
EP.external_step(owner::CircularOrbitPropagator) = owner.step_s

function EP.external_preflight(owner::CircularOrbitPropagator, args)
    owner.step_s > 0 || throw(ArgumentError("CircularOrbitPropagator: step_s must be positive."))
    owner.mu > 0 || throw(ArgumentError("CircularOrbitPropagator: mu must be positive."))
    return nothing
end

EP.external_prepare(owner::CircularOrbitPropagator, args, u0) = CircleRuntime(owner.step_s, 0,
    [circle(EP.external_initial_state(u0, index), owner.mu) for index in owner.spacecraft])

function EP.external_sync!(runtime::CircleRuntime, t::Float64)
    target = floor(Int, t / runtime.step_s + 1e-9)
    taken = max(target - runtime.steps, 0)
    runtime.steps += taken
    return taken
end

EP.external_time(runtime::CircleRuntime) = runtime.steps * runtime.step_s

# The analytic orbit is exact at any time, so the lag after the last whole step
# needs no separate model.
function EP.external_state(runtime::CircleRuntime, t::Float64, slot::Int)
    c = runtime.circles[slot]
    return (pos=position(c, t), vel=velocity(c, t), q=c.q, ω=SVector(0.0, 0.0, 0.0))
end

EP.external_acceleration(runtime::CircleRuntime, t::Float64, slot::Int) =
    -runtime.circles[slot].rate^2 * position(runtime.circles[slot], t)

end # module StarterName
