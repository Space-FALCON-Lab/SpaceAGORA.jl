"""
    SpaceAGORA.ContactLoadServices

Experimental v0.2.0 contract for same-trial-state contact loads. Core owns all
body states and integration; providers return contact COM wrenches and optional
accepted-transition velocity corrections. Each run owns a separate runtime,
closed on success and failure. Providers reset trial scratch on every query,
commit hybrid modes only through accepted transitions, and never integrate time.
"""
module ContactLoadServices
using LinearAlgebra, StaticArrays

const CONTRACT_VERSION = v"0.2.0"
struct CompatibilityError <: Exception
    message::String
end
Base.showerror(io::IO, e::CompatibilityError) = print(io, e.message)
"""Require the exact experimental contact contract before defining provider methods."""
function require_version(expected::VersionNumber)
    expected == CONTRACT_VERSION || throw(CompatibilityError(
        "Contact adapter requires $expected; this build provides $CONTRACT_VERSION."))
    nothing
end

"""Immutable provider configuration; it must contain no native resource or run state."""
abstract type AbstractContactModel end
"""One run's provider state. Reset trial scratch on every evaluate call; never integrate."""
abstract type AbstractContactRuntime end

"""Contact model whose discrete contact modes change only at accepted core transitions."""
abstract type AbstractHybridContactModel <: AbstractContactModel end
"""Per-run hybrid contact state; trial evaluations must not commit mode changes."""
abstract type AbstractHybridContactRuntime <: AbstractContactRuntime end

"""Numerical velocity projection at an accepted transition, in world linear
and body angular coordinates. This is recorded separately from contact work."""
struct ContactVelocityCorrection
    linear_world_m_s::SVector{3,Float64}
    angular_body_rad_s::SVector{3,Float64}
end

"""Positive inside the current contact mode; a downward zero crossing requests
an accepted transition. Must be a pure query at the supplied complete state."""
contact_event_value(::AbstractHybridContactRuntime,states,time_s::Float64)=
    throw(ArgumentError("Hybrid provider does not implement contact_event_value."))
"""Commit contact-mode changes at an accepted core state and return one numerical
velocity correction per body. Called at initialization, localized events and
accepted steps; never during rejected stages. Corrections must conserve total
linear/angular momentum and cannot add kinetic energy. Trials must not mutate
the committed mode, and this hook must recompute from the supplied state."""
contact_accept!(::AbstractHybridContactRuntime,states,time_s::Float64)=
    throw(ArgumentError("Hybrid provider does not implement contact_accept!."))

"""Physical body at its COM. Slot 1 is the rigid bus; remaining slots are the
attachment or compliant-arm links in state-column order."""
struct ContactBody
    mass_kg::Float64
    inertia_body_kg_m2::SMatrix{3,3,Float64,9}
    function ContactBody(mass, inertia)
        m = Float64(mass); Ibody = SMatrix{3,3,Float64,9}(inertia)
        isfinite(m) && m > 0 && all(isfinite, Ibody) && issymmetric(Ibody) && isposdef(Ibody) ||
            throw(ArgumentError("Contact bodies require positive finite mass and symmetric positive-definite inertia."))
        new(m, Ibody)
    end
end

"""Read-only trial observation: SI, COM position/velocity in inertial-parallel axes,
scalar-last active body-to-world quaternion, body angular velocity, and all smooth
physical COM force/torque in world axes. A fixed or constant-velocity origin is
used for positions/velocities. Smooth torques exclude gyroscopic pseudo-torques.
"""
struct ContactState
    position_m::SVector{3,Float64}
    velocity_m_s::SVector{3,Float64}
    quaternion::SVector{4,Float64}
    angular_velocity_body_rad_s::SVector{3,Float64}
    smooth_force_n::SVector{3,Float64}
    smooth_torque_nm::SVector{3,Float64}
end

"""Contact-only force and torque about this body's COM, both in world axes (N, N m)."""
struct ContactWrench
    force_n::SVector{3,Float64}
    torque_nm::SVector{3,Float64}
end

"""Validate provider geometry/settings before allocating native resources."""
contact_preflight(model::AbstractContactModel, bodies::NTuple{N,ContactBody}) where {N} =
    throw(ArgumentError("No contact provider implements preflight for $(typeof(model))."))
"""Allocate one runtime; clean up partial allocations if preparation throws."""
contact_prepare(model::AbstractContactModel, bodies::NTuple{N,ContactBody}) where {N} =
    throw(ArgumentError("No contact provider implements prepare for $(typeof(model))."))
"""Return one ContactWrench per body for this complete trial state/time.
Reject nonfinite results, warnings, capacity exhaustion and unconverged solves;
never reuse a previous force or advance native time. Inputs cannot be mutated.
"""
contact_evaluate!(runtime::AbstractContactRuntime, states::NTuple{N,ContactState}, time_s::Float64) where {N} =
    throw(ArgumentError("No contact evaluator implements $(typeof(runtime))."))
"""Release this run's native resources. Must be safe after a failed query and idempotent."""
contact_close!(runtime::AbstractContactRuntime) =
    throw(ArgumentError("No contact cleanup implements $(typeof(runtime))."))
"""Return a shallow configuration copy with the one experimental contact model installed."""
function with_contact_model end

# Internal ownership record; deliberately absent from the provider inventory.
struct ContactRunState
    runtime::AbstractContactRuntime
    origin_m::SVector{3,Float64}
    origin_velocity_m_s::SVector{3,Float64}
    busy::Threads.Atomic{Bool}
    observations::Base.RefValue{Any}
    projection_work_J::Base.RefValue{Float64}
    projection_count::Base.RefValue{Int}
    projections::Vector{NamedTuple{(:time_s,:work_J),Tuple{Float64,Float64}}}
end
const IMPLEMENTED_FUNCTIONS = (:contact_preflight, :contact_prepare, :contact_evaluate!, :contact_close!,
    :contact_event_value, :contact_accept!)
const CONSUMED_SERVICES = (:AbstractContactModel, :AbstractContactRuntime, :ContactBody,
    :ContactState, :ContactWrench, :with_contact_model, :AbstractHybridContactModel,
    :AbstractHybridContactRuntime, :ContactVelocityCorrection)
export AbstractContactModel, AbstractContactRuntime, ContactBody, ContactState, ContactWrench
export contact_preflight, contact_prepare, contact_evaluate!, contact_close!, with_contact_model
export AbstractHybridContactModel, AbstractHybridContactRuntime, ContactVelocityCorrection
export contact_event_value, contact_accept!
end
