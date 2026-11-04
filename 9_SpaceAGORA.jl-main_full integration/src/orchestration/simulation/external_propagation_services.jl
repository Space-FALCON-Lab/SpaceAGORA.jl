"""
    SpaceAGORA.ExternalPropagationServices

Versioned services for optional external spacecraft integrators. Call
`require_version(v"1.0.0")` before defining adapter hook methods, and again at
runtime initialization. `IMPLEMENTED_FUNCTIONS` lists the extension points;
`CONSUMED_SERVICES` lists the types and core helpers adapters may use. Existing
hook/type aliases retain their implementation owners and method identity.

Owners coexist per configuration; this facade has no global provider registry
and loads no optional engine. Each run prepares its own mutable owner runtime.
The documented step, frame, ownership and validation contracts apply equally to
all adapters. A change to signatures, result conventions or lifecycle semantics
requires a reviewed contract version. Internal module paths are unsupported.
"""
module ExternalPropagationServices

"""Exact version of this external-propagation adapter contract."""
const CONTRACT_VERSION = v"1.0.0"

"""Error raised when the requested external-propagation service version is unsupported."""
struct CompatibilityError <: Exception
    message::String
end
Base.showerror(io::IO, error::CompatibilityError) = print(io, error.message)

"""
    require_version(expected::VersionNumber)

Require the exact service version before defining adapter methods and during
adapter initialization. Throw `CompatibilityError` with requested/provided
versions on mismatch. No provider or simulation state is changed.
"""
function require_version(expected::VersionNumber)
    expected == CONTRACT_VERSION || throw(CompatibilityError(
        "External propagation adapter requires SpaceAGORA service contract $expected; this build provides $CONTRACT_VERSION."))
    return nothing
end

# Supported adapter boundary. The hook identities and engine calculations remain
# owned by their existing modules; no optional engine is loaded here.
using ..SimulationModel.ExternalPropagation: AbstractExternalPropagator
using ..SimulationModel.ExternalPropagation: external_spacecraft, external_step, external_preflight, external_prepare
using ..SimulationModel.ExternalPropagation: external_sync!, external_time, external_state, external_acceleration
using ..SimulationModel: gravity_backbone_relative_acceleration_ii
using ..SimulationEngine: sample_planet_frame_with_lpi
import ..SimulationModel
import ..SimulationEngine
import StaticArrays

export AbstractExternalPropagator, external_spacecraft, external_step, external_preflight, external_prepare
export external_sync!, external_time, external_state, external_acceleration
export with_external_propagators, validate_external_propagators, external_runtimes, external_initial_state
export sample_planet_frame_with_lpi, gravity_backbone_relative_acceleration_ii

@doc (@doc SimulationModel.ExternalPropagation.AbstractExternalPropagator) AbstractExternalPropagator
@doc (@doc SimulationModel.ExternalPropagation.external_spacecraft) external_spacecraft
@doc (@doc SimulationModel.ExternalPropagation.external_step) external_step
@doc (@doc SimulationModel.ExternalPropagation.external_preflight) external_preflight
@doc (@doc SimulationModel.ExternalPropagation.external_prepare) external_prepare
@doc (@doc SimulationModel.ExternalPropagation.external_sync!) external_sync!
@doc (@doc SimulationModel.ExternalPropagation.external_time) external_time
@doc (@doc SimulationModel.ExternalPropagation.external_state) external_state
@doc (@doc SimulationModel.ExternalPropagation.external_acceleration) external_acceleration
@doc (@doc SimulationModel.gravity_backbone_relative_acceleration_ii) gravity_backbone_relative_acceleration_ii
@doc (@doc SimulationEngine.sample_planet_frame_with_lpi) sample_planet_frame_with_lpi

"""
    with_external_propagators(config::SimulationConfiguration, owners::Tuple)

Return a configuration with exactly `owners` as its external propagators; `()`
removes all owners. Every other field and reference is preserved. This is a
shallow update: owner templates may be shared across runs, while
[`external_prepare`](@ref) must create a separate mutable runtime for each run.
The input configuration is unchanged. Ownership and solver compatibility are
checked by [`validate_external_propagators`](@ref) and again during run setup.
"""
function with_external_propagators(config::SimulationModel.SimulationConfiguration, owners::Tuple)
    return SimulationModel.SimConfig._with_configuration(config; external_propagators=owners)
end

"""
    validate_external_propagators(config::SimulationConfiguration, solver_mode::Symbol)

Check external ownership against the requested solver mode; return `nothing`
or throw `ArgumentError`. This calls each owner's [`external_preflight`](@ref),
but does not prepare a runtime or advance a simulation. Run setup still performs
these checks using its resolved solver mode.

Supported modes are `:tsit5`, `:auto_stiff`, `:rodas5p` and `:dp8`. Checkpointing,
resume, continuous events on owned spacecraft, engine articulation or control
of owned spacecraft, and callback rates incompatible with an owner's fixed
step are refused. An empty owner tuple requires no external-propagation checks.
Ensemble/process-pool routing and the flat-RHS environment guard are checked
elsewhere during run setup. Owner-defined preflight may perform its own checks;
this function is not a guarantee of complete run readiness or side-effect freedom.
"""
function validate_external_propagators(config::SimulationModel.SimulationConfiguration, solver_mode::Symbol)
    return SimulationEngine._validate_external_propagators!(config, solver_mode)
end

"""
    external_runtimes(integrator)

Return an iterable of the current run's external runtimes, in configuration
owner order, for use inside a `SaveField` getter or simulation callback. It is
empty when the run has no external owners. The iterable and runtime objects are
borrowed: callers must not change their contents, step them, or retain them for
another run. Their concrete storage type is not part of this interface.

After an accepted solver step, external synchronization precedes other discrete
callbacks. An owner may lag the solver time by less than one fixed owner step;
use [`external_state`](@ref) to obtain its spacecraft state at the callback time.
"""
external_runtimes(integrator) = integrator.p.shared_buffers.external_runtimes

"""
    external_initial_state(u0, spacecraft_index::Integer) -> (pos, vel, q, ω)

Read one spacecraft's initial state from the `u0` passed to
[`external_prepare`](@ref). The index is 1-based in the configuration's
spacecraft list, not local to the owner. Returns immutable copies: inertial
position `pos` [m], inertial velocity `vel` [m/s], scalar-last body-to-inertial
quaternion `q`, and body angular velocity `ω` [rad/s]. If orientation simulation
is disabled, `q` is the identity `(0, 0, 0, 1)` and `ω` is zero.
"""
function external_initial_state(u0, spacecraft_index::Integer)
    sc = u0.sc[spacecraft_index]
    pos = StaticArrays.SVector{3, Float64}(sc.pos)
    vel = StaticArrays.SVector{3, Float64}(sc.vel)
    q = hasproperty(sc, :q) ? StaticArrays.SVector{4, Float64}(sc.q) : StaticArrays.SVector(0.0, 0.0, 0.0, 1.0)
    ω = hasproperty(sc, :ω) ? StaticArrays.SVector{3, Float64}(sc.ω) : StaticArrays.SVector(0.0, 0.0, 0.0)
    return (; pos, vel, q, ω)
end

"""Hook identities to which an adapter may add methods for its own types."""
const IMPLEMENTED_FUNCTIONS = (
    :external_spacecraft, :external_step, :external_preflight, :external_prepare,
    :external_sync!, :external_time, :external_state, :external_acceleration,
)

"""Core types and read/configuration services supplied to external adapters."""
const CONSUMED_SERVICES = (
    :AbstractExternalPropagator, :with_external_propagators, :validate_external_propagators,
    :external_runtimes, :external_initial_state, :sample_planet_frame_with_lpi,
    :gravity_backbone_relative_acceleration_ii,
)

public CONTRACT_VERSION, CompatibilityError, require_version, IMPLEMENTED_FUNCTIONS, CONSUMED_SERVICES
end # module ExternalPropagationServices
