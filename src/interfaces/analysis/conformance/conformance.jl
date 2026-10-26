"""
    SpaceAGORA.Conformance

Checks that code built around SpaceAGORA (a force/torque model, a service facade)
behaves as the engine's contracts require. Each `check_*` function runs the engine's
own machinery (its pre-run hook check and `run_simulation`) on the object and returns
a [`ConformanceReport`](@ref); it needs no test framework. With the `Test` standard
library loaded, each `test_*` function records the same checks as tests:

```julia
using Test, SpaceAGORA
using SpaceAGORA.Conformance

@testset "MyModels" begin
    test_force_torque_model(MyDisturbanceTorque(gain=1e-6))
end
```

A conformance pass certifies interface behaviour only. Scientific validation of a
model, its parameters and the claims made from it stay with its author.
"""
module Conformance

import Logging
import Random
import SciMLBase
import TOML
import ..SimulationModel
import ..SimulationEngine
import ..SimulationLifecycle
import ..GNCRuntime
import ..GNCAlgorithms
import ..RPOPlannerInterfaces
import ..RPOPlannerLifecycle
import ..AlgorithmLifecycle
import ..BurnPlanning
import ..JointReferencePlanning

export ConformanceCheck, ConformanceReport, ConformanceFixture, passed, conformance_fixture
export check_force_torque_model, check_service_facade, check_sensor
export check_guidance_model, check_navigation_model, check_control_model
export check_algorithm, check_algorithm_family, rpo_conformance_requests, rpo_trusted_queries, rpo_conformance_simulation
export burn_plan_conformance_requests, burn_plan_conformance_simulation
export joint_reference_conformance_requests, joint_reference_conformance_simulation
export test_report, test_force_torque_model, test_service_facade, test_sensor
export test_guidance_model, test_navigation_model, test_control_model
export test_algorithm, test_algorithm_family

include(joinpath(@__DIR__, "reports.jl"))
include(joinpath(@__DIR__, "fixtures.jl"))
include(joinpath(@__DIR__, "walker.jl"))
include(joinpath(@__DIR__, "checkpoints.jl"))
include(joinpath(@__DIR__, "effectors.jl"))
include(joinpath(@__DIR__, "facades.jl"))
include(joinpath(@__DIR__, "sensors.jl"))
include(joinpath(@__DIR__, "models.jl"))
include(joinpath(@__DIR__, "algorithms.jl"))
include(joinpath(@__DIR__, "rpo.jl"))
include(joinpath(@__DIR__, "burn_plan.jl"))
include(joinpath(@__DIR__, "joint_reference.jl"))

# The check ids of each report kind, in report order. Core's negative controls
# fail when an id here has no control.
const CHECK_IDS = (
    force_torque_model = (:hooks_present, :supported_sample_signature, :repeatable_run,
        :deepcopy_equivalent, :global_rng_independent, :concurrent_equals_sequential, :output_hook_pure,
        :environment_declared),
    service_facade = (:facade_shape,),
    sensor = (:run_state_fields, :requirements_supported, :measurement_shape, :repeatable_run, :deepcopy_equivalent,
        :global_rng_independent, :concurrent_equals_sequential, :seed_repeatable, :custom_state_reset,
        :checkpoint_round_trip),
    guidance_model = (:hooks_present, :repeatable_run, :deepcopy_equivalent, :global_rng_independent,
        :concurrent_equals_sequential, :initialize_resets, :bind_fresh_state, :checkpoint_round_trip),
    navigation_model = (:hooks_present, :repeatable_run, :deepcopy_equivalent, :global_rng_independent,
        :concurrent_equals_sequential, :initialize_resets, :bind_fresh_state, :checkpoint_round_trip),
    control_model = (:hooks_present, :supported_sample_signature, :repeatable_run,
        :deepcopy_equivalent, :global_rng_independent, :concurrent_equals_sequential, :initialize_resets,
        :output_hook_pure, :positive_consumption, :bind_fresh_state, :checkpoint_round_trip, :outputs_in_scope),
    algorithm = _ALGORITHM_CHECK_IDS,
    algorithm_family = _ALGORITHM_FAMILY_CHECK_IDS,
)

# Defined by ext/SpaceAGORATestExt.jl when Test is loaded. Today Test is always
# loaded with SpaceAGORA (MathOptInterface, a dependency of OSQP, loads it), so the
# hint below matters only if that changes.
const _TEST_FUNCTIONS = (test_report, test_force_torque_model, test_service_facade, test_sensor, test_guidance_model,
    test_navigation_model, test_control_model, test_algorithm, test_algorithm_family)

function _test_extension_hint(io::IO, exc::MethodError, extension_loaded::Bool)
    extension_loaded && return nothing
    any(f -> exc.f === f, _TEST_FUNCTIONS) || return nothing
    print(io, "\nThe test_* conformance functions come from SpaceAGORA's Test extension: ",
        "load the Test standard library first (`using Test`). Without it, the matching ",
        "check_* function returns the same checks as a ConformanceReport.")
    return nothing
end

function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        loaded = Base.get_extension(parentmodule(@__MODULE__), :SpaceAGORATestExt) !== nothing
        _test_extension_hint(io, exc, loaded)
    end
    return nothing
end

end # module Conformance
