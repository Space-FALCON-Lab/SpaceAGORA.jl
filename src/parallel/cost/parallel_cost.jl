module ParallelCost

using TOML
using SHA

using ..EnvironmentModels
using ..ParallelPolicy
using Polyester: @batch
import Polyester
using ..DynamicEffectors: ConstantGravityModel, InverseSquaredGravityModel
using ..DynamicEffectors: InverseSquaredJ2GravityModel
using ..DynamicEffectors: NBodyGravityModel, GravitationalHarmonicsModel
using ..DynamicEffectors: SolarRadiationPressureModel

export WorkCounts, MachineConstants, RateCurve, rate_at, effector_cost_terms
export constellation_work_counts, model_in_cost_domain, flat_queue_node_effector
export EffectorProbe, validate_declaration, probe_ns_map
export PlanCandidate, PlanPrediction, predict_plan_ns, plan_candidates, select_plan
export calibrate_machine, machine_fingerprint, machine_constants_path
export save_machine_constants, load_machine_constants, constants_are_current
export timed_min, paired_compare, PairedComparison
export calib_batch_kernel!, calib_touch_kernel!, calib_table, calib_scalar_kernel!
export reference_kernel_ns, reference_memory_kernel_ns
export timing_noise_floor

include(joinpath(@__DIR__, "cost_terms.jl"))
include(joinpath(@__DIR__, "robust_timing.jl"))
include(joinpath(@__DIR__, "calibration_kernels.jl"))
include(joinpath(@__DIR__, "work_counts.jl"))
include(joinpath(@__DIR__, "effector_probe.jl"))
include(joinpath(@__DIR__, "machine_calibration.jl"))
include(joinpath(@__DIR__, "cost_predictor.jl"))

end # module ParallelCost
