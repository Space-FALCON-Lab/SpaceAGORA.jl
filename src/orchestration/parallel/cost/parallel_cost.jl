module ParallelCost

using TOML
using SHA
import ...RuntimeServices

using ..ParallelPolicy

export MachineConstants
export ContentionInputs, alloc_bytes_per_second, lock_duty_cycle, lock_wait_hold_ratio
export calibrate_machine, machine_fingerprint, machine_constants_path
export save_machine_constants, load_machine_constants
export ensure_machine_constants!, machine_constants_generation
export timed_min, paired_compare, PairedComparison
export StreamingPairedTrial, next_arm, observe!, trial_active, trial_verdict
export trial_rounds, trial_speedups
export calib_batch_kernel!, calib_table, calib_alloc_kernel!
export usl_speedup, usl_peak_workers
export timing_noise_floor

include(joinpath(@__DIR__, "cost_terms.jl"))
include(joinpath(@__DIR__, "robust_timing.jl"))
include(joinpath(@__DIR__, "calibration_kernels.jl"))
include(joinpath(@__DIR__, "machine_calibration.jl"))

end # module ParallelCost
