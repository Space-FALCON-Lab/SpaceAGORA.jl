Base.@kwdef struct ParallelConfig
    profile::String = ""
    outer_parallel_active::Bool = false
    parallel_policy_adaptive::Bool = false
    effector_parallel_mode::String = "auto"
    rhs_batch_parallel_mode::String = "auto"
    density_callback_parallel_mode::String = "auto"
    control_callback_parallel_mode::String = "auto"
    thermal_callback_parallel_mode::String = "auto"
end
