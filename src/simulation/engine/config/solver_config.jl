"""
    SolverConfig

Typed runtime configuration for solver selection, split-IMEX settings, and
multirate integration policy.

`split_imex` uses the atmosphere-implicit IMEX partition. `multirate` keeps
the control-focused split path. `gravity_backbone_split` is a translational
fixed-step gravity-backbone split mode built on a symplectic gravity core with
explicit SRP / N-body velocity kicks; it is not a fully symplectic whole-system
solve.
"""
Base.@kwdef struct SolverConfig
    mode::String = ""
    maxiters::Union{Nothing, Int} = nothing
    split_imex_solver::String = "kencarp4"
    gravity_backbone_dt_s::Union{Nothing, Float64} = nothing
    multirate_fast_substeps::Int = 8
    multirate_slow_dt_s::Union{Nothing, Float64} = nothing
    multirate_slow_solver::String = "tsit5"
    multirate_fast_solver::String = "auto_stiff"
end
