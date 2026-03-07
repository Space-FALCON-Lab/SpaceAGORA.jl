"""
    SolverConfig

Typed runtime configuration for solver selection, split-IMEX settings, and
multirate integration policy.
"""
Base.@kwdef struct SolverConfig
    mode::String = ""
    maxiters::Union{Nothing, Int} = nothing
    split_imex_solver::String = "kencarp4"
    multirate_fast_substeps::Int = 8
    multirate_slow_dt_s::Union{Nothing, Float64} = nothing
    multirate_slow_solver::String = "tsit5"
    multirate_fast_solver::String = "auto_stiff"
end
