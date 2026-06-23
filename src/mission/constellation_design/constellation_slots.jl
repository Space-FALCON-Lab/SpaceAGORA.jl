module ConstellationSlots

using Base.Threads
using ProgressMeter
using DataFrames, CSV
using LinearAlgebra
using Statistics
using Random

const DEFAULT_MAX_SHELL_ROWS = 2_000_000

function _angular_span_rad(min_deg, max_deg)
    min_rad = deg2rad(Float64(min_deg))
    max_rad = deg2rad(Float64(max_deg))
    return max_rad - min_rad
end

function _validate_shell_row(row::DataFrameRow, required_cols::Vector{Symbol})
    for col in required_cols
        if !haskey(row, col)
            error("Missing required column: $col")
        end
    end
    return true
end

function _filter_shell_by_bounds(df::DataFrame, 
                                  a_min::Real, a_max::Real,
                                  e_min::Real, e_max::Real,
                                  i_min::Real, i_max::Real,
                                  raan_min::Real, raan_max::Real,
                                  arg_p_min::Real, arg_p_max::Real,
                                  ta_min::Real, ta_max::Real)
    
    filtered = filter(row -> 
        a_min <= row.a <= a_max &&
        e_min <= row.e <= e_max &&
        i_min <= row.inc <= i_max &&
        raan_min <= row.raan <= raan_max &&
        arg_p_min <= row.arg_p <= arg_p_max &&
        ta_min <= row.ta <= ta_max,
        df)
    
    return filtered
end

function _prune_shell(df::DataFrame, keep_indices::Vector{Int})
    if isempty(keep_indices)
        return DataFrame()
    end
    return df[keep_indices, :]
end

function _sample_shell(df::DataFrame, n_samples::Integer; rng::AbstractRNG=Random.default_rng())
    n_rows = nrow(df)
    if n_samples >= n_rows
        return df
    end
    sample_indices = rand(rng, 1:n_rows, n_samples)
    return df[sample_indices, :]
end

export _angular_span_rad, _validate_shell_row, _filter_shell_by_bounds
export _prune_shell, _sample_shell

end # module ConstellationSlots
