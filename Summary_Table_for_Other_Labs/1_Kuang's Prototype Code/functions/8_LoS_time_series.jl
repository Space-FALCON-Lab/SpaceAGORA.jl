"""
This module computes line-of-sight metrics over time for all satellite pairs.
"""

"""
    Compute line-of-sight metrics over time for all satellite pairs.

    Inputs:
        sol: solution object from ODE solver (contains time series of states)
        R_atm: atmosphere radius for line-of-sight blocking (default R_ATMDEF)

    Returns:
        t: time vector
        metrics: dictionary of metrics keyed by (i,j) tuples, where metrics[(i,j)] = Dict(:clearance=>..., :range=>...)
"""
function los_time_series(sol; R_atm=R_ATMDEF)

    N = Int(length(sol.u[1]) ÷ 6) # number of satellites
    t = sol.t
    metrics = Dict{Tuple{Int,Int}, Dict{Symbol, Vector{Float64}}}()
    # metrics is a dictionary where the keys are tuples of satellite indices (i,j)
    # and the values are dictionaries with keys :clearance and :range, each mapping to a vector of Float64 values
    for i in 1:N-1, j in i+1:N # for each unique pair of satellites
        clearance = zeros(length(t))
        range     = zeros(length(t))
        for k in eachindex(t) # for each time step
            uk = sol.u[k] # state vector at time step k
            ri = @SVector [uk[idx(i,1)], uk[idx(i,2)], uk[idx(i,3)]] # position of satellite i
            rj = @SVector [uk[idx(j,1)], uk[idx(j,2)], uk[idx(j,3)]] # position of satellite j
            met = los_metrics(ri, rj; R_atm=R_atm) # compute metrics
            clearance[k] = met.clearance # atmosphere clearance
            range[k]     = met.slant_range
        end
        metrics[(i,j)] = Dict(:clearance=>clearance, :range=>range) # store results in dictionary
    end
    return t, metrics
end

"""
    Print a summary of line-of-sight metrics over time for all satellite pairs.

    Inputs:
        sol: solution object from ODE solver (contains time series of states)
        R_atm: atmosphere radius for line-of-sight blocking (default R_ATMDEF)
"""
function print_los_summary(sol; R_atm=R_ATMDEF)

    t, metrics = los_time_series(sol; R_atm=R_atm)
    for (pair, data) in metrics
        clear = data[:clearance]; dist = data[:range] # clear is atmosphere clearance, dist is slant range
        minclear, kc = findmin(clear)
        minrange, kr = findmin(dist)
        @printf("\nLoS pair %d-%d:\n", pair[1], pair[2])
        @printf("  min atmosphere clearance = %.1f m at t=%.1f s\n", minclear, t[kc])
        @printf("  min slant range          = %.1f km at t=%.1f s\n", minrange/1e3, t[kr])
    end
end