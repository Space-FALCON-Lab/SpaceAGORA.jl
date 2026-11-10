# The spacecraft thermal-node record, mounted on the root link with
# `ThermalNodes.add_thermal!`. Included by `Components`
# (src/models/vehicle/spacecraft/components.jl), which owns and exports the
# record, as it does the power and data records.

"""
    ThermalNode(; heat_capacity_j_k, conductance_w_k, initial_k, minimum_k, maximum_k,
                optical_absorptivity=0.0, optical_reflectivity=0.0)

The spacecraft as one lumped thermal node, mounted on the root link (at most
one). With internal heat `Q` (W) and an environment sink temperature `T_env` (K),
its temperature `T` (K) follows `C dT/dt = Q - H (T - T_env)`, where `C` is
`heat_capacity_j_k` and `H` the effective linear heat-rejection
`conductance_w_k`. `minimum_k` and `maximum_k` are its allowed range and
`initial_k` its temperature at the start of a run.

`optical_absorptivity` and `optical_reflectivity` are the fractions of an incident
laser beam the spacecraft absorbs as heat and reflects; their sum is at most 1.

This is a single-node model, not a radiative or conductive network. The engine
does not integrate it; the constellation resource scheduler
(`plan_resource_schedule`) reads it.
"""
struct ThermalNode
    heat_capacity_j_k::Float64
    conductance_w_k::Float64
    initial_k::Float64
    minimum_k::Float64
    maximum_k::Float64
    optical_absorptivity::Float64
    optical_reflectivity::Float64

    function ThermalNode(heat_capacity_j_k::Real, conductance_w_k::Real, initial_k::Real, minimum_k::Real,
                         maximum_k::Real, optical_absorptivity::Real, optical_reflectivity::Real)
        (isfinite(heat_capacity_j_k) && heat_capacity_j_k > 0.0) || throw(ArgumentError(
            "ThermalNode heat_capacity_j_k must be finite and positive, got $(heat_capacity_j_k)."))
        (isfinite(conductance_w_k) && conductance_w_k > 0.0) || throw(ArgumentError(
            "ThermalNode conductance_w_k must be finite and positive, got $(conductance_w_k)."))
        (isfinite(minimum_k) && isfinite(maximum_k) && 0.0 < minimum_k <= maximum_k) || throw(ArgumentError(
            "ThermalNode needs 0 < minimum_k <= maximum_k, got $(minimum_k) and $(maximum_k)."))
        (minimum_k <= initial_k <= maximum_k) || throw(ArgumentError(
            "ThermalNode initial_k must be in [minimum_k, maximum_k], got $(initial_k)."))
        (0.0 <= optical_absorptivity <= 1.0 && 0.0 <= optical_reflectivity <= 1.0 &&
            optical_absorptivity + optical_reflectivity <= 1.0 + 1e-12) || throw(ArgumentError(
            "ThermalNode optical_absorptivity and optical_reflectivity must lie in [0, 1] with a sum of at most 1."))
        return new(Float64(heat_capacity_j_k), Float64(conductance_w_k), Float64(initial_k), Float64(minimum_k),
                   Float64(maximum_k), Float64(optical_absorptivity), Float64(optical_reflectivity))
    end
end

ThermalNode(; heat_capacity_j_k::Real, conductance_w_k::Real, initial_k::Real, minimum_k::Real, maximum_k::Real,
            optical_absorptivity::Real=0.0, optical_reflectivity::Real=0.0) =
    ThermalNode(heat_capacity_j_k, conductance_w_k, initial_k, minimum_k, maximum_k,
                optical_absorptivity, optical_reflectivity)
