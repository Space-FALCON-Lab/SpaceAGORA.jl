using OrdinaryDiffEq
using StaticArrays
using .SimulationModel

function reaction_wheel_model!(
    link::Link,
    τ::MVector{3, Float64},
    dt::Float64
)
    """
    Update the angular velocity of the reaction wheels of a given link based on a desired torque.
    # Arguments
    - `link`: Link object containing the reaction wheel properties.
    - `τ`: Desired torque vector applied to link, expressed in body frame, SVector{3, Float64}.
    - `dt`: Time step for the update, Float64.
    # Returns
    - `ω_new`: Updated angular velocity vector, SVector{3, Float64}.
    """
    function ω_dot!(du, u, p, t)
        du[:] .= p[1]*p[2]  # du = J_rw \ τ
    end

    # Integrate the constant wheel acceleration over one control interval.
    prob = ODEProblem(ω_dot!, link.rw_assembly.h_wheels, (0.0, dt), [pinv(link.rw_assembly.J_rw), τ])
    link.rw_assembly.h_wheels .= solve(prob, Tsit5()).u[end]
end
