# Test RW model
# include("../config.jl")
# using Plots
# import .config
using OrdinaryDiffEq
# using SimpleDiffEq
using StaticArrays
# using BenchmarkTools

function reaction_wheel_model!(
    link,
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
        # u: current angular velocity vector
        # p: parameters (link properties)
        # t: time (not used here)
        du[:] .= p[1]*p[2]  # du = J_rw \ τ

        # println("du: $du")
    end

    prob = ODEProblem(ω_dot!, link.rw, (0.0, dt), [pinv(link.J_rw), τ])
    # println("link.rw: $(link.rw)")
    link.rw .= solve(prob, BS3()).u[end]
    # println("Updated link.rw: $(link.rw)")
    # println("link.rw post: $(link.rw)")
end



# spacecraft = config.SpacecraftModel()
# # Add bodies to the spacecraft model
# main_bus = config.Link(root=true, 
#                         r=SVector{3, Float64}(0.0, 0.0, 0.0), 
#                         q=SVector{4, Float64}([0, 0, 0, 1]),
#                         ṙ=SVector{3, Float64}([0,0,0]), 
#                         dims=SVector{3, Float64}([2.2,2.6,1.7]), 
#                         ref_area=2.6*1.7,
#                         m=391.0, 
#                         gyro=3,
#                         J_rw=MMatrix{3, 3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])) # Reaction wheel inertia

# L_panel = config.Link(r=SVector{3, Float64}(0.0, -2.6/2 - 3.89/4, 0.0), 
#                         q=SVector{4, Float64}([0, 0, 0, 1]),
#                         ṙ=SVector{3, Float64}([0,0,0]), 
#                         dims=SVector{3, Float64}([0.01, 3.89/2, 1.7]), 
#                         ref_area=3.89*1.7/2,
#                         m=10.0, 
#                         gyro=0)
# R_panel = config.Link(r=SVector{3, Float64}(0.0, 2.6/2 + 3.89/4, 0.0),
#                         q=SVector{4, Float64}([0, 0, 0, 1]),
#                         ṙ=SVector{3, Float64}([0,0,0]), 
#                         dims=SVector{3, Float64}([0.01, 3.89/2, 1.7]), 
#                         ref_area=3.89*1.7/2,
#                         m=10.0, 
#                         gyro=0)

# config.add_body!(spacecraft, main_bus, prop_mass=50.0)
# config.add_body!(spacecraft, L_panel)
# config.add_body!(spacecraft, R_panel)

# L_panel_joint = config.Joint(main_bus, L_panel)
# R_panel_joint = config.Joint(R_panel, main_bus)
# config.add_joint!(spacecraft, L_panel_joint)
# config.add_joint!(spacecraft, R_panel_joint)

# println("Spacecraft model initialized with $(length(spacecraft.links)) bodies.")
# # println("Spacecraft roots: $spacecraft.roots")
# println("Spacecraft COM: $(config.get_COM(spacecraft, main_bus))")
# println("Spacecraft MOI: $(config.get_inertia_tensor(spacecraft, main_bus))")

# function test_reaction_wheel_model()
#     """
#     Test the reaction wheel model by applying a torque and updating the angular velocity.
#     """
#     dt = 0.1  # Time step
#     τ = MVector{3, Float64}([0.1; 0.2; 0.3])  # Desired torque vector

#     # Update the reaction wheel model
#     reaction_wheel_model!(main_bus, τ, dt)

#     # Print the updated angular velocity of the main bus
#     # println("Updated angular momentum of main bus: $(main_bus.rw)")
# end

# @benchmark test_reaction_wheel_model()


