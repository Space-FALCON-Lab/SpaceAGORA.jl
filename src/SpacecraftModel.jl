# include("utils/quaternion_utils.jl")
# import .quaternion_utils
# using StaticArrays
# using LinearAlgebra
# using CSV
# using DataFrames

# const I3 = SMatrix{3, 3, Float64}(diagm(ones(3)))

# @kwdef mutable struct Facet
#     area::Float64 = 0.0 # Area of the facet
#     attitude::MVector{4, Float64} = @MVector [0.0, 0.0, 0.0, 1.0] # Attitude of the facet relative to the body frame
#     normal_vector::MVector{3, Float64} = @MVector [1.0, 0.0, 0.0] # Normal vector in the facet body frame, for flat plate this is in the x-direction
#     cp::MVector{3, Float64} = @MVector [0.0, 0.0, 0.0] # Center of pressure with respect to the center of mass of the Link containing the facet, in the Link frame
#     ρ::Float64 = 0.0 # Diffuse coefficient
#     δ::Float64 = 0.0 # Specular coefficient
#     name::String = "" # Identifying name
# end

# @kwdef mutable struct Thruster
#     max_thrust::Float64 = 1.0 # Maximum thrust, N
#     location::MVector{3, Float64} = MVector{3, Float64}(zeros(3)) # Location in the link frame, relative to the CoM of the link, m
#     direction::MVector{3, Float64} = MVector{3, Float64}(zeros(3)) # Unit vector direction of thrust in the link frame, n/d
#     Isp::Float64 = 0.0 # Specific impulse of the thruster, s
#     cutoff_frequency::Float64 = 1.0 # Governing parameter for thruster ramp-up/ramp-down, rad/s. See https://avslab.github.io/basilisk/Documentation/simulation/dynamics/Thrusters/thrusterStateEffector/thrusterStateEffector.html for details
#     min_firing_time::Float64 = 0.0 # Minimum time for the thruster to be firing, s.
#     level_on::Float64 = 0.75 # Upper threshold for the Schmitt trigger, nd
#     level_off::Float64 = 0.25 # Lower threshold for the Schmitt trigger, nd
#     κ::Float64 = 0.0 # Thrust factor, nd. This value is integrated during ramp up/down as part of the first-order filter describing the true thrust from the thruster.
#     thrust::Float64 = 0.0 # Current thrust magnitude, to be updated during simulation. N
#     stop_firing_time::Float64 = 0.0 # Time at which the thruster should stop firing, s. This is used to determine how long to continue firing the thruster into the next timestep
# end

# @kwdef mutable struct Magnet
#     m::MVector{3, Float64} = MVector{3, Float64}(zeros(3)) # Magnetic dipole moment in the body frame, nT
#     location::MVector{3, Float64} = MVector{3, Float64}(zeros(3)) # Location in the link frame, relative to the CoM of the link, m
# end

# mutable struct Link
#     root::Bool # Whether this link is a root link (i.e., the main bus or core body of the spacecraft). 
#     #^ This will have orientation expressed relative to the inertial frame.
#     r::MVector{3, Float64} # Position of COM (Body frame for non-root, inertial frame for root)
#     q::MVector{4, Float64} # Orientation (Body frame for non-root, inertial frame for root)
#     ṙ::MVector{3, Float64} # velocity (body frame for non-root, inertial frame for root)
#     ω::MVector{3, Float64} # angular velocity (body frame for non-root, inertial frame for root)
#     dims::MVector{3, Float64} # Size[x,y,z] Box= x=thikness, y=z=width, height
#     ref_area::Float64 # Reference area for aerodynamic calculations
#     m::Float64 # Mass
#     mass::SMatrix{3, 3, Float64} # Mass Matrix
#     inertia::SMatrix{3, 3, Float64} # Inertial matrix
#     aᵇ::SVector{3, Float64} # Left extent (Body frame)
#     bᵇ::SVector{3, Float64} # Right extent (Body frame)
#     α::Float64 # Angle of attack, rad
#     β::Float64 # Sideslip angle, rad
#     θ::Float64 # Flow angle, rad
#     gyro::Int64 # Number of Gyroscope
#     max_torque::Float64 # Maximum torque that can be applied by the reaction wheels
#     max_h::Float64 # Maximum angular momentum that can be stored in the reaction wheels
#     rw::Vector{Float64} # angular momentum reaction wheels
#     J_rw::Matrix{Float64} # reaction wheel jacobian
#     rw_τ::MVector{3, Float64} # Reaction wheel torque vector, to be updated at each simulation step
#     net_force::MVector{3, Float64} # Net force acting on the link, to be updated at each simulation step
#     net_torque::MVector{3, Float64} # Net torque acting on the link, to be updated at each simulation step
#     attitude_control_function::Function # Function to control the attitude of the link (apply a torque via reaction wheels, thrusters, etc.), assumes rigid joints
#     actuation_function::Function # Function to apply actuation forces/torques at the link (change orientation of solar panels, etc.)  
#     attitude_control_rate::Float64 # Rate at which the attitude control function is called, in seconds
#     ω_wheel_derivatives::Vector{Float64} # Angular momentum derivatives of the reaction wheels, to be updated at each simulation step
#     SRP_facets::Vector{Facet}
#     J_thruster::Matrix{Float64} # Thruster Jacobian matrix
#     thrusters::Vector{Thruster}
#     thrust_calculation_function::Function # Function to calculate the average thrust over the control period, e.g., Schmitt trigger
#     magnets::Vector{Magnet} # List of magnetic dipoles attached to the link
#     function Link(;root=false,
#                     r=SVector{3, Float64}([0, 0, 0]), 
#                     q=SVector{4, Float64}([0,0,0,1]), 
#                     ṙ=SVector{3, Float64}([0,0,0]), 
#                     ω=SVector{3, Float64}([0,0,0]), 
#                     dims=SVector{3, Float64}([0.5, 0.5, 0.1]), 
#                     ref_area=1.0,
#                     m = 3.0,
#                     mass=SMatrix{3, 3, Float64}(m*I3), 
#                     inertia=SMatrix{3, 3, Float64}(1 / 12 * m * diagm([dims[2]^2 + dims[3]^2; dims[1]^2 + dims[3]^2; dims[1]^2 + dims[2]^2])),
#                     a=SVector{3, Float64}([-0.5*dims[1], 0, 0]),
#                     b=SVector{3, Float64}([0.5*dims[1], 0, 0]),
#                     α=pi/2.0,
#                     β=0.0,
#                     θ=0.0,
#                     gyro = 0,
#                     max_torque = 0.25,
#                     max_h = 70.0,
#                     rw = MVector{gyro, Float64}(zeros(gyro)),
#                     J_rw=MMatrix{3, gyro, Float64}(zeros(3, gyro)),
#                     rw_τ=MVector{3, Float64}(zeros(3)),
#                     net_force=MVector{3, Float64}(zeros(3)),
#                     net_torque=MVector{3, Float64}(zeros(3)),
#                     attitude_control_function=()->0.0,
#                     actuation_function=()->0.0,
#                     attitude_control_rate=0.1,
#                     ω_wheel_derivatives=MVector{gyro, Float64}(zeros(gyro)),
#                     SRP_facets=Facet[],
#                     J_thruster=Matrix{Float64}(zeros(3, 1)),
#                     thrusters=Thruster[],
#                     thrust_calculation_function=()->0.0,
#                     magnets=Magnet[])#SMatrix{3,Int(gyro)}(1.0I))
#                     println(length(rw))
#         new(root, r, q, ṙ, ω, dims, ref_area, m, mass, inertia, a, b, α, β, θ, gyro, max_torque, max_h, rw, J_rw, rw_τ, net_force, net_torque, attitude_control_function, actuation_function, attitude_control_rate, ω_wheel_derivatives, SRP_facets, J_thruster, thrusters, thrust_calculation_function, magnets)
#     end

#     function Link(link::Link)
#         new(link.root, link.r, link.q, link.ṙ, link.ω, link.dims, link.ref_area, link.m, link.mass, link.inertia, link.aᵇ, link.bᵇ, link.gyro, link.max_torque, link.max_h, copy(link.rw), copy(link.J_rw), copy(link.rw_τ), copy(link.net_force), copy(link.net_torque), 
#             link.attitude_control_function, link.actuation_function, link.attitude_control_rate, copy(link.ω_wheel_derivatives), copy(link.SRP_facets), copy(link.J_thruster), copy(link.thrusters), link.thrust_calculation_function)
#     end
# end

# mutable struct Joint
#     link1::Link
#     link2::Link
#     p1ᵇ::SVector{3, Float64} # Position of joint in body frame of link1
#     p2ᵇ::SVector{3, Float64} # Position of joint in body frame of link2
#     Kx::SMatrix{3, 3} # Stiffness Matrix
#     Kt::SMatrix{3, 3} # Stiffness Matrix
#     Cx::SMatrix{3, 3} # Damping Matrix
#     Ct::SMatrix{3, 3} # Damping Matrix
#     translational_displacement::SVector{3, Float64} # Translational displacement vector
#     rotational_displacement::SVector{4, Float64} # Rotational displacement vector (quaternion)
#     function Joint(link1::Link, p1ᵇ::SVector{3, Float64}, 
#                     link2::Link, p2ᵇ::SVector{3, Float64}, 
#                     Kx=SMatrix{3,3, Float64}(0.0I), 
#                     Kt=SMatrix{3,3, Float64}(0.0I), 
#                     Cx=zeros(SMatrix{3,3, Float64}),
#                     Ct=zeros(SMatrix{3,3, Float64}))
#         new(link1, link2, p1ᵇ, p2ᵇ, Kx, Kt, Cx, Ct,
#             SVector{3, Float64}(0.0, 0.0, 0.0), SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))
#     end

#     function Joint(link1::Link, link2::Link; p1=link1.bᵇ, 
#                                  p2=link2.aᵇ, 
#                                  Kx=SMatrix{3,3, Float64}(0.0I), 
#                                  Kt=SMatrix{3,3, Float64}(0.0I), 
#                                  Cx=zeros(SMatrix{3,3, Float64}),
#                                  Ct=zeros(SMatrix{3,3, Float64}),
#                                  translational_displacement=SVector{3, Float64}(0.0, 0.0, 0.0),
#                                  rotational_displacement=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))
        
#         new(link1, link2, p1, p2, Kx, Kt, Cx, Ct, 
#             translational_displacement, rotational_displacement)
#     end

#     function Joint(;link1=Link(), link2=Link(), p1=link1.bᵇ, 
#         p2=link2.aᵇ, 
#         Kx=SMatrix{3,3, Float64}(1.0I), 
#         Kt=SMatrix{3,3, Float64}(1.0I), 
#         Cx=zeros(SMatrix{3,3, Float64}),
#         Ct=zeros(SMatrix{3,3, Float64}),
#         translational_displacement=SVector{3, Float64}(0.0, 0.0, 0.0),
#         rotational_displacement=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))

#         new(link1, link2, p1, p2, Kx, Kt, Cx, Ct, 
#             translational_displacement, rotational_displacement)
#     end

#     function Joint(joint::Joint)
#         new(joint.link1, joint.link2, joint.p1ᵇ, joint.p2ᵇ, 
#             joint.Kx, joint.Kt, joint.Cx, joint.Ct,
#             joint.translational_displacement, joint.rotational_displacement)
#     end
# end

# mutable struct SpacecraftModel
#     joints::Vector{Joint} # List of joints
#     links::Vector{Link} # List of links (bodies)
#     roots::Vector{Link} # Vector of root links (main bus or core bodies)
#     instant_actuation::Bool # Whether control inputs (e.g., solar panel angles) are applied instantly
#     # dry_mass::Float64 # Dry mass of the spacecraft
#     prop_mass::Vector{Float64} # Fuel mass available for maneuvers
#     inertia_tensors::Vector{SMatrix{3, 3, Float64}} # Inertia tensors of the spacecraft bodies in the body frame
#     n_reaction_wheels::Int64 # Number of reaction wheels in the spacecraft model
#     n_thrusters::Int64 # Number of thrusters in the spacecraft model
#     dynamic_effectors::Tuple{AbstractForceTorqueModel} # List of dynamic effector models (gravity, drag, etc.)
#     # inertia_tensor::MMatrix{3, 3, Float64} # Inertia tensor of the spacecraft in the inertial frame
#     # COM::SVector{3, Float64} # Center of mass in the body frame (relative to origin of the root body)

#     function SpacecraftModel(;joints=Joint[], links=Link[], roots=Link[], 
#                         instant_actuation=true, 
#                         prop_mass=Float64[],
#                         inertia_tensors=SMatrix{3, 3, Float64}[],
#                         n_reaction_wheels=0,
#                         n_thrusters=0)
#         new(joints, links, roots, instant_actuation, prop_mass, inertia_tensors, n_reaction_wheels, n_thrusters)
#     end
# end

# # Function to add a body to the spacecraft model
# function add_body!(model::SpacecraftModel, 
#                 body::Link; prop_mass::Union{Nothing, Float64}=nothing)
#     """
#     Adds a new body to the spacecraft model.
#     - `model`: The spacecraft model to which the body will be added.
#     - `body`: The rigid body to be added.
#     """

#     push!(model.links, body) # Add the body to the links vector
#     if body.root
#         # If the body is a root body, it should be added to the roots vector
#         @assert prop_mass !== nothing "Propellant mass must be provided for root body"
#         push!(model.roots, body) # Add the body to the roots vector
#         push!(model.prop_mass, prop_mass) # Add propellant mass for root body
#         # Initialize inertia tensor for the root body
#     end
#     model.n_reaction_wheels += body.gyro # Increment the number of reaction wheels
#     if body.root
#         update_inertia_tensor!(model, body) # Calculate inertia tensor for the body
#     end
# end

# # Function to add a joint to the spacecraft model
# function add_joint!(model::SpacecraftModel, joint::Joint)
#     """
#     Adds a new joint to the spacecraft model.
#     - `model`: The spacecraft model to which the joint will be added.
#     - `joint`: The joint to be added.
#     """
#     push!(model.joints, joint) # Add the joint to the joints vector
#     update_inertia_tensor!(model, joint.link1) # Calculate inertia tensor
# end

# function add_facet!(link::Link, facet::Facet)
#     """
#     Adds a facet to the Link
#     - `link` : The Link to which the facet is added
#     - `facet` : The facet to add
#     """
#     push!(link.SRP_facets, facet)
# end

# function add_facet!(link::Link, facets::Vector{Facet})
#     """
#     Adds a facet to the Link
#     - `link` : The Link to which the facet is added
#     - `facet` : The facet to add
#     """
#     push!(link.SRP_facets, facets...)
# end

# function add_magnet!(link::Link, magnet::Magnet)
#     """
#     Adds a magnetic dipole to the Link
#     - `link` : The Link to which the magnet is added
#     - `magnet` : The magnet to add
#     """
#     push!(link.magnets, magnet)
# end

# function add_magnet!(link::Link, magnets::Vector{Magnet})
#     """
#     Adds a magnetic dipole to the Link
#     - `link` : The Link to which the magnet is added
#     - `magnet` : The magnet to add
#     """
#     push!(link.magnets, magnets...)
# end

# function add_dynamic_effector!(model::SpacecraftModel, effector::Any)
#     """
#     Adds a dynamic effector model to the spacecraft model.
#     - `model`: The spacecraft model to which the effector will be added.
#     - `effector`: The dynamic effector model to be added.
#     """
#     model.dynamic_effectors = (model.dynamic_effectors..., effector)
# end
# ##########################
# ### Thruster functions ###
# ##########################
# function add_thruster!(model::SpacecraftModel, link::Link, thruster::Thruster)
#     """
#     Adds a thruster to the Link
#     - `link` : The Link to which the thruster is added
#     - `thruster` : The thruster to add
#     """
#     # Add the thruster torque to the Jacobian matrix
#     normalize!(thruster.direction) # Ensure the thruster direction is a unit vector
#     # TODO: The if statement here to add thrusters to the jacobian matrix may not be necessary. 
#     # This can instead be calculated at each step in Complete_passage since it will have to account for articulated joints
#     if isempty(link.thrusters)
#         link.J_thruster .= cross(thruster.location, thruster.direction) # If this is the first thruster to be added, simply set the Jacobian equal to the r x F vector
#     else
#         link.J_thruster = hcat(link.J_thruster, cross(thruster.location, thruster.direction)) # If there are already thrusters in the link, append the new r x F vector to the Jacobian
#     end

#     # Append the thruster to the list of thrusters for future reference
#     push!(link.thrusters, thruster)
#     model.n_thrusters += 1 # Increment the number of thrusters in the spacecraft model

# end

# function update_thrusters!(link::Link, torque::AbstractVector{Float64}, t::Float64)
#     """
#     Updates the thrusters of the Link with the given torque vector.
#     - `link`: The Link whose thrusters are to be updated.
#     - `torque`: The torque vector to apply to the thrusters.
#     """
#     # Update the thruster Jacobian matrix to account for possible articulated joints
#     link.J_thruster = MMatrix{3, length(link.thrusters), Float64}(zeros(3, length(link.thrusters))) # Reset the Jacobian matrix
#     rot_to_body = rotate_to_body(link) # Get the rotation matrix to convert from inertial to body frame
#     # println("Rotation matrix: $rot_to_body")
#     for (i, thruster) in enumerate(link.thrusters)
#         normalize!(thruster.direction) # Ensure the thruster direction is a unit vector
#         link.J_thruster[:, i] = cross(rot_to_body * thruster.location + link.r, rot_to_body * thruster.direction) # Update the Jacobian matrix with the r x F vector in the body frame
#     end
#     # Calculate the thrust vector from the torque vector
#     # println("J_thruster: $(link.J_thruster)")
#     # println("Torque req: $torque")
#     thrust_vector = pinv(link.J_thruster) * torque # Solve for the thrust vector using the Jacobian matrix
#     thrust_vector .-= minimum(thrust_vector) # Ensure no negative thrust values
#     # data = CSV.read("basilisk_thruster_force.csv", DataFrame)
#     # time = data[:, 1]
#     # idx = findmin(abs.(time .- t))[2]
#     # thrust_vector = data[idx, 2:end] # Get the thrust vector from the CSV file
#     # println("Thrust vector: $thrust_vector")
#     for (i, thruster) in enumerate(link.thrusters)
#         # println("thruster $i: Requested thrust: $(thrust_vector[i]) N")
#         thruster.thrust = thrust_vector[i] # Update the requested thrust in the thruster
#         link.thrust_calculation_function(link, thruster, thrust_vector[i], t) # Call the function to calculate the average thrust over the control period
#     end
# end

# ####################################
# ### Thrust Calculation functions ###
# ####################################
# function thrust_calculation_schmitt_trigger!(link::Link, thruster::Thruster, thrust::Float64, time::Float64)
#     """
#     Applies a Schmitt trigger to the thrusters of the Link.
#     - `link`: The Link whose thrusters are to be updated.
#     - `thruster`: The thruster to apply the thrust to.
#     - `thrust`: The thrust value to apply to the thruster.
#     """
#     # Interate over each thrust value and determine the firing time
#     ti = min(thrust / thruster.max_thrust * link.attitude_control_rate, link.attitude_control_rate) # Calculate the time interval for the thrust
#     if ti < thruster.min_firing_time # If the time interval is less than the minimum firing time, set it to the minimum firing time or 0, based on schmitt trigger
#         # println("Schmitt triggering thrust time: $ti")
#         ti = schmitt_trigger(ti, thruster.level_on, thruster.level_off) * thruster.min_firing_time # Use Schmitt trigger to determine if the thruster should fire
#         # println("Updated ti: $ti")
#     end
#     # println("TI: $ti")
#     # Integrate the thrust over the time interval to determine the total impulse, accounting for ramp-up/ramp-down
#     CSV.write("thruster_debug.csv", DataFrame(time=time, on_time_request=ti, thrust_req=thrust), append=true)   
#     total_integrated_thrust = integrate_impulse!(link, thruster, ti, time) # Integrate the impulse over the time interval
#     # thruster.stop_firing_time = time + ti # Set the stop firing time to the current time plus the control period. This is used in the next timestep to determine how long to continue firing the thruster into the next timestep
#     thruster.thrust = total_integrated_thrust / link.attitude_control_rate # Update the average thrust value in the thruster
#     # Update the thruster state with the integrated impulse
#     # thruster.thrust = thruster.max_thrust # Update the average thrust value in the thruster
# end
# #####################################
# #####################################
# #####################################

# function create_facet_list(area_list::Vector{Float64}, attitude_list::Vector{SVector{4, Float64}}, normal_vector_list::Vector{SVector{3, Float64}}, 
#                             cp_loc_list::Vector{SVector{3, Float64}}, diffuse_coeffs_list::Vector{Float64}, specular_coeffs_list::Vector{Float64}, facet_names_list::Vector{String})
#     list_length = length(area_list)
#     lists = Vector[area_list, attitude_list, normal_vector_list, cp_loc_list, specular_coeffs_list, diffuse_coeffs_list]
#     @assert all(l -> length(l) == list_length, lists) "Not all lists are the same length"
#     facet_vector = Vector{Facet}(undef, list_length)
#     for i in eachindex(area_list)
#         facet_vector[i] = Facet(area_list[i], attitude_list[i], normal_vector_list[i], cp_loc_list[i], diffuse_coeffs_list[i], specular_coeffs_list[i], facet_names_list[i])
#     end
#     return facet_vector
# end

# function traverse_bodies(model::SpacecraftModel, body::Link)
#     """
#     Traverses the spacecraft model starting from the given body and returns all bodies connected to it.
#     - `model`: The spacecraft model.
#     - `body`: The body from which to start the traversal.
#     """
#     visited = Set{Link}() # Set to keep track of visited bodies
#     root_index = 0 # Index of the root body, if needed
#     queue = [body] # Initialize the queue with the starting body
#     push!(visited, body) # Mark the initial body as visited
#     while !isempty(queue)
#         current_body = popfirst!(queue) # Get the next body in the queue
#         if current_body.root
#             root_index = findfirst(isequal(current_body), model.roots) # Find the index of the root body
#         end
#         # Add children bodies to the queue
#         for joint in model.joints
#             if joint.link1 == current_body && !in(joint.link2, visited)
#                 push!(queue, joint.link2) # Add the second link of the joint to the queue
#                 push!(visited, joint.link2) # Mark it as visited
#             elseif joint.link2 == current_body && !in(joint.link1, visited)
#                 push!(queue, joint.link1) # Add the first link of the joint to the queue
#                 push!(visited, joint.link1) # Mark it as visited
#             end
#         end
#     end
#     @assert root_index != 0 "Root body not found in the model" # Ensure a root body was found
#     return collect(visited), root_index # Return all visited bodies as a vector
# end
# function get_COM(model::SpacecraftModel, body::Link)
#     """
#     Returns the center of mass of the spacecraft that the body is a part of.
#     - `model`: The spacecraft model.
#     - `body`: The body for which to get the center of mass.
#     """
#     # Initialize COM to zero
#     COM = MVector{3, Float64}(0.0, 0.0, 0.0)
#     total_mass = 0.0
#     # BFS starting from body to find all bodies attached to the current body
#     bodies, root_index = traverse_bodies(model, body)
#     println("Bodies found: ", length(bodies))
#     COM = sum([b.r * b.m for b in bodies]) # Sum the position vectors weighted by mass
#     total_mass = sum([b.m for b in bodies]) # Sum the mass of all bodies
#     return COM / total_mass # Return the center of mass
# end

# function get_COM(bodies::Vector{Link})
#     """
#     Returns the center of mass of a collection of bodies.
#     - `bodies`: A vector of bodies for which to calculate the center of mass.
#     """
#     COM = MVector{3, Float64}(0.0, 0.0, 0.0)
#     total_mass = 0.0
#     for body in bodies
#         COM += body.r * body.m # Sum the position vectors weighted by mass
#         total_mass += body.m # Sum the mass
#     end
#     return COM / total_mass # Return the center of mass
# end

# # Function to get the inertia tensor of a body in the spacecraft model
# function update_inertia_tensor!(model::SpacecraftModel, body::Link)
#     """
#     Returns the inertia tensor of the spacecraft that the body is a part of.
#     - `model`: The spacecraft model.
#     - `body`: The body for which to get the inertia tensor.
#     """
#     # Initialize inertia tensor to zero
#     inertia_tensor = SMatrix{3, 3, Float64}(zeros(3, 3))
#     # total_mass = 0.0
#     # BFS starting from body to find all bodies attached to the current body
#     bodies, root_index = traverse_bodies(model, body)
#     for b in bodies
#         # Apply parallel axis theorem to update inertia tensor
#         R = b.root ? SMatrix{3, 3, Float64}(I(3)) : rot(b.q) # Rotation matrix from quaternion
#         I_body = R * b.inertia * R' # Transform inertia tensor to the body frame
#         r = SVector{3, Float64}(b.r) # Position vector of the body
#         fuel_mass = b.root ? model.prop_mass[root_index] : 0.0 # Get propellant mass if root body
#         inertia_tensor += I_body + (b.m+fuel_mass) * hat(r) * hat(r)' # Parallel axis theorem
#     end
#     if root_index > length(model.inertia_tensors)
#         # If the root index is out of bounds, extend the inertia tensors vector
#         push!(model.inertia_tensors, inertia_tensor) # Add the inertia tensor for the root body
#     else
#         model.inertia_tensors[root_index] = inertia_tensor # Update the inertia tensor for the root body
#     end
#     return inertia_tensor # Return the inertia tensor
# end

# function update_inertia_tensor(bodies::Vector{Link})
#     """
#     Returns the inertia tensor of a collection of bodies.
#     - `bodies`: A vector of bodies for which to calculate the inertia tensor.
#     """
#     inertia_tensor = SMatrix{3, 3, Float64}(zeros(3, 3))
#     total_mass = 0.0
#     for body in bodies
#         # Apply parallel axis theorem to update inertia tensor
#         R = body.root ? SMatrix{3, 3, Float64}(I(3)) : rot(body.q) # Rotation matrix from quaternion
#         I_body = R * body.inertia * R' # Transform inertia tensor to the body frame
#         # fuel_mass = body.root ? body.m : 0.0 # Get propellant mass if root body
#         inertia_tensor += I_body + body.m * hat(body.r) * hat(body.r)' # Parallel axis theorem
#         total_mass += body.m # Sum the mass
#     end
#     return inertia_tensor # Return the inertia tensor
# end

# function get_inertia_tensor(model::SpacecraftModel, body::Link)
#     """
#     Returns the inertia tensor of the spacecraft that the body is a part of, in the body frame.
#     - `model`: The spacecraft model.
#     - `body`: The body for which to get the inertia tensor.
#     """
#     # Initialize inertia tensor to zero
#     # inertia_tensor = SMatrix{3, 3, Float64}(zeros(3, 3))
#     # total_mass = 0.0
#     # BFS starting from body to find all bodies attached to the current body
#     bodies, root_index = traverse_bodies(model, body)
#     return model.inertia_tensors[root_index] # Return the inertia tensor of the root body    
# end

# function get_inertia_tensor(model::SpacecraftModel, root_index::Int)
#     """
#     Returns the inertia tensor of the spacecraft that the body is a part of.
#     - `model`: The spacecraft model.
#     - `root_index`: The index of the root body for which to get the inertia tensor.
#     """
#     @assert root_index <= length(model.inertia_tensors) "Root index out of bounds"
#     return model.inertia_tensors[root_index] # Return the inertia tensor of the root body
# end

# function set_inertia_tensor!(model::SpacecraftModel, body::Link, inertia_tensor::SMatrix{3, 3, Float64})
#     """
#     Sets the inertia tensor of a body in the spacecraft model.
#     - `model`: The spacecraft model.
#     - `body`: The body for which to set the inertia tensor.
#     - `inertia_tensor`: The inertia tensor to set.
#     """
#     # Find the index of the body in the links vector
#     bodies, root_index = traverse_bodies(model, body)
#     # @assert index !== nothing "Body not found in the model"
#     model.inertia_tensors[root_index] = inertia_tensor # Set the inertia tensor for the body
# end


# function get_spacecraft_mass(model::SpacecraftModel, body::Link; dry=false)
#     """
#     Calculates the total mass of the spacecraft model by summing the masses of all bodies connected to body.
#     This is useful for dynamics calculations.
#     - `model`: The spacecraft model.
#     - `body`: The body for which to calculate the total mass.
#     """
#     total_mass = 0.0
#     # BFS starting from body to find all bodies attached to the current body
#     bodies, root_index = traverse_bodies(model, body)
#     for b in bodies
#         total_mass += b.m # Sum the mass of each body
#     end
#     return dry ? total_mass : total_mass + model.prop_mass[root_index] # Return the total mass
# end

# function get_spacecraft_mass(model::SpacecraftModel, bodies::Vector{Link}, root_index::Int;dry=false)
#     """
#     Calculates the total mass of a collection of bodies by summing their masses.
#     - `bodies`: A vector of bodies for which to calculate the total mass.
#     """
#     total_mass = sum([b.m for b in bodies]) # Sum the mass of each body
#     return dry ? total_mass : total_mass + model.prop_mass[root_index] # Return the total mass
# end

# function get_spacecraft_mass(model::SpacecraftModel;dry=false)
#     """
#     Calculates the total mass of the spacecraft model by summing the masses of all bodies.
#     This is useful for dynamics calculations.
#     - `model`: The spacecraft model.
#     """
#     masses = Float64[]
#     for root in model.roots
#         # BFS starting from root to find all bodies attached to the root body
#         bodies, root_index = traverse_bodies(model, root)
#         total_mass = sum([b.m for b in bodies]) # Sum the mass of each body
#         push!(masses, dry ? total_mass : total_mass + model.prop_mass[root_index]) # Return the total mass
#     end
#     return length(model.roots) == 1 ? masses[1] : masses # Return the total mass of all roots
# end

# function get_spacecraft_mass(model::SpacecraftModel, bodies::Vector{Link}, root_index::Int; dry=false)
#     """
#     Calculates the total mass of a collection of bodies by summing their masses.
#     - `bodies`: A vector of bodies for which to calculate the total mass.
#     """
#     total_mass = sum([b.m for b in bodies]) # Sum the mass of each body
#     return dry ? total_mass : total_mass + model.prop_mass[root_index] # Return the total mass
# end

# function get_spacecraft_reference_area(model::SpacecraftModel)
#     """
#     Calculates the total reference area of the spacecraft model by summing the areas of all bodies.
#     This is useful for aerodynamic calculations.
#     - `model`: The spacecraft model.
#     """
#     total_area = Float64[]
#     for root in model.roots
#         area = 0.0
#         # BFS starting from root to find all bodies attached to the root body
#         bodies, root_index = traverse_bodies(model, root)
#         for body in bodies
#             if body.ref_area > 0.0 # Only consider bodies with defined area
#                 area += body.ref_area # Sum the reference areas of each body
#             end
#         end
#         push!(total_area, area) # Add the total area for this root
#     end
#     return length(total_area) == 1 ? total_area[1] : total_area # Return the total area of all roots
# end

# function get_spacecraft_reference_area(model::SpacecraftModel, body::Link)
#     """
#     Calculates the total reference area of a collection of bodies by summing their areas.
#     - `bodies`: A vector of bodies for which to calculate the total reference area.
#     """
#     total_area = 0.0
#     # BFS starting from body to find all bodies attached to the current body
#     bodies, root_index = traverse_bodies(model, body)
#     for body in bodies
#         if body.ref_area > 0.0 # Only consider bodies with defined area
#             total_area += body.ref_area # Sum the reference areas of each body
#         end
#     end
#     return total_area # Return the total area
# end

# function get_spacecraft_reference_area(bodies::Vector{Link})
#     """
#     Calculates the total reference area of a collection of bodies by summing their areas.
#     - `bodies`: A vector of bodies for which to calculate the total reference area.
#     """
#     total_area = 0.0
#     for body in bodies
#         if body.ref_area > 0.0 # Only consider bodies with defined area
#             total_area += body.ref_area # Sum the reference areas of each body
#         end
#     end
#     return total_area # Return the total area
# end

# function get_spacecraft_length(model::SpacecraftModel, body::Link)
#     """
#     Calculates the total length of the spacecraft model by summing the lengths of all bodies.
#     This is useful for determining the overall size of the spacecraft.
#     - `model`: The spacecraft model.
#     - `body`: The body for which to calculate the total length.
#     """
#     max_length = 0.0
#     # BFS starting from body to find all bodies attached to the current body
#     bodies, root_index = traverse_bodies(model, body)
#     for b in bodies
#         if b.dims[1] > max_length # Use the maximum dimension as the length
#             max_length = b.dims[1] # Update the maximum length
#         end
#     end
#     return max_length # Return the maximum length found
# end

# function get_SA_area(model::SpacecraftModel, body::Link)
#     """
#     Calculates the total surface area of the solar array by summing the areas of all flat plates.
#     This is useful for aerodynamic calculations.
#     - `model`: The spacecraft model.
#     - `body`: The body for which to calculate the total surface area of the solar array.
#     """
#     total_area = 0.0
#     # BFS starting from body to find all bodies attached to the current body
#     bodies, root_index = traverse_bodies(model, body)
#     for b in bodies
#         if !b.root # Only consider flat plates
#             total_area += b.ref_area # Sum the areas of flat plates
#         end
#     end
#     return total_area # Return the total area of the solar array
# end

# function get_SA_area(bodies::Vector{Link})
#     """
#     Calculates the total surface area of a collection of bodies by summing the areas of all flat plates.
#     - `bodies`: A vector of bodies for which to calculate the total surface area of the solar array.
#     """
#     total_area = 0.0
#     for b in bodies
#         if !b.root # Only consider flat plates
#             total_area += b.ref_area # Sum the areas of flat plates
#         end
#     end
#     return total_area # Return the total area of the solar array
# end

# function get_SC_area(model::SpacecraftModel, body::Link)
#     """
#     Calculates the total surface area of the spacecraft bus by summing the areas of all boxes.
#     This is useful for aerodynamic calculations.
#     - `model`: The spacecraft model.
#     - `body`: The body for which to calculate the total surface area of the spacecraft bus.
#     """
#     total_area = 0.0
#     # BFS starting from body to find all bodies attached to the current body
#     bodies, root_index = traverse_bodies(model, body)
#     for b in bodies
#         if b.root # Only consider boxes
#             total_area += b.ref_area # Sum the areas of boxes
#         end
#     end
#     return total_area # Return the total area of the spacecraft bus
# end

# function get_SC_area(bodies::Vector{Link})
#     """
#     Calculates the total surface area of a collection of bodies by summing the areas of all boxes.
#     - `bodies`: A vector of bodies for which to calculate the total surface area of the spacecraft bus.
#     """
#     total_area = 0.0
#     for b in bodies
#         if b.root # Only consider boxes
#             total_area += b.ref_area # Sum the areas of boxes
#         end
#     end
#     return total_area # Return the total area of the spacecraft bus
# end

# function get_normal_vector(model::SpacecraftModel, body::Link, root_index::Int; normalized=false)
#     """
#     Returns the normal vector of the body in the inertial frame.
#     - `body`: The body for which to get the normal vector.
#     """
#     # Assuming the body is a flat plate, the normal vector is along the z-axis in the body frame
#     # Convert to inertial frame using the orientation quaternion
#     R = rotate_to_inertial(model, body, root_index) # Get the rotation matrix to convert from body frame to inertial frame

#     # Return the normal vector in the inertial frame
#     if normalized
#         normal = R * SVector{3, Float64}(1.0, 0.0, 0.0) # Normal vector in inertial frame
#         return normal / norm(normal) # Normalize the vector
#     else
#         return R * SVector{3, Float64}(1.0, 0.0, 0.0) # Normal vector in inertial frame
#     end
# end

# function get_tangent_vector(model::SpacecraftModel, body::Link, root_index::Int; normalized=false)
#     """
#     Returns the tangent vector of the body in the inertial frame (CN direction, up).
#     - `body`: The body for which to get the normal vector.
#     """
#     # Assuming the body is a flat plate, the normal vector is along the z-axis in the body frame
#     # Convert to inertial frame using the orientation quaternion
#     R = rotate_to_inertial(model, body, root_index) # Get the rotation matrix to convert from body frame to inertial frame

#     # Return the normal vector in the inertial frame
#     if normalized
#         normal = R * SVector{3, Float64}(0.0, 0.0, 1.0) # Normal vector in inertial frame
#         return normal / norm(normal) # Normalize the vector
#     else
#         return R * SVector{3, Float64}(0.0, 0.0, 1.0) # Normal vector in inertial frame
#     end
# end

# function rotate_to_inertial(model::SpacecraftModel, body::Link, root_index::Int)
#     """
#     Returns the rotation matrix to convert from the link frame to the inertial frame.
#     - `model`: The spacecraft model.
#     - `body`: The body for which to get the rotation matrix.
#     """
#     if body.root
#         return rot(body.q)' # Rotation matrix from quaternion
#     else
#         return rot(model.roots[root_index].q)' * rot(body.q)' # Rotation matrix from quaternion
#     end
# end

# function rotate_to_body(body::Link)
#     """
#     Returns the rotation matrix to convert from the link frame to the body frame
#     - `model` : The spacecraft model.
#     - `body` : The body for which to get the rotation matrix
#     """
#     if body.root
#         return I(3)
#     else
#         return rot(body.q)'
#     end
# end
# function rotate_link(body::Link, q::SVector{4, Float64})
#     """
#     Rotates the link to a new orientation defined by the quaternion `q`.
#     - `body`: The body to be rotated.
#     - `root_index`: The index of the root body in the model.
#     - `q`: The new orientation quaternion.
#     """
#     @assert !body.root "Cannot rotate a root body directly"
#     # Update the orientation of the body
#     body.q .= q # Set the new orientation quaternion
# end

# function rotate_link(body::Link, dcm::SMatrix{3, 3, Float64})
#     """
#     Rotates the link to a new orientation defined by the direction cosine matrix `dcm`.
#     - `body`: The body to be rotated.
#     - `dcm`: The direction cosine matrix representing the new orientation.
#     """
#     @assert !body.root "Cannot rotate a root body directly"
#     # Update the orientation of the body
#     body.q .= dcm_to_quaternion(dcm)
# end

# function rotate_link(body::Link, axis::SVector{3, Float64},  θ::Float64)
#     """
#     Rotates the link to a new orientation defined by the Euler angles `θ`.
#     - `body`: The body to be rotated.
#     - `axis`: The rotation axis as a 3D vector.
#     - `θ`: The rotation angle in radians.

#     """
#     @assert !body.root "Cannot rotate a root body directly"
#     if norm(axis) <= 1e-6
#         # @warn "Rotation axis norm is too small, using default axis (0, 1, 0)"
#         axis = SVector{3, Float64}(0.0, 1.0, 0.0) # Default axis if norm is too small
#     end
#     axis = axis / norm(axis) # Normalize the rotation axis
#     # Update the orientation of the body
#     body.q .= SVector{4, Float64}([axis .* sin(θ/2); cos(θ/2)]) # Convert Euler angles to quaternion
# end

# function copy(model::SpacecraftModel)
#     """
#     Creates a deep copy of the spacecraft model.
#     - `model`: The spacecraft model to be copied.
#     """
#     new_model = SpacecraftModel(
#         joints = Joint[Joint(j.link1, j.p1ᵇ, j.link2, j.p2ᵇ, j.Kx, j.Kt, j.Cx, j.Ct) for j in model.joints],
#         links = Link[Link(j) for j in model.links],
#         roots = Link[Link(r) for r in model.roots],
#         instant_actuation = model.instant_actuation,
#         prop_mass = copy(model.prop_mass),
#         inertia_tensors = [copy(it) for it in model.inertia_tensors],
#         n_reaction_wheels = model.n_reaction_wheels
#     )
#     return new_model # Return the deep copy of the spacecraft model
# end

# #################################
# ### Thruster helper functions ###
# #################################
# function schmitt_trigger(input::Float64, level_on::Float64, level_off::Float64)
#     """
#     Implements a Schmitt trigger.
#     - `input`: The input signal to the Schmitt trigger.
#     - `level_on`: The threshold above which the output is activated.
#     - `level_off`: The threshold below which the output is deactivated.
#     """
#     state = 0.0
#     if input > level_on
#         state = 1.0
#     elseif input < level_off
#         state = 0.0
#     end
#     return state
# end

# function integrate_impulse!(link::Link, thruster::Thruster, on_time_request::Float64, time::Float64)
#     """
#     Integrates the impulse and thrust factor over the on-time request period.
#     """
#     ω = thruster.cutoff_frequency # Get the cutoff frequency of the thruster
#     # Calculate the impulse using the first-order filter
#     κ = thruster.κ # Get the current thrust factor
#     # carryover_time = max(0.0, thruster.stop_firing_time - time) # Calculate the carryover time from the previous timestep
#     # # if carryover_time > 0.0
#     # impulse = thruster.max_thrust * κ / ω * (1 - exp(-ω * 0.1)) # Initial non-firing time
#     # κ *= exp(-ω * 0.1) # Update the thrust factor after initial non-firing time

#     # on_time = min(on_time_request, link.attitude_control_rate) # Ensure the on-time does not exceed the control period
#     # on_time_request -= 0.001
#     # println("On-time request: $on_time_request s, κ: $κ")
#     total_integrated_thrust = thruster.max_thrust * (on_time_request + (κ - 1)/ω * (1-exp(-ω * on_time_request))) # Calculate the total impulse
#     κ = 1 + (κ - 1) * exp(-ω * on_time_request) # Calculate the final thrust factor
#     if on_time_request < link.attitude_control_rate
#         total_integrated_thrust += thruster.max_thrust * κ / ω * (1 - exp(-ω * (link.attitude_control_rate - on_time_request))) # Add the impulse from ramp-down if applicable
#         κ *= exp(-ω * (link.attitude_control_rate - on_time_request)) # Update the final thrust factor after ramp-down
#     end
#     # else
#     #     total_integrated_thrust = thruster.max_thrust * κ / ω * (1 - exp(-ω * (link.attitude_control_rate))) # If no on-time, just do ramp-down
#     #     κ *= exp(-ω * (link.attitude_control_rate)) # Update the final thrust factor after ramp-down
#     #     println("No on-time, ramping down")
#     # end
    
#     # if κ < 1.0e-10
#     #     κ = 0.0 # Set to zero if below threshold
#     # end
#     thruster.κ = κ # Update the thrust factor in the thruster
#     return total_integrated_thrust # Return the total impulse integrated over the control period
# end