module Structure

using StaticArrays
using LinearAlgebra

using ..SpacecraftModels: SpacecraftModel, Link, Joint
using ..Kinematics: rotate_to_inertial

export traverse_bodies,
       get_COM,
       update_inertia_tensor!,
       update_inertia_tensor,
       get_inertia_tensor,
       set_inertia_tensor!,
       get_spacecraft_mass,
       get_spacecraft_reference_area,
       get_spacecraft_length,
       get_SA_area,
       get_SC_area,
       get_normal_vector,
       get_tangent_vector

include(joinpath(@__DIR__, "assembly_graph.jl"))
include(joinpath(@__DIR__, "mass_properties.jl"))
include(joinpath(@__DIR__, "geometry_properties.jl"))

end # module Structure
