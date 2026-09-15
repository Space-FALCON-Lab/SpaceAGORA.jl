module Structure

using StaticArrays
using LinearAlgebra
using JSON
using Base64
using Random

using ..SpacecraftModels: SpacecraftModel, Link, Joint
using ..Kinematics: rotate_to_inertial, rot

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
       get_tangent_vector,
       load_model_triangles,
       model_bounding_box,
       model_bounding_box_center,
       sample_model_pointcloud,
       model_format,
       gltf_required_extensions,
       GLTF_UNSUPPORTED_REQUIRED

include(joinpath(@__DIR__, "assembly_graph.jl"))
include(joinpath(@__DIR__, "mass_properties.jl"))
include(joinpath(@__DIR__, "geometry_properties.jl"))
include(joinpath(@__DIR__, "mesh_geometry.jl"))

end # module Structure
