module Structure

using StaticArrays
using LinearAlgebra
using JSON
using Base64
using Random

using ..SpacecraftModels: SpacecraftModel, Link, Joint
using ..Kinematics: rotate_to_inertial, rot

export get_COM,
       get_inertia_tensor,
       get_spacecraft_reference_area,
       get_SA_area,
       get_SC_area,
       get_normal_vector,
       get_tangent_vector,
       load_model_triangles,
       articulate_triangles,
       articulation_payload,
       model_bounding_box,
       model_bounding_box_center,
       sample_model_pointcloud,
       model_format,
       gltf_required_extensions,
       GLTF_UNSUPPORTED_REQUIRED

include(joinpath(@__DIR__, "mass_properties.jl"))
include(joinpath(@__DIR__, "geometry_properties.jl"))
include(joinpath(@__DIR__, "mesh_geometry.jl"))

end # module Structure
