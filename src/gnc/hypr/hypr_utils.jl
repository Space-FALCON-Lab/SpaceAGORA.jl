"""Shared path/RRT helpers and compatibility bindings for optional HYPR swarm policy."""
module HYPRUtils

using LinearAlgebra

export hypr_path_length, hypr_bezier_point, hypr_bezier_point!, hypr_sample_count_path
export hypr_iteration_weights, hypr_material_improvement, hypr_protected_particle_mask
export hypr_rrt_nearest_index, hypr_rrt_near_indices, hypr_rrt_steer
export hypr_rrt_tree_path, hypr_rrt_join_paths, hypr_rrt_refresh_subtree_costs!

# Keep the compatibility module and exports; source owners are shared paths,
# RRT tree operations and HYPR swarm policy.
include(joinpath(@__DIR__, "..", "shared", "path_geometry.jl"))
include(joinpath(@__DIR__, "..", "rrt", "tree_operations.jl"))


using ..HYPRSupport
function hypr_iteration_weights end
function hypr_material_improvement end
function hypr_protected_particle_mask end

end
