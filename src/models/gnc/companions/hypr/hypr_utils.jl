"""Shared path and RRT helpers for the optional HYPR planners."""
module HYPRUtils

using LinearAlgebra

export hypr_path_length, hypr_bezier_point, hypr_bezier_point!, hypr_sample_count_path
export hypr_rrt_nearest_index, hypr_rrt_near_indices, hypr_rrt_steer
export hypr_rrt_tree_path, hypr_rrt_join_paths, hypr_rrt_refresh_subtree_costs!

# Keep the compatibility module and exports; source owners are shared paths
# and RRT tree operations. HYPR owns its swarm policy.
include(joinpath(@__DIR__, "..", "..", "planning", "common", "path_geometry.jl"))
include(joinpath(@__DIR__, "..", "..", "planning", "common", "rrt_tree_operations.jl"))

end
