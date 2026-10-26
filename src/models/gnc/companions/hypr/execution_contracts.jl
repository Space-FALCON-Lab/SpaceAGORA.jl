# Compatibility function identities only. The optional companion owns execution.
# HYPR's own helpers are defined in its extension, not declared here.
using ..HYPRSupport
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_adaptive_sampling_min_ds_m end
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_fuel_proxy_from_samples end
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_normalized_path_cost_components end
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_path_cost_normalization_refs end
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_post_refine_path end
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_pso_plan_path end
rpo_pso_plan_path(args...; kwargs...) = HYPRSupport.unavailable(rpo_pso_plan_path, args)
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_reference_from_path end
rpo_reference_from_path(args...; kwargs...) = HYPRSupport.unavailable(rpo_reference_from_path, args)
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_retime_profile end
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_retime_sampling_ds_m end
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_rrt_connect_bezier_plan_path end
rpo_rrt_connect_bezier_plan_path(args...; kwargs...) = HYPRSupport.unavailable(rpo_rrt_connect_bezier_plan_path, args)
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_rrt_connect_plan_path end
rpo_rrt_connect_plan_path(args...; kwargs...) = HYPRSupport.unavailable(rpo_rrt_connect_plan_path, args)
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_rrt_star_plan_path end
rpo_rrt_star_plan_path(args...; kwargs...) = HYPRSupport.unavailable(rpo_rrt_star_plan_path, args)
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_sample_path end
rpo_sample_path(args...; kwargs...) = HYPRSupport.unavailable(rpo_sample_path, args)
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_sample_path_bezier_adaptive end
rpo_sample_path_bezier_adaptive(args...; kwargs...) = HYPRSupport.unavailable(rpo_sample_path_bezier_adaptive, args)
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_sample_path_bezier_adaptive_with_params end
rpo_sample_path_bezier_adaptive_with_params(args...; kwargs...) = HYPRSupport.unavailable(rpo_sample_path_bezier_adaptive_with_params, args)
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_sample_path_polyline_adaptive end
rpo_sample_path_polyline_adaptive(args...; kwargs...) = HYPRSupport.unavailable(rpo_sample_path_polyline_adaptive, args)
"""HYPR compatibility entry point; load `SpaceAGORA` and `HYPR` for execution."""
function rpo_sample_path_with_params end
rpo_sample_path_with_params(args...; kwargs...) = HYPRSupport.unavailable(rpo_sample_path_with_params, args)
