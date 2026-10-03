# Compatibility function identities only. The optional companion owns execution.
# Private names remain available qualified for existing research clients.
using ..HYPRSupport
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _rpo_rrt_configured_result end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function _rpo_sampling_settings end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_adaptive_pso_config end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_adaptive_sampling_min_ds_m end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_estimate_geometry_complexity end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_fit_bezier_fixed_endpoints end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_fuel_proxy_dt_s end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_fuel_proxy_from_samples end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_manuscript_adaptive_pso_config end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_manuscript_exploration_score end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_manuscript_path_cost_components end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_normalized_path_cost_components end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_path_cost end
rpo_path_cost(args...; kwargs...) = HYPRSupport.unavailable(rpo_path_cost, args)
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_path_cost_normalization_refs end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_position_to_path end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_post_refine_path end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_probe_geometry_metrics end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_bounds end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_cull_swarm! end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_early_stopping_feasible end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_effective_safe_distance end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_empty_warmstart_diagnostics end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_iteration_weights end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_material_improvement end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_plan_path end
rpo_pso_plan_path(args...; kwargs...) = HYPRSupport.unavailable(rpo_pso_plan_path, args)
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_project_to_segment end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_protected_particle_mask end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_rrt_warmstart_path end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_stagnation_count_after_learning end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_station_bounds end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_tapered_noise_scale end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_pso_warmstart_bounds end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_reference_from_path end
rpo_reference_from_path(args...; kwargs...) = HYPRSupport.unavailable(rpo_reference_from_path, args)
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_refine_lower_degree end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_refine_shortcut_refit end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_refine_tighten_handles end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_refinement_bernstein end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_refinement_better end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_refinement_clamp_path end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_refinement_config end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_refinement_project_to_segment end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_refinement_sample_params end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_refinement_segment_is_safe end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_refinement_segment_samples end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_refinement_shortcut_samples end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_retime_available_distance end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_retime_path end
rpo_retime_path(args...; kwargs...) = HYPRSupport.unavailable(rpo_retime_path, args)
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_retime_pointwise_speed end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_retime_profile end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_retime_sampling_ds_m end
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_retimed_reference end
rpo_retimed_reference(args...; kwargs...) = HYPRSupport.unavailable(rpo_retimed_reference, args)
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_rrt_connect_bezier_plan_path end
rpo_rrt_connect_bezier_plan_path(args...; kwargs...) = HYPRSupport.unavailable(rpo_rrt_connect_bezier_plan_path, args)
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_rrt_connect_plan_path end
rpo_rrt_connect_plan_path(args...; kwargs...) = HYPRSupport.unavailable(rpo_rrt_connect_plan_path, args)
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_rrt_star_plan_path end
rpo_rrt_star_plan_path(args...; kwargs...) = HYPRSupport.unavailable(rpo_rrt_star_plan_path, args)
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_sample_path end
rpo_sample_path(args...; kwargs...) = HYPRSupport.unavailable(rpo_sample_path, args)
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_sample_path_bezier_adaptive end
rpo_sample_path_bezier_adaptive(args...; kwargs...) = HYPRSupport.unavailable(rpo_sample_path_bezier_adaptive, args)
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_sample_path_bezier_adaptive_with_params end
rpo_sample_path_bezier_adaptive_with_params(args...; kwargs...) = HYPRSupport.unavailable(rpo_sample_path_bezier_adaptive_with_params, args)
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_sample_path_polyline_adaptive end
rpo_sample_path_polyline_adaptive(args...; kwargs...) = HYPRSupport.unavailable(rpo_sample_path_polyline_adaptive, args)
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_sample_path_with_params end
rpo_sample_path_with_params(args...; kwargs...) = HYPRSupport.unavailable(rpo_sample_path_with_params, args)
"""HYPR compatibility entry point; load `SpaceAGORAHYPR` for execution."""
function rpo_try_accept_refinement end
