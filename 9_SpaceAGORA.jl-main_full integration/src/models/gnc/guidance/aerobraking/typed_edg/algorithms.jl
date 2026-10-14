# Private compatibility bindings for existing core and legacy callers.
# The standalone companion owns every numerical kernel in this inventory.
module EDGAlgorithms
using ..EDGServices: _invoke_kernel, EDGControlDecision
const COMPATIBILITY_FUNCTIONS = (
    :_edg_at_control_tick,
    :_edg_balanced_tpbvp_heat_load_window,
    :_edg_base_alpha,
    :_edg_candidate_constraint_sample,
    :_edg_candidate_energy_mapping,
    :_edg_certify_targeting_candidates,
    :_edg_chebyshev4,
    :_edg_chebyshev44,
    :_edg_chebyshev_t0_t4,
    :_edg_closed_form_correction_fit,
    :_edg_closed_form_density,
    :_edg_closed_form_density_fit,
    :_edg_closed_form_f1_f2,
    :_edg_closed_form_heat_load_trajectory,
    :_edg_command_alpha!,
    :_edg_constrained_heat_load_alpha_profile,
    :_edg_constraint_min_alpha,
    :_edg_constraint_status,
    :_edg_controlled_panel_area,
    :_edg_controlled_ref_area,
    :_edg_corrected_target_energy_from_apoapsis,
    :_edg_disable_uncertified_targeting!,
    :_edg_drag_passage_duration,
    :_edg_dynamic_spacecraft_aero_coefficients,
    :_edg_eccentric_anomaly_from_true,
    :_edg_exact_targeting_exit,
    :_edg_first_low_alpha_interval_indices,
    :_edg_first_two_switch_alpha_profile,
    :_edg_heat_budget_status,
    :_edg_heat_load_alpha_profile,
    :_edg_heat_load_coefficients,
    :_edg_heat_load_costates,
    :_edg_heat_load_lambdas,
    :_edg_heat_load_low_drag_active,
    :_edg_heat_load_profile_for_k,
    :_edg_heat_load_scale_height,
    :_edg_heat_load_security_required,
    :_edg_heat_load_shooting_profile,
    :_edg_heat_load_shooting_stage,
    :_edg_heat_load_track_env,
    :_edg_heat_rate_alpha,
    :_edg_integrate_heat_load_shooting,
    :_edg_integrate_series,
    :_edg_integrated_heat_load_trajectory,
    :_edg_integrated_max_energy_depletion_trajectory,
    :_edg_integrated_targeting_trajectory,
    :_edg_interpolate_bracket_value,
    :_edg_journal_plate_coefficients,
    :_edg_journal_spacecraft_aero_coefficients,
    :_edg_journal_targeting_plan,
    :_edg_low_alpha_switch_window,
    :_edg_max_heat_load_for_links,
    :_edg_maxwellian_heat_rate,
    :_edg_mean_anomaly_from_true,
    :_edg_orbit_metrics_from_rv,
    :_edg_ordered_feasible_prefix,
    :_edg_outbound_exit_index,
    :_edg_padded_heat_load_window,
    :_edg_pass_heat_load_for_links,
    :_edg_predict_mass,
    :_edg_predict_max_energy_depletion_outcome,
    :_edg_predict_targeting_outcome,
    :_edg_prediction_advance,
    :_edg_prediction_event,
    :_edg_prediction_gravity_acceleration,
    :_edg_prediction_step,
    :_edg_prediction_time_grid,
    :_edg_profile_heat_load,
    :_edg_profile_heat_rates,
    :_edg_recompute_second_heat_load_switch,
    :_edg_recompute_switches!,
    :_edg_remaining_closed_form_heat_load_trajectory,
    :_edg_remaining_closed_form_reference,
    :_edg_replan_due,
    :_edg_run_target_energy_bracketing!,
    :_edg_select_mapped_candidate,
    :_edg_set_targeting_fallback!,
    :_edg_solve_heat_load_switches,
    :_edg_solve_targeting_switch,
    :_edg_structural_alpha,
    :_edg_switch_due,
    :_edg_target_energy_from_apoapsis,
    :_edg_target_energy_from_reachable_bracket,
    :_edg_targeting_aero_acceleration,
    :_edg_targeting_bracket_outcomes,
    :_edg_targeting_candidate_certificate,
    :_edg_targeting_constrained_alpha,
    :_edg_targeting_control_period,
    :_edg_targeting_outcome_with_heat_load,
    :_edg_targeting_prediction_time_grid,
    :_edg_targeting_switch_outcomes,
    :_edg_targeting_switch_solution,
    :_edg_total_ref_area,
    :_edg_update_heat_load_security!,
    :_edg_vacuum_apoapsis_correction,
    :_edg_vacuum_drag_passage_exit,
    :_edg_weighted_aero_coefficients,
    :_energy_depletion_heat_rate_calc,
    :_energy_depletion_heatrate_root_alpha,
    :_energy_depletion_struct_drag_area,
    :_energy_depletion_struct_load_root_alpha,
    :control_decision!,
    :guidance_decision!,
)

_edg_at_control_tick(args...; kwargs...) = _invoke_kernel(:_edg_at_control_tick, args...; kwargs...)
_edg_balanced_tpbvp_heat_load_window(args...; kwargs...) = _invoke_kernel(:_edg_balanced_tpbvp_heat_load_window, args...; kwargs...)
_edg_base_alpha(args...; kwargs...) = _invoke_kernel(:_edg_base_alpha, args...; kwargs...)
_edg_candidate_constraint_sample(args...; kwargs...) = _invoke_kernel(:_edg_candidate_constraint_sample, args...; kwargs...)
_edg_candidate_energy_mapping(args...; kwargs...) = _invoke_kernel(:_edg_candidate_energy_mapping, args...; kwargs...)
_edg_certify_targeting_candidates(args...; kwargs...) = _invoke_kernel(:_edg_certify_targeting_candidates, args...; kwargs...)
_edg_chebyshev4(args...; kwargs...) = _invoke_kernel(:_edg_chebyshev4, args...; kwargs...)
_edg_chebyshev44(args...; kwargs...) = _invoke_kernel(:_edg_chebyshev44, args...; kwargs...)
_edg_chebyshev_t0_t4(args...; kwargs...) = _invoke_kernel(:_edg_chebyshev_t0_t4, args...; kwargs...)
_edg_closed_form_correction_fit(args...; kwargs...) = _invoke_kernel(:_edg_closed_form_correction_fit, args...; kwargs...)
_edg_closed_form_density(args...; kwargs...) = _invoke_kernel(:_edg_closed_form_density, args...; kwargs...)
_edg_closed_form_density_fit(args...; kwargs...) = _invoke_kernel(:_edg_closed_form_density_fit, args...; kwargs...)
_edg_closed_form_f1_f2(args...; kwargs...) = _invoke_kernel(:_edg_closed_form_f1_f2, args...; kwargs...)
_edg_closed_form_heat_load_trajectory(args...; kwargs...) = _invoke_kernel(:_edg_closed_form_heat_load_trajectory, args...; kwargs...)
_edg_command_alpha!(args...; kwargs...) = _invoke_kernel(:_edg_command_alpha!, args...; kwargs...)
_edg_constrained_heat_load_alpha_profile(args...; kwargs...) = _invoke_kernel(:_edg_constrained_heat_load_alpha_profile, args...; kwargs...)
_edg_constraint_min_alpha(args...; kwargs...) = _invoke_kernel(:_edg_constraint_min_alpha, args...; kwargs...)
_edg_constraint_status(args...; kwargs...) = _invoke_kernel(:_edg_constraint_status, args...; kwargs...)
_edg_controlled_panel_area(args...; kwargs...) = _invoke_kernel(:_edg_controlled_panel_area, args...; kwargs...)
_edg_controlled_ref_area(args...; kwargs...) = _invoke_kernel(:_edg_controlled_ref_area, args...; kwargs...)
_edg_corrected_target_energy_from_apoapsis(args...; kwargs...) = _invoke_kernel(:_edg_corrected_target_energy_from_apoapsis, args...; kwargs...)
_edg_disable_uncertified_targeting!(args...; kwargs...) = _invoke_kernel(:_edg_disable_uncertified_targeting!, args...; kwargs...)
_edg_drag_passage_duration(args...; kwargs...) = _invoke_kernel(:_edg_drag_passage_duration, args...; kwargs...)
_edg_dynamic_spacecraft_aero_coefficients(args...; kwargs...) = _invoke_kernel(:_edg_dynamic_spacecraft_aero_coefficients, args...; kwargs...)
_edg_eccentric_anomaly_from_true(args...; kwargs...) = _invoke_kernel(:_edg_eccentric_anomaly_from_true, args...; kwargs...)
_edg_exact_targeting_exit(args...; kwargs...) = _invoke_kernel(:_edg_exact_targeting_exit, args...; kwargs...)
_edg_first_low_alpha_interval_indices(args...; kwargs...) = _invoke_kernel(:_edg_first_low_alpha_interval_indices, args...; kwargs...)
_edg_first_two_switch_alpha_profile(args...; kwargs...) = _invoke_kernel(:_edg_first_two_switch_alpha_profile, args...; kwargs...)
_edg_heat_budget_status(args...; kwargs...) = _invoke_kernel(:_edg_heat_budget_status, args...; kwargs...)
_edg_heat_load_alpha_profile(args...; kwargs...) = _invoke_kernel(:_edg_heat_load_alpha_profile, args...; kwargs...)
_edg_heat_load_coefficients(args...; kwargs...) = _invoke_kernel(:_edg_heat_load_coefficients, args...; kwargs...)
_edg_heat_load_costates(args...; kwargs...) = _invoke_kernel(:_edg_heat_load_costates, args...; kwargs...)
_edg_heat_load_lambdas(args...; kwargs...) = _invoke_kernel(:_edg_heat_load_lambdas, args...; kwargs...)
_edg_heat_load_low_drag_active(args...; kwargs...) = _invoke_kernel(:_edg_heat_load_low_drag_active, args...; kwargs...)
_edg_heat_load_profile_for_k(args...; kwargs...) = _invoke_kernel(:_edg_heat_load_profile_for_k, args...; kwargs...)
_edg_heat_load_scale_height(args...; kwargs...) = _invoke_kernel(:_edg_heat_load_scale_height, args...; kwargs...)
_edg_heat_load_security_required(args...; kwargs...) = _invoke_kernel(:_edg_heat_load_security_required, args...; kwargs...)
_edg_heat_load_shooting_profile(args...; kwargs...) = _invoke_kernel(:_edg_heat_load_shooting_profile, args...; kwargs...)
_edg_heat_load_shooting_stage(args...; kwargs...) = _invoke_kernel(:_edg_heat_load_shooting_stage, args...; kwargs...)
_edg_heat_load_track_env(args...; kwargs...) = _invoke_kernel(:_edg_heat_load_track_env, args...; kwargs...)
_edg_heat_rate_alpha(args...; kwargs...) = _invoke_kernel(:_edg_heat_rate_alpha, args...; kwargs...)
_edg_integrate_heat_load_shooting(args...; kwargs...) = _invoke_kernel(:_edg_integrate_heat_load_shooting, args...; kwargs...)
_edg_integrate_series(args...; kwargs...) = _invoke_kernel(:_edg_integrate_series, args...; kwargs...)
_edg_integrated_heat_load_trajectory(args...; kwargs...) = _invoke_kernel(:_edg_integrated_heat_load_trajectory, args...; kwargs...)
_edg_integrated_max_energy_depletion_trajectory(args...; kwargs...) = _invoke_kernel(:_edg_integrated_max_energy_depletion_trajectory, args...; kwargs...)
_edg_integrated_targeting_trajectory(args...; kwargs...) = _invoke_kernel(:_edg_integrated_targeting_trajectory, args...; kwargs...)
_edg_interpolate_bracket_value(args...; kwargs...) = _invoke_kernel(:_edg_interpolate_bracket_value, args...; kwargs...)
_edg_journal_plate_coefficients(args...; kwargs...) = _invoke_kernel(:_edg_journal_plate_coefficients, args...; kwargs...)
_edg_journal_spacecraft_aero_coefficients(args...; kwargs...) = _invoke_kernel(:_edg_journal_spacecraft_aero_coefficients, args...; kwargs...)
_edg_journal_targeting_plan(args...; kwargs...) = _invoke_kernel(:_edg_journal_targeting_plan, args...; kwargs...)
_edg_low_alpha_switch_window(args...; kwargs...) = _invoke_kernel(:_edg_low_alpha_switch_window, args...; kwargs...)
_edg_max_heat_load_for_links(args...; kwargs...) = _invoke_kernel(:_edg_max_heat_load_for_links, args...; kwargs...)
_edg_maxwellian_heat_rate(args...; kwargs...) = _invoke_kernel(:_edg_maxwellian_heat_rate, args...; kwargs...)
_edg_mean_anomaly_from_true(args...; kwargs...) = _invoke_kernel(:_edg_mean_anomaly_from_true, args...; kwargs...)
_edg_orbit_metrics_from_rv(args...; kwargs...) = _invoke_kernel(:_edg_orbit_metrics_from_rv, args...; kwargs...)
_edg_ordered_feasible_prefix(args...; kwargs...) = _invoke_kernel(:_edg_ordered_feasible_prefix, args...; kwargs...)
_edg_outbound_exit_index(args...; kwargs...) = _invoke_kernel(:_edg_outbound_exit_index, args...; kwargs...)
_edg_padded_heat_load_window(args...; kwargs...) = _invoke_kernel(:_edg_padded_heat_load_window, args...; kwargs...)
_edg_pass_heat_load_for_links(args...; kwargs...) = _invoke_kernel(:_edg_pass_heat_load_for_links, args...; kwargs...)
_edg_predict_mass(args...; kwargs...) = _invoke_kernel(:_edg_predict_mass, args...; kwargs...)
_edg_predict_max_energy_depletion_outcome(args...; kwargs...) = _invoke_kernel(:_edg_predict_max_energy_depletion_outcome, args...; kwargs...)
_edg_predict_targeting_outcome(args...; kwargs...) = _invoke_kernel(:_edg_predict_targeting_outcome, args...; kwargs...)
_edg_prediction_advance(args...; kwargs...) = _invoke_kernel(:_edg_prediction_advance, args...; kwargs...)
_edg_prediction_event(args...; kwargs...) = _invoke_kernel(:_edg_prediction_event, args...; kwargs...)
_edg_prediction_gravity_acceleration(args...; kwargs...) = _invoke_kernel(:_edg_prediction_gravity_acceleration, args...; kwargs...)
_edg_prediction_step(args...; kwargs...) = _invoke_kernel(:_edg_prediction_step, args...; kwargs...)
_edg_prediction_time_grid(args...; kwargs...) = _invoke_kernel(:_edg_prediction_time_grid, args...; kwargs...)
_edg_profile_heat_load(args...; kwargs...) = _invoke_kernel(:_edg_profile_heat_load, args...; kwargs...)
_edg_profile_heat_rates(args...; kwargs...) = _invoke_kernel(:_edg_profile_heat_rates, args...; kwargs...)
_edg_recompute_second_heat_load_switch(args...; kwargs...) = _invoke_kernel(:_edg_recompute_second_heat_load_switch, args...; kwargs...)
_edg_recompute_switches!(args...; kwargs...) = _invoke_kernel(:_edg_recompute_switches!, args...; kwargs...)
_edg_remaining_closed_form_heat_load_trajectory(args...; kwargs...) = _invoke_kernel(:_edg_remaining_closed_form_heat_load_trajectory, args...; kwargs...)
_edg_remaining_closed_form_reference(args...; kwargs...) = _invoke_kernel(:_edg_remaining_closed_form_reference, args...; kwargs...)
_edg_replan_due(args...; kwargs...) = _invoke_kernel(:_edg_replan_due, args...; kwargs...)
_edg_run_target_energy_bracketing!(args...; kwargs...) = _invoke_kernel(:_edg_run_target_energy_bracketing!, args...; kwargs...)
_edg_select_mapped_candidate(args...; kwargs...) = _invoke_kernel(:_edg_select_mapped_candidate, args...; kwargs...)
_edg_set_targeting_fallback!(args...; kwargs...) = _invoke_kernel(:_edg_set_targeting_fallback!, args...; kwargs...)
_edg_solve_heat_load_switches(args...; kwargs...) = _invoke_kernel(:_edg_solve_heat_load_switches, args...; kwargs...)
_edg_solve_targeting_switch(args...; kwargs...) = _invoke_kernel(:_edg_solve_targeting_switch, args...; kwargs...)
_edg_structural_alpha(args...; kwargs...) = _invoke_kernel(:_edg_structural_alpha, args...; kwargs...)
_edg_switch_due(args...; kwargs...) = _invoke_kernel(:_edg_switch_due, args...; kwargs...)
_edg_target_energy_from_apoapsis(args...; kwargs...) = _invoke_kernel(:_edg_target_energy_from_apoapsis, args...; kwargs...)
_edg_target_energy_from_reachable_bracket(args...; kwargs...) = _invoke_kernel(:_edg_target_energy_from_reachable_bracket, args...; kwargs...)
_edg_targeting_aero_acceleration(args...; kwargs...) = _invoke_kernel(:_edg_targeting_aero_acceleration, args...; kwargs...)
_edg_targeting_bracket_outcomes(args...; kwargs...) = _invoke_kernel(:_edg_targeting_bracket_outcomes, args...; kwargs...)
_edg_targeting_candidate_certificate(args...; kwargs...) = _invoke_kernel(:_edg_targeting_candidate_certificate, args...; kwargs...)
_edg_targeting_constrained_alpha(args...; kwargs...) = _invoke_kernel(:_edg_targeting_constrained_alpha, args...; kwargs...)
_edg_targeting_control_period(args...; kwargs...) = _invoke_kernel(:_edg_targeting_control_period, args...; kwargs...)
_edg_targeting_outcome_with_heat_load(args...; kwargs...) = _invoke_kernel(:_edg_targeting_outcome_with_heat_load, args...; kwargs...)
_edg_targeting_prediction_time_grid(args...; kwargs...) = _invoke_kernel(:_edg_targeting_prediction_time_grid, args...; kwargs...)
_edg_targeting_switch_outcomes(args...; kwargs...) = _invoke_kernel(:_edg_targeting_switch_outcomes, args...; kwargs...)
_edg_targeting_switch_solution(args...; kwargs...) = _invoke_kernel(:_edg_targeting_switch_solution, args...; kwargs...)
_edg_total_ref_area(args...; kwargs...) = _invoke_kernel(:_edg_total_ref_area, args...; kwargs...)
_edg_update_heat_load_security!(args...; kwargs...) = _invoke_kernel(:_edg_update_heat_load_security!, args...; kwargs...)
_edg_vacuum_apoapsis_correction(args...; kwargs...) = _invoke_kernel(:_edg_vacuum_apoapsis_correction, args...; kwargs...)
_edg_vacuum_drag_passage_exit(args...; kwargs...) = _invoke_kernel(:_edg_vacuum_drag_passage_exit, args...; kwargs...)
_edg_weighted_aero_coefficients(args...; kwargs...) = _invoke_kernel(:_edg_weighted_aero_coefficients, args...; kwargs...)
_energy_depletion_heat_rate_calc(args...; kwargs...) = _invoke_kernel(:_energy_depletion_heat_rate_calc, args...; kwargs...)
_energy_depletion_heatrate_root_alpha(args...; kwargs...) = _invoke_kernel(:_energy_depletion_heatrate_root_alpha, args...; kwargs...)
_energy_depletion_struct_drag_area(args...; kwargs...) = _invoke_kernel(:_energy_depletion_struct_drag_area, args...; kwargs...)
_energy_depletion_struct_load_root_alpha(args...; kwargs...) = _invoke_kernel(:_energy_depletion_struct_load_root_alpha, args...; kwargs...)
control_decision!(args...; kwargs...) = _invoke_kernel(:control_decision!, args...; kwargs...)
guidance_decision!(args...; kwargs...) = _invoke_kernel(:guidance_decision!, args...; kwargs...)

end # module EDGAlgorithms
