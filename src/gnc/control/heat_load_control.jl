# Compatibility bindings; typed EDG calculations have one owner.
using ..EDGAlgorithms: _edg_heat_load_scale_height,
    _edg_total_ref_area,
    _edg_predict_mass,
    _edg_max_heat_load_for_links,
    _edg_weighted_aero_coefficients,
    _edg_heat_load_coefficients,
    _edg_eccentric_anomaly_from_true,
    _edg_mean_anomaly_from_true,
    _edg_drag_passage_duration,
    _edg_prediction_time_grid,
    _edg_closed_form_heat_load_trajectory,
    _edg_integrated_heat_load_trajectory,
    _edg_heat_load_lambdas,
    _edg_heat_load_alpha_profile,
    _edg_heat_load_track_env,
    _edg_constrained_heat_load_alpha_profile,
    _edg_profile_heat_rates,
    _edg_integrate_series,
    _edg_profile_heat_load,
    _edg_first_low_alpha_interval_indices,
    _edg_first_two_switch_alpha_profile,
    _edg_low_alpha_switch_window,
    _edg_balanced_tpbvp_heat_load_window,
    _edg_padded_heat_load_window,
    _edg_heat_load_profile_for_k,
    _edg_solve_heat_load_switches
using ..EDGServices: _edg_sample_prediction_atmosphere

@inline function _edg_heat_load_low_drag_active(model, t::Float64, i::Int)::Bool
    return EDGAlgorithms._edg_heat_load_low_drag_active(model.config, model.state, t, i)
end
