# Historical function identities remain in GuidanceHooks for existing exports,
# qualified clients and serialized function references. Analysis owns the methods.
function rpo_flatten_planner_results end
function rpo_write_namedtuple_csv end
function rpo_group_metric_mean end
function rpo_comparison_metric_specs end
function rpo_comparison_metric_value end
function rpo_comparison_axis_name end
function rpo_comparison_trace_axis_name end
function rpo_comparison_metric_summary_plot end
function rpo_comparison_station_mesh_trace end
function rpo_comparison_station_trace end
function rpo_comparison_path_family_plot end
function rpo_comparison_single_path_plot end
function rpo_comparison_failed_paths_plot end
function rpo_comparison_cost_iteration_xvalues end
function rpo_comparison_cost_iteration_plot end
function rpo_comparison_artifact_slug end
function rpo_write_failed_path_outputs end
function rpo_write_planner_comparison_outputs end
function rpo_comparison_metric_includes_result end

include(joinpath(@__DIR__, "..", "..", "..", "..", "..", "interfaces", "analysis", "planning", "rpo_comparison_reports.jl"))
# Preserve historical qualified access without a second PlotlyJS import owner.
using .RPOComparisonReports: PlotlyJS
