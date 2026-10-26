# A completed two-planner comparison, small enough to read: RRT-Connect with a
# successful and a failed case, CHOMP with one case, and a station point cloud.
# Included by rpo_comparison_report_data_tests.jl and by its child processes.
ok(planner, label, fuel) = (planner=planner, case_label=label, success=true, planner_compute_time=0.5,
    fuel_used_pct=fuel, control_effort_total=2.0, thrust_saturation_fraction=0.0, min_clearance=3.0,
    final_pos_error=0.1, actual_travel_duration=4.0, keepout_violations=0)
failed(planner, label) = merge(ok(planner, label, 99.0), (success=false, keepout_violations=2))
path(label; raw=nothing, costs=nothing) = begin
    ref = Float64[0 1 2; 0 0 1; 0 0 0]
    nt = (case=(label=label, start_rtn=[0.0, 0.0, 0.0], goal_rtn=[2.0, 1.0, 0.0]),
          tracking=(r_ref_rtn=ref, x_hist=vcat(ref .+ 0.5, zeros(3, 3))), planner_iteration_count=6)
    raw === nothing || (nt = merge(nt, (raw_path=raw,)))
    costs === nothing || (nt = merge(nt, (cost_history=costs,)))
    nt
end
batch = (planner_types=[:rrt_connect, :chomp],
    results_by_planner=Dict(:rrt_connect => [ok(:rrt_connect, "a", 1.0), failed(:rrt_connect, "b")],
                            :chomp => [ok(:chomp, "a", 3.0)]),
    plans_by_planner=Dict(:rrt_connect => [path("a"; costs=[5.0, NaN, 3.0]), path("b"; raw=[0 2; 0 1; 0 0])],
                          :chomp => [path("a"; costs=[4.0])]),
    geometry=(station=(points_body=Float64[1 4 7 10; 2 5 8 11; 3 6 9 12],),))
