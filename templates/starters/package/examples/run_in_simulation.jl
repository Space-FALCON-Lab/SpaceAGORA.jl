# The planner in a short simulation of SpaceAGORA's RPO pilot: a chaser 3 m from
# the station plans to 5 m. The workflow runs it after the tests.
using SpaceAGORA, StarterName

configuration = make_rpo_configuration(planner=StraightLinePlanner(), mission_time_s=1.0)
solution = run_simulation(configuration; return_solution=true)
report = rpo_run_report(solution)
installed = only(filter(record -> record.event === :installed, only(report.planners).records))
println("run_in_simulation: t = ", last(solution.t), " s, reference samples: ", length(installed.reference.t_ref_s),
    ", validation: ", installed.validation.reason)
