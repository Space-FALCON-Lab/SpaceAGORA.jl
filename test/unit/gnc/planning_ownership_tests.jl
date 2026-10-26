module PlanningOwnershipTests
using Test, Serialization, SpaceAGORA
import PlotlyJS # the report figures come from its package extension
const SM = SpaceAGORA.SimulationModel
const G = SM.GuidanceHooks
const Reports = G.RPOComparisonReports

# Independent inventory of historical qualified report entry points. Keeping the
# defining function objects lets existing imports and serialized references work.
const REPORT_FUNCTIONS = (
    :rpo_flatten_planner_results, :rpo_write_namedtuple_csv, :rpo_group_metric_mean,
    :rpo_comparison_metric_specs, :rpo_comparison_metric_value,
    :rpo_comparison_metric_includes_result, :rpo_comparison_axis_name,
    :rpo_comparison_trace_axis_name, :rpo_comparison_metric_summary_plot,
    :rpo_comparison_station_mesh_trace, :rpo_comparison_station_trace,
    :rpo_comparison_path_family_plot, :rpo_comparison_single_path_plot,
    :rpo_comparison_failed_paths_plot, :rpo_comparison_cost_iteration_xvalues,
    :rpo_comparison_cost_iteration_plot, :rpo_comparison_artifact_slug,
    :rpo_write_failed_path_outputs, :rpo_write_planner_comparison_outputs,
)

@testset "Historical report identities, analysis-owned methods" begin
    for name in REPORT_FUNCTIONS
        f = getproperty(G, name)
        @test getproperty(Reports, name) === f
        @test parentmodule(f) === G
        @test all(m -> m.module === Reports, methods(f))
        io = IOBuffer()
        serialize(io, f)
        seekstart(io)
        @test deserialize(io) === f
    end
    for name in (:rpo_flatten_planner_results, :rpo_comparison_cost_iteration_plot,
                 :rpo_write_planner_comparison_outputs)
        @test getproperty(SM, name) === getproperty(G, name)
    end
    @test !isdefined(Reports, :RPOPSOConfig)
    @test !isdefined(Reports, :HYPRSupport)
    @test SpaceAGORA.HYPRServices.RPO.RPOPSOConfig === G.RPOPSOConfig
    @test SpaceAGORA.HYPRServices.RPO.rpo_retime_samples === G.rpo_retime_samples
    @test SpaceAGORA.HYPRServices.RobotArm.RobotArmPlan === SM.RobotArmPlanning.RobotArmPlan
    @test parentmodule(SM.RPOReferenceGeometry) === SM.NavigationHooks
end

# These function objects stay in GuidanceHooks for existing qualified callers and
# serialized references; only their implementations move to the sibling owner.
const TRACKING_FUNCTIONS = (
    :rpo_lqmpc_reference_preview, :rpo_lqmpc_tracking_fuel_used_pct,
    :rpo_track_retimed_path_lqmpc, :rpo_comparison_segment_clearance,
    :rpo_track_two_body_path, :rpo_execute_two_body_reference,
)

@testset "Historical tracking identities, explicit implementation owner" begin
    tracking_owner = SM.RPOComparisonTracking
    @test parentmodule(tracking_owner) === SM
    @test !isdefined(G, :RPOComparisonTracking)
    @test !isdefined(SpaceAGORA, :RPOComparisonTracking)
    for name in TRACKING_FUNCTIONS
        f = getproperty(G, name)
        @test getproperty(tracking_owner, name) === f
        @test parentmodule(f) === G
        @test length(methods(f)) == 1
        @test all(m -> m.module === tracking_owner, methods(f))
        io = IOBuffer()
        serialize(io, f)
        seekstart(io)
        @test deserialize(io) === f
        # Only the existing exported tracker is bound in SimulationModel. The
        # helper generics were not root API and must not become new root names.
        if name === :rpo_track_retimed_path_lqmpc
            @test getproperty(SM, name) === f
            @test Base.isexported(G, name)
            @test Base.isexported(SM, name)
        else
            @test !isdefined(SM, name)
            @test !Base.isexported(G, name)
        end
        @test !isdefined(SpaceAGORA, name)
    end

    # Types remain defined by their original owners, not aliases to relocated
    # types. Do not construct or serialize the controller's native solver state.
    for name in (:RPOLQMPCTrackingSettings, :RPOPlannerComparisonCase,
                 :RPOPlannerComparisonConfig, :RPOTwoBodyPlant)
        T = getproperty(G, name)
        @test parentmodule(T) === G
        if name === :RPOTwoBodyPlant
            @test !isdefined(SM, name)
        else
            @test getproperty(SM, name) === T
        end
        @test !isdefined(SpaceAGORA, name)
        io = IOBuffer()
        serialize(io, T)
        seekstart(io)
        @test deserialize(io) === T
    end
    # Each dependency must resolve in this module world. A binding accidentally
    # imported from the package's other composition would preserve its spelling
    # while using the wrong function/type identity.
    for (source, names) in (
        (SM.ControlHooks, (:init_rpo_lqmpc, :rpo_lqmpc_control,
            :rpo_allocate_six_axis_thrusters, :rpo_thruster_wrench_body)),
        (G, (:RPOPSOConfig, :RPOLQMPCTrackingSettings, :rpo_pso_config,
            :rpo_reference_from_path, :rpo_init_two_body_plant, :rpo_step_two_body!)),
        (SM.NavigationHooks, (:rpo_clearance_distance_to_station,)),
        (SM.FrameTransforms, (:rtn_dcm_from_inertial,)),
    )
        for name in names
            @test getproperty(tracking_owner, name) === getproperty(source, name)
        end
    end
    @test tracking_owner.HYPRSupport === SM.HYPRSupport
    @test tracking_owner.RPOLQMPCTrackingSettings === G.RPOLQMPCTrackingSettings
    @test parentmodule(SM.RpoLQMPCController) === SM.ControlHooks
    @test SM.RpoLQMPCController === SM.ControlHooks.RpoLQMPCController
    @test parentmodule(SM.SixAxisThrusterModel) === SM.ThrusterModels
    @test SM.ControlHooks.SixAxisThrusterModel === SM.ThrusterModels.SixAxisThrusterModel
    @test tracking_owner.SixAxisThrusterModel === SM.ThrusterModels.SixAxisThrusterModel
    @test !isdefined(SpaceAGORA, :RpoLQMPCController)
    @test !isdefined(SpaceAGORA, :SixAxisThrusterModel)
end

_call_name(x::Symbol) = x
_call_name(x::Expr) = x.head === :. ? _call_name(last(x.args)) : nothing
_call_name(x::QuoteNode) = _call_name(x.value)
_call_name(x) = nothing

_mentions_control_module(x) = x === :ControlHooks || x == "ControlHooks"
_mentions_control_module(x::QuoteNode) = _mentions_control_module(x.value)
_mentions_control_module(x::Expr) = any(_mentions_control_module, x.args)

function _has_reflective_control_lookup(x)
    x === :_control_module && return true
    x isa QuoteNode && return _has_reflective_control_lookup(x.value)
    x isa Expr || return false
    if x.head === :call && !isempty(x.args) &&
       _call_name(first(x.args)) in (:getfield, :getproperty, :eval) &&
       any(_mentions_control_module, x.args[2:end])
        return true
    end
    return any(_has_reflective_control_lookup, x.args)
end

@testset "Comparison tracking does not discover ControlHooks at runtime" begin
    root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    comparison_dir = joinpath(root, "src", "models", "gnc", "planning", "rpo", "comparison")
    paths = String[]
    for (dir, _, files) in walkdir(comparison_dir)
        append!(paths, [joinpath(dir, name) for name in files if endswith(name, ".jl")])
    end
    push!(paths, joinpath(root, "src", "interfaces", "analysis", "planning", "rpo_comparison_tracking.jl"))
    for path in sort(paths)
        @testset "$(relpath(path, root))" begin
            @test !_has_reflective_control_lookup(Meta.parseall(read(path, String)))
        end
    end
    # Put the previous dependency pattern back into an actual controller call
    # in memory. Exactly the altered comparison source must be rejected; no
    # repository source is changed and the legitimate imports still pass.
    tracking_path = joinpath(root, "src", "interfaces", "analysis", "planning", "rpo_comparison_tracking.jl")
    original = read(tracking_path, String)
    original_call = "ctrl = init_rpo_lqmpc("
    @test occursin(original_call, original)
    for lookup in ("_control_module()",
                   "getfield(parentmodule(@__MODULE__), :ControlHooks)",
                   "Base.getproperty(parentmodule(@__MODULE__), Symbol(\"ControlHooks\"))")
        changed = replace(original, original_call => "ctrl = ($lookup).init_rpo_lqmpc("; count=1)
        @test changed != original
        violations = filter(paths) do path
            source = path == tracking_path ? changed : read(path, String)
            _has_reflective_control_lookup(Meta.parseall(source))
        end
        @test violations == [tracking_path]
    end
    @test !_has_reflective_control_lookup(Meta.parseall("import ..ControlHooks: init_rpo_lqmpc"))
    # Guidance as a whole is checked below.
end

_is_parent_call(x) = x isa Expr && x.head === :call && _call_name(first(x.args)) === :parentmodule

# The names a source binds to its parent module or an ancestor
# (`const _PARENT = parentmodule(@__MODULE__)`), or to a name already in `found`
# (`const P = _PARENT`). One pass; `_all_parent_aliases` repeats it.
function _parent_aliases(x, found=Set{Symbol}())
    if x isa Expr
        x.head === :(=) && x.args[1] isa Symbol && (_is_parent_call(x.args[2]) || x.args[2] in found) &&
            push!(found, x.args[1])
        foreach(arg -> _parent_aliases(arg, found), x.args)
    end
    return found
end

# The parent-module aliases of a set of sources, with simple assignments from an
# alias followed to a fixed point, across files: an alias one file defines and
# another renames counts.
function _all_parent_aliases(exprs)
    found = Set{Symbol}()
    while true
        before = length(found)
        foreach(ex -> _parent_aliases(ex, found), exprs)
        length(found) == before && return found
    end
end

# Each name a source looks up at run time through its parent module, by name or
# through a call to `parentmodule`: `getfield` and `getproperty` calls and dot
# access (`_PARENT.Name`), in source order.
function _parent_module_lookups(x, aliases, found=Any[])
    through(owner) = owner in aliases || _is_parent_call(owner)
    if x isa Expr
        if x.head === :call && length(x.args) >= 3 && _call_name(first(x.args)) in (:getfield, :getproperty) &&
                through(x.args[2])
            push!(found, x.args[3] isa QuoteNode ? x.args[3].value : x.args[3])
        elseif x.head === :. && length(x.args) == 2 && x.args[2] isa QuoteNode && through(x.args[1])
            push!(found, x.args[2].value)
        end
        foreach(arg -> _parent_module_lookups(arg, aliases, found), x.args)
    end
    return found
end

@testset "Guidance does not discover ControlHooks at runtime" begin
    root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    guidance_dir = joinpath(root, "src", "models", "gnc", "guidance")
    paths = String[]
    for (dir, _, files) in walkdir(guidance_dir)
        append!(paths, [joinpath(dir, name) for name in files if endswith(name, ".jl")])
    end
    sort!(paths)
    sources = Dict(path => read(path, String) for path in paths)
    for path in paths
        @testset "$(relpath(path, root))" begin
            @test !_has_reflective_control_lookup(Meta.parseall(sources[path]))
        end
    end
    @test !isdefined(SpaceAGORA.SimulationModel.GuidanceHooks, :_control_module)
    # Guidance looks nothing up through its parent module at run time: typed EDG
    # guidance imports EDGAlgorithms, which loads before GuidanceHooks. Any such
    # lookup is new and needs review.
    allowed = Dict{String, Vector{Any}}()
    function lookups(sources)
        parsed = Dict(path => Meta.parseall(source) for (path, source) in sources)
        aliases = _all_parent_aliases(values(parsed))
        return Dict(relpath(path, root) => found for (path, ex) in parsed
            for found in (_parent_module_lookups(ex, aliases),) if !isempty(found))
    end
    @test :_CORE in _all_parent_aliases(Meta.parseall.(values(sources)))
    @test lookups(sources) == allowed
    # A ControlHooks lookup put back into an actual guidance source, in memory, is
    # caught; no repository file changes. So is the EDGAlgorithms lookup guidance
    # no longer needs.
    edg_path = joinpath(guidance_dir, "aerobraking", "energy_depletion_lifecycle.jl")
    for lookup in ("getfield(_PARENT, :ControlHooks)", "getproperty(parentmodule(@__MODULE__), :ControlHooks)",
                   "_PARENT.ControlHooks", "_CORE.ControlHooks", "getfield(_PARENT, :EDGAlgorithms)")
        changed = sources[edg_path] * "\nconst _PARENT = parentmodule(@__MODULE__)\n_late_control() = $lookup\n"
        @test lookups(merge(sources, Dict(edg_path => changed))) != allowed
        startswith(lookup, "get") && occursin("ControlHooks", lookup) &&
            @test _has_reflective_control_lookup(Meta.parseall(changed))
    end
    # So is one through a renamed alias of the parent: in the same file, through a
    # second renaming, or with the alias defined in one guidance file and used in
    # another.
    hooks_path = joinpath(guidance_dir, "hooks.jl")
    for (alias_path, alias_lines, lookup) in (
            (edg_path, "const P = _CORE", "P.ControlHooks"),
            (edg_path, "const P = _CORE\nconst Q = P", "Q.ControlHooks"),
            (edg_path, "const P = _CORE", "getfield(P, :ControlHooks)"),
            (hooks_path, "const P = _CORE", "P.ControlHooks"))
        changed = merge(sources, Dict(alias_path => sources[alias_path] * "\n$alias_lines\n"))
        changed[edg_path] = changed[edg_path] * "\n_late_control() = $lookup\n"
        @test lookups(changed) != allowed
    end
end

@testset "Reporting preserves result selection and file content" begin
    success = (planner=:rrt_connect, case_label="one,two", success=true,
        planner_compute_time=0.5, fuel_used_pct=1.0, control_effort_total=2.0,
        thrust_saturation_fraction=0.0, min_clearance=3.0, final_pos_error=0.1,
        actual_travel_duration=4.0, keepout_violations=0)
    failure = merge(success, (case_label="failed", success=false, fuel_used_pct=99.0))
    batch = (planner_types=[:rrt_connect],
        results_by_planner=Dict(:rrt_connect=>[success, failure]),
        plans_by_planner=Dict(:rrt_connect=>[(case=(label="one,two",),
            cost_history=[5.0, NaN, 3.0], planner_iteration_count=4)]),)
    @test G.rpo_flatten_planner_results(batch) == [success, failure]
    @test G.rpo_group_metric_mean([(v=1.0,), (v=NaN,), (v=3.0,)], :v) == 2.0
    @test isnan(G.rpo_group_metric_mean([(v=NaN,)], :v))
    @test G.rpo_comparison_cost_iteration_xvalues((iterations=7,), [3.0]) == [7]
    @test G.rpo_comparison_cost_iteration_xvalues((iterations=7,), [5.0, 3.0]) == [1, 2]
    @test G.rpo_comparison_cost_iteration_xvalues((;), Float64[]) == Int[]
    plot = G.rpo_comparison_metric_summary_plot(batch)
    boxes = filter(t -> get(t.fields, :type, "") == "box", plot.data)
    @test length(boxes) == 9
    @test boxes[1].fields[:y] == [100.0, 0.0]
    @test boxes[3].fields[:y] == [1.0] # failed runs do not enter the fuel aggregate
    costs = G.rpo_comparison_cost_iteration_plot(batch; planner=:rrtconnect)
    @test length(costs.data) == 1
    @test costs.data[1].fields[:x] == [1, 2]
    @test costs.data[1].fields[:y] == [5.0, 3.0]
    mktempdir() do dir
        output = G.rpo_write_planner_comparison_outputs(batch; output_dir=dir,
            write_plotly_outputs=false, write_failed_path_outputs=false)
        @test basename(output.csv) == "rpo_planner_comparison_results.csv"
        lines = readlines(output.csv)
        @test lines[1] == join(string.(keys(success)), ",")
        @test lines[2] == "rrt_connect,one;two,true,0.5,1.0,2.0,0.0,3.0,0.1,4.0,0"
        @test lines[3] == "rrt_connect,failed,false,0.5,99.0,2.0,0.0,3.0,0.1,4.0,0"
        @test readdir(dir) == [basename(output.csv)]
        @test output.metrics_plot === output.path_plot === nothing
        @test isempty(output.method_path_plots) && isempty(output.cost_iteration_plots)
        @test isempty(output.failed_path_plots)
        @test !output.plotly_outputs && !output.failed_path_outputs
    end
end

@testset "Station geometry preserves coordinates and mesh connectivity" begin
    points = Float64[1 4 7 10 13 16 19; 2 5 8 11 14 17 20; 3 6 9 12 15 18 21]
    batch = (geometry=(station=(points_body=points,),),)
    cloud = G.rpo_comparison_station_trace(batch; max_points=2)
    @test cloud.fields[:type] == "scatter3d"
    @test cloud.fields[:mode] == "markers"
    @test cloud.fields[:x] == [1.0, 13.0]
    @test cloud.fields[:y] == [2.0, 14.0]
    @test cloud.fields[:z] == [3.0, 15.0]
    # An explicitly absent mesh uses the same point-cloud representation.
    no_mesh = G.rpo_comparison_station_trace(merge(batch, (station_triangles=nothing,)))
    @test no_mesh.fields[:x] == points[1, :]
    @test no_mesh.fields[:y] == points[2, :]
    @test no_mesh.fields[:z] == points[3, :]
    triangles = Float64[0 1 0 2 3 2; 0 0 1 0 0 1; 0 0 0 1 1 1]
    mesh = G.rpo_comparison_station_trace(merge(batch, (station_triangles=triangles,)))
    @test mesh.fields[:type] == "mesh3d"
    @test mesh.fields[:x] == [0.0, 1.0, 0.0, 2.0, 3.0, 2.0]
    @test mesh.fields[:y] == [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
    @test mesh.fields[:z] == [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    # Plotly indices are zero-based: the two faces must remain disconnected.
    @test mesh.fields[:i] == [0, 3]
    @test mesh.fields[:j] == [1, 4]
    @test mesh.fields[:k] == [2, 5]
end

@testset "Path reports distinguish reference, tracked, raw and failed cases" begin
    reference = Float64[10 11 12; 20 21 22; 30 31 32]
    actual = Float64[40 41 42; 50 51 52; 60 61 62]
    raw = Float64[70 71; 80 81; 90 91]
    # Distinct velocity rows catch an accidental use of the whole state history.
    tracking = (r_ref_rtn=reference, x_hist=vcat(actual, fill(-100.0, 3, 3)))
    selected = (case=(label="selected", start_rtn=[1.0, 2.0, 3.0],
                      goal_rtn=[4.0, 5.0, 6.0]), tracking=tracking, raw_path=raw)
    other = (case=(label="other", start_rtn=[-1.0, -2.0, -3.0],
                   goal_rtn=[-4.0, -5.0, -6.0]),
             tracking=(r_ref_rtn=-reference, x_hist=-tracking.x_hist))
    batch = (geometry=(station=(points_body=zeros(3, 1),),),
             planner_types=[:rrt_connect, :chomp],
             plans_by_planner=Dict(:rrt_connect=>[other, selected], :chomp=>[other]))
    single = G.rpo_comparison_single_path_plot(batch, selected; planner=:rrtconnect)
    @test [t.fields[:name] for t in single.data] ==
        ["station point cloud", "desired", "tracked", "raw planner path", "start", "goal"]
    for (idx, coordinates) in ((2, reference), (3, actual), (4, raw))
        @test single.data[idx].fields[:x] == coordinates[1, :]
        @test single.data[idx].fields[:y] == coordinates[2, :]
        @test single.data[idx].fields[:z] == coordinates[3, :]
    end
    @test single.data[4].fields[:mode] == "lines+markers"
    @test single.data[5].fields[:x] == [1.0]
    @test single.data[5].fields[:y] == [2.0]
    @test single.data[5].fields[:z] == [3.0]
    @test single.data[6].fields[:x] == [4.0]
    @test single.data[6].fields[:y] == [5.0]
    @test single.data[6].fields[:z] == [6.0]
    no_raw = G.rpo_comparison_single_path_plot(batch, other; planner=:rrt_connect)
    @test [t.fields[:name] for t in no_raw.data] ==
        ["station point cloud", "desired", "tracked", "start", "goal"]
    @test no_raw.data[2].fields[:x] == -reference[1, :]

    failed = G.rpo_comparison_failed_paths_plot(batch; planner=:rrtconnect, failed_indices=[2])
    @test [t.fields[:name] for t in failed.data] ==
        ["station point cloud", "case selected desired", "case selected tracked",
         "case selected raw planner path", "case selected start", "case selected goal"]
    for (idx, coordinates) in ((2, reference), (3, actual), (4, raw))
        @test failed.data[idx].fields[:x] == coordinates[1, :]
        @test failed.data[idx].fields[:y] == coordinates[2, :]
        @test failed.data[idx].fields[:z] == coordinates[3, :]
    end
    @test failed.data[5].fields[:x] == [1.0]
    @test failed.data[6].fields[:z] == [6.0]
    failed_no_raw = G.rpo_comparison_failed_paths_plot(batch; planner=:rrt_connect, failed_indices=[1])
    @test [t.fields[:name] for t in failed_no_raw.data] ==
        ["station point cloud", "case other desired", "case other tracked", "case other start", "case other goal"]
    @test failed_no_raw.data[3].fields[:y] == -actual[2, :]
    @test length(G.rpo_comparison_failed_paths_plot(batch; planner=:rrt_connect, failed_indices=Int[]).data) == 1

    family = G.rpo_comparison_path_family_plot(batch; planner=:rrtconnect)
    @test length(family.data) == 9
    @test all(t -> !occursin("CHOMP", t.fields[:name]), family.data)
    @test family.data[2].fields[:x] == -reference[1, :]
    @test family.data[6].fields[:x] == reference[1, :]
    @test family.data[7].fields[:z] == actual[3, :]
    @test family.data[8].fields[:y] == [2.0]
    @test family.data[9].fields[:z] == [6.0]
    all_planners = G.rpo_comparison_path_family_plot(batch)
    @test length(all_planners.data) == 13
    @test all_planners.data[10].fields[:name] == "CHOMP case other desired"
end
end # module PlanningOwnershipTests
