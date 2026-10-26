module ConstellationResultsTests
# `run_constellation_ensemble(args; merge_results=true)`: the members table built
# by the long builder, identity, order across routes, unchanged trajectories and
# sample values, members with different columns, failures, files, refusals,
# run-wide columns beside per-spacecraft ones, and a resumed ensemble.
using Test
using LinearAlgebra
using StaticArrays
using Logging
using DataFrames
using Arrow
using CSV
using TOML
using Base64
import OrdinaryDiffEq
using SpaceAGORA
using SpaceAGORA.SimulationModel
import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

const SM = SpaceAGORA.SimulationModel
const SC = SpaceAGORA.SimulationCampaigns
const SE = SpaceAGORA.SimulationEngine
const PLANET = make_no_gram_planet(:earth)
const T0 = InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
const Z3 = SVector(0.0, 0.0, 0.0)
const IDENTITY = ["spacecraft", "spacecraft_id", "plane", "slot", "time"]

quiet(f) = with_logger(f, ConsoleLogger(stderr, Logging.Error))

# A model that fails when evaluated, and one that takes time on every evaluation.
struct Fails <: SpaceAGORA.AbstractForceTorqueModel end
SpaceAGORA.wrench(::Fails, x::StateSample, env::EnvironmentSample, t::Float64) = error("this member's model failed")
struct Slow <: SpaceAGORA.AbstractForceTorqueModel end
SpaceAGORA.wrench(::Slow, x::StateSample, env::EnvironmentSample, t::Float64) = (Libc.systemsleep(2.0e-4); (Z3, Z3))
# A model bound to a run position by its own field, which an ensemble refuses.
struct Bound <: SpaceAGORA.AbstractForceTorqueModel
    spacecraft_index::Int
end
SpaceAGORA.wrench(::Bound, x::StateSample, env::EnvironmentSample, t::Float64) = (Z3, Z3)

function three_body(index; id=index, ν=30.0 * index)
    return make_three_body_spacecraft(bus_dims=(1.0, 1.0, 1.0), panel_dims=(0.01, 1.0, 0.5), bus_mass=200.0,
        panel_mass_each=5.0, panel_offset_y=1.0,
        ic=InitialCondition(ra=7.0e6, rp=7.0e6, i=53.0, ω=0.0, Ω=0.0, ν=ν), prop_mass=0.0, id=id)
end

function one_link(index; id=index)
    root = SM.Link(root=true, m=300.0, ref_area=2.0)
    return SM.SpacecraftModel(SM.Joint[], [root], root, true, 300.0, 0.0, root.inertia, 0, 0,
        SM.InitialCondition(a=PLANET.Rp_e + 6.0e5, e=0.0, i=53.0, Ω=0.0, ω=0.0, ν=40.0 * index), id)
end

# A Walker 6/3/1 constellation with ids 10k + 3, and one spacecraft without a slot.
function fleet()
    walker = walker_constellation((plane, slot, index) -> three_body(index);
        total=6, planes=3, phasing=1, altitude=6.0e5, planet=PLANET, i=53)
    foreach(((k, sc),) -> sc.id = 10 * k + 3, enumerate(walker))
    return [walker; three_body(7; id=99)]
end

function configuration(spacecraft; effectors=(SM.InverseSquaredJ2GravityModel(),), dir=mktempdir(), results=true,
                       solver_mode=:tsit5, mission_time=300.0, settings...)
    base = make_example_config(planet=PLANET, spacecraft=spacecraft[1], mission_time=mission_time, initial_time=T0,
        dynamic_effectors=effectors, density_model=SM.NoAtmosphereModel(), ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=false, keplerian=true, EI_km=120.0, verbose=false, results=results, results_directory=dir,
        solver_config=SM.SolverConfig(solver_mode=solver_mode))
    simulation_settings = isempty(settings) ? base.simulation_settings :
        SM.SimConfig._with_fields(base.simulation_settings; settings...)
    return SM.SimConfig._with_configuration(base; simulation_settings=simulation_settings,
        dynamics_model=SM.DynamicsModel(spacecraft, effectors))
end

# The configuration member `k` runs, as the ensemble builds it.
member(args, k) = SC._ensemble_member_configuration(SE._resolve_effector_scopes(args), args.dynamics_model.spacecraft[k],
    k, length(args.dynamics_model.spacecraft))

# What member `k`'s rows in the members table must be, from its own run with the
# default (wide) layout, rearranged by `long_table` rather than by the long
# builder: its spacecraft rows with `spacecraft` set to `k`, and its run-wide
# columns beside them.
function expected_rows(args, k; kwargs...)
    own = quiet(() -> run_simulation(member(args, k); return_results=true, kwargs...))
    rows = long_table(own)
    rows[!, "spacecraft"] = fill(Int32(k), nrow(rows))
    for name in names(own.table)
        (name == "time" || startswith(name, "sc1_")) || (rows[!, name] = own.table[!, name])
    end
    return rows
end

same_columns(a::DataFrame, b::DataFrame) =
    names(a) == names(b) && all(name -> isequal(a[!, name], b[!, name]), names(a))

@testset "summary status cannot be overwritten by saved fields" begin
    args = configuration([three_body(1)]; results=false)
    for name in ("success", "error", "elapsed_s", "rows", "member_tag", "results_subdirectory"), per_satellite in (false, true)
        field = SM.SaveField(:custom, (u, t, i) -> nothing; per_satellite, column_prefix=name)
        value = per_satellite ? [999.0] : 999.0
        tables = SM.IOOutputs._build_long_results_tables([0.0, 10.0],
            [Dict(:custom => value), Dict(:custom => value)], [field], args)
        outcome = (value=nothing, spacecraft=tables.spacecraft, run=tables.run)
        sample = SC.MonteCarloSampleResult(index=1, seed=1, success=true, elapsed_s=0.5, value=outcome)
        campaign = SC.MonteCarloResult([sample], 0.5, 1)
        @test_throws ArgumentError SC._constellation_ensemble_result(args, [args], campaign)
    end
end

@testset "merged field ownership follows constellation positions" begin
    for checkpoints in (false, true)
        args = configuration([three_body(1), one_link(2)]; results_layout=:long,
            checkpoint_enabled=checkpoints, checkpoint_interval_s=100.0)
        gap = SM.SaveField(:gap, (u, t, integrator) -> [missing]; per_satellite=true)
        ragged = SM.SaveField(:ragged, (u, t, integrator) ->
            [u.sc[1].mass < 250.0 ? Float64[] : [t]]; per_satellite=true)
        clock = SM.SaveField(:clock, (u, t, integrator) -> t)
        result = quiet(() -> run_constellation_ensemble(args; merge_results=true, return_results=true, fail_fast=true,
            save_fields=vcat(SM.default_save_fields(member(args, 1)), [gap, ragged, clock])))
        saved = SM.IOOutputs._read_long_results_table(joinpath(args.simulation_settings.results_directory,
            "constellation_members.feather"))
        disk = read_constellation_results(args.simulation_settings.results_directory; overwrite=true)
        disk_saved = SM.IOOutputs._read_long_results_table(disk.files[1])
        for table in (result.members, saved, disk.members, disk_saved)
            wide = try SM.SceneVisualization._wide_view(table, nothing, 2) catch; DataFrame() end
            @test "sc1_gap" in names(wide) && "sc2_gap" in names(wide)
            @test !("sc1_ragged_1" in names(wide)) && "sc2_ragged_1" in names(wide)
            @test "sc1_clock" in names(wide) && "sc2_clock" in names(wide)
        end
        for sample in result.campaign.samples
            own = sample.value
            wide = SM.SceneVisualization._wide_view(own.table, run_table(own), 1)
            @test "sc1_gap" in names(wide) && all(ismissing, wide.sc1_gap)
        end
    end
end

@testset "disk validation preserves both existing merged files" begin
    for layout in (:wide, :long), per_satellite in (false, true), reserved in ("read_status", "success", "results_subdirectory")
        args = configuration([one_link(1)]; mission_time=2.0, results_layout=layout)
        field = SM.SaveField(:custom, (u, t, i) -> per_satellite ? [999.0] : 999.0;
            per_satellite, column_prefix=reserved)
        quiet(() -> run_constellation_ensemble(args; save_fields=[field], fail_fast=true))
        directory = args.simulation_settings.results_directory
        members_path = joinpath(directory, "constellation_members.feather")
        summary_path = joinpath(directory, "constellation_summary.csv")
        write(members_path, "previous merged members")
        CSV.write(summary_path, DataFrame(spacecraft=[1], spacecraft_id=[1], plane=[missing], slot=[missing],
            member_tag=["sat_1_id_1"], success=[true], error=[missing], elapsed_s=[0.0]))
        before = [(read(path), stat(path).mtime) for path in (members_path, summary_path)]
        @test_throws ArgumentError read_constellation_results(directory; overwrite=true)
        @test [(read(path), stat(path).mtime) for path in (members_path, summary_path)] == before
    end
end

@testset "disk reload follows a recovery member's actual output folder" begin
    for layout in (:wide, :long)
        directory = mktempdir()
        first_leg = configuration([one_link(1), one_link(2)]; dir=directory, mission_time=2.0,
            results_layout=layout, save_visualization_scene=true, checkpoint_enabled=true, checkpoint_interval_s=1.0)
        quiet(() -> run_constellation_ensemble(first_leg; merge_results=true, fail_fast=true))
        for k in 1:2
            rm(SE._checkpoint_paths(member(first_leg, k)).rows)
        end
        resume = configuration([one_link(1), one_link(2)]; dir=directory, mission_time=4.0,
            results_layout=layout, save_visualization_scene=true, checkpoint_enabled=true, checkpoint_interval_s=1.0, resume_from_checkpoint=true)
        live = quiet(() -> run_constellation_ensemble(resume; merge_results=true, fail_fast=true))
        @test live.constellation.results_subdirectory == fill("resumed_at_2.0s", 2)
        disk = read_constellation_results(directory; overwrite=true)
        @test same_columns(disk.members, live.members)
        @test first(disk.members.time) == 2.0 && last(disk.members.time) == 4.0
        @test all(==("read"), disk.constellation.read_status)
        samples, _ = SM.SceneVisualization.discover_ensemble_samples(directory)
        @test length(samples) == 2
        @test all(sample -> basename(sample.directory) == "resumed_at_2.0s", samples)
        @test [sample.label for sample in samples] == ["id 1", "id 2"]
        scenes = [SM.SceneVisualization.read_visualization_scene(joinpath(sample.directory, "simulation_results_scene.json")) for sample in samples]
        tables = [() -> SM.SceneVisualization._saved_results_table(scenes[i], joinpath(sample.directory, scenes[i].results_feather))
            for (i, sample) in enumerate(samples)]
        frames, _, _ = SM.SceneVisualization.build_ensemble_frames(samples, scenes, tables)
        @test last(reinterpret(Float64, base64decode(frames["t_s"]))) == 4.0
        # Older summaries cannot choose between an original and a recovery bundle.
        summary_path = joinpath(directory, "constellation_summary.csv")
        CSV.write(summary_path, select(disk.constellation, Not(:results_subdirectory)))
        before = read(joinpath(directory, "constellation_members.feather"))
        @test_throws "no recorded results_subdirectory" read_constellation_results(directory; overwrite=true)
        @test read(joinpath(directory, "constellation_members.feather")) == before
        @test_throws "no recorded results_subdirectory" SM.SceneVisualization.discover_ensemble_samples(directory)
    end
end

@testset "merged constellation results" begin
    @testset "the default ensemble is unchanged" begin
        dir = mktempdir()
        result = quiet(() -> run_constellation_ensemble(configuration(fleet()[1:2]; dir=dir)))
        @test result isa MonteCarloResult
        @test all(sample -> sample.value === nothing, result.samples)
        @test sort(readdir(dir)) == ["sat_1_id_13", "sat_2_id_23"]
    end

    @testset "the members table holds each member's own rows" begin
        args = configuration(fleet(); save_visualization_scene=true)
        res = quiet(() -> run_constellation_ensemble(args; merge_results=true))
        @test res isa ConstellationEnsembleResult
        @test names(res.members)[1:5] == IDENTITY
        @test eltype(res.members.spacecraft) == Int32 && eltype(res.members.plane) == Union{Missing, Int16}
        # A run-wide field (the Sun direction the scene saves) sits beside the
        # per-spacecraft ones, and no name keeps a `sc<i>_` prefix.
        @test "sun_dir_1" in names(res.members) && "pos_1" in names(res.members)
        @test !any(name -> occursin(r"^sc\d+_", name), names(res.members))
        for k in 1:7
            rows = res.members[res.members.spacecraft .== k, :]
            @test same_columns(rows, expected_rows(args, k))
        end
        # Rows are in position order, then saved time.
        @test issorted(collect(zip(res.members.spacecraft, res.members.time)))
    end

    @testset "identity" begin
        args = configuration(fleet())
        res = quiet(() -> run_constellation_ensemble(args; merge_results=true))
        first_rows = res.members[res.members.time .== 0.0, :]
        @test first_rows.spacecraft == Int32.(1:7)
        @test first_rows.spacecraft_id == [13, 23, 33, 43, 53, 63, 99]
        @test isequal(first_rows.plane, Union{Missing, Int16}[1, 1, 2, 2, 3, 3, missing])
        @test isequal(first_rows.slot, Union{Missing, Int16}[1, 2, 1, 2, 1, 2, missing])
        @test isequal(res.constellation.plane, first_rows.plane) && res.constellation.spacecraft_id == first_rows.spacecraft_id
        # A repeated id (deprecated): the members stay apart by position.
        repeated = [three_body(1; id=5), three_body(2; id=5)]
        same_id = quiet(() -> run_constellation_ensemble(configuration(repeated); merge_results=true))
        @test same_id.members.spacecraft_id == fill(5, nrow(same_id.members))
        @test sort(unique(same_id.members.spacecraft)) == Int32[1, 2]
        @test same_id.constellation.member_tag == ["sat_1_id_5", "sat_2_id_5"]
    end

    @testset "the tables do not depend on the route" begin
        # Member 1 runs a slow model, so with worker threads it finishes last.
        spacecraft = [three_body(k) for k in 1:4]
        effectors = (SM.InverseSquaredJ2GravityModel(), scoped(Slow(), BoundSpacecraft((1,))))
        run(; kwargs...) = quiet(() -> run_constellation_ensemble(configuration(spacecraft; effectors=effectors);
            merge_results=true, kwargs...))
        serial = run()
        others = Any[run(threads=:auto, route_state=SC.OuterRouteState())]
        Threads.nthreads() >= 2 && push!(others, run(threads=2))
        for other in others
            @test isequal(other.members, serial.members)
            @test isequal(select(other.constellation, Not(:elapsed_s)), select(serial.constellation, Not(:elapsed_s)))
        end
        if Threads.nthreads() >= 2
            threaded = last(others)
            finished = sortperm([sample.finished_ns for sample in threaded.campaign.samples])
            @test last(finished) == 1                       # it finished last, and still comes first
            @test threaded.members.spacecraft[1] == 1
        end
    end

    @testset "trajectories and sample values are unchanged" begin
        for solver_mode in (:tsit5, :rodas5p)
            args = configuration([three_body(k) for k in 1:3]; results=false, solver_mode=solver_mode)
            plain = quiet(() -> run_constellation_ensemble(args; return_solution=true))
            merged = quiet(() -> run_constellation_ensemble(args; return_solution=true, merge_results=true))
            for (a, b) in zip(plain.samples, merged.campaign.samples)
                @test isequal(a.value.t, b.value.t) && isequal(collect.(a.value.u), collect.(b.value.u))
            end
            @test nrow(merged.members) > 0 && isempty(merged.files)
        end
        args = configuration([three_body(k) for k in 1:2]; results=false)
        none = quiet(() -> run_constellation_ensemble(args; merge_results=true))
        @test all(sample -> sample.value === nothing, none.campaign.samples)
        solutions = quiet(() -> run_constellation_ensemble(args; merge_results=true, return_solution=true))
        @test all(sample -> sample.value isa OrdinaryDiffEq.ODESolution, solutions.campaign.samples)
        results = quiet(() -> run_constellation_ensemble(args; merge_results=true, return_results=true))
        @test all(sample -> sample.value isa SimulationResults, results.campaign.samples)
        # The value is run_simulation's own: its wide table, as without the keyword.
        @test all(sample -> "sc1_pos_1" in names(sample.value.table), results.campaign.samples)
        @test length(results.campaign.successful) == 2 && isempty(results.campaign.failed)
    end

    @testset "returned results in the long layout are left as they are" begin
        for results in (false, true)
            spacecraft = [three_body(k) for k in 1:3]
            setup() = configuration(spacecraft; results=results, results_layout=:long, save_visualization_scene=true)
            plain = quiet(() -> run_constellation_ensemble(setup(); return_results=true))
            merged = quiet(() -> run_constellation_ensemble(setup(); return_results=true, merge_results=true))
            for (a, b) in zip(plain.samples, merged.campaign.samples)
                @test isequal(a.value.table, b.value.table) && isequal(run_table(a.value), run_table(b.value))
                @test b.value.table.spacecraft == fill(Int32(1), nrow(b.value.table))
                @test !("sun_dir_1" in names(b.value.table)) && "sun_dir_1" in names(run_table(b.value))
                @test b.value.table !== merged.members
                @test !isempty(sprint(show, run_table(b.value)))                   # not left inconsistent
            end
            @test nrow(merged.members) == sum(s -> nrow(s.value.table), merged.campaign.samples)
            @test sort(unique(merged.members.spacecraft)) == Int32[1, 2, 3]
        end
    end

    @testset "each member's rows are released as they are merged" begin
        args = SE._resolve_effector_scopes(configuration([three_body(k) for k in 1:3]; results=false))
        configs = [SC._ensemble_member_configuration(args, args.dynamics_model.spacecraft[k], k, 3) for k in 1:3]
        campaign = quiet(() -> SC._run_constellation_members(args, configs, 1, false, nothing, nothing,
            pairs((;)), true))
        tables = [WeakRef(sample.value.spacecraft) for sample in campaign.samples]
        res = SC._constellation_ensemble_result(args, configs, campaign)
        GC.gc(true)
        # The campaign's samples hold only the caller's values now, and the tables
        # of the members appended after the first are gone.
        @test all(sample -> sample.value === nothing, campaign.samples)
        @test all(sample -> sample.value === nothing, campaign.successful)
        @test tables[1].value === res.members
        @test tables[2].value === nothing && tables[3].value === nothing
    end

    @testset "members with different columns" begin
        # With the scene saved, a member with links has link poses; one without has none.
        spacecraft = [three_body(1), one_link(2), three_body(3)]
        res = quiet(() -> run_constellation_ensemble(configuration(spacecraft; save_visualization_scene=true);
            merge_results=true))
        poses = [name for name in names(res.members) if startswith(name, "link_pose")]
        @test !isempty(poses)
        @test all(name -> all(ismissing, res.members[res.members.spacecraft .== 2, name]), poses)
        @test all(name -> !any(ismissing, res.members[res.members.spacecraft .!= 2, name]), poses)
        # First seen: member 1's columns, in its order.
        own = names(expected_rows(configuration(spacecraft; save_visualization_scene=true), 1))
        @test names(res.members)[1:length(own)] == own
    end

    @testset "a failed member" begin
        spacecraft = [three_body(k) for k in 1:4]
        effectors = (SM.InverseSquaredJ2GravityModel(), scoped(Fails(), BoundSpacecraft((3,))))
        dir = mktempdir()
        res = quiet(() -> run_constellation_ensemble(configuration(spacecraft; effectors=effectors, dir=dir);
            merge_results=true))
        row = res.constellation[3, :]
        @test !row.success && occursin("this member's model failed", row.error) && row.rows == 0
        @test ismissing(row.time) && ismissing(row.pos_1)
        @test !(3 in res.members.spacecraft) && sort(unique(res.members.spacecraft)) == Int32[1, 2, 4]
        @test all(res.constellation.success[[1, 2, 4]]) && all(ismissing, res.constellation.error[[1, 2, 4]])
        @test res.constellation.rows[[1, 2, 4]] == [count(==(k), res.members.spacecraft) for k in (1, 2, 4)]
        @test length(res.campaign.failed) == 1 && res.campaign.failed[1].index == 3
        # The summary's last row is the member's last row.
        last_two = res.members[res.members.spacecraft .== 2, :][end, :]
        @test res.constellation.time[2] == last_two.time && res.constellation.pos_1[2] == last_two.pos_1
        # With fail_fast the call throws as before, and writes no merged file.
        failing_dir = mktempdir()
        @test_throws Exception quiet(() -> run_constellation_ensemble(
            configuration(spacecraft; effectors=effectors, dir=failing_dir); merge_results=true, fail_fast=true))
        @test !isfile(joinpath(failing_dir, "constellation_summary.csv"))
    end

    @testset "files" begin
        dir = mktempdir()
        args = configuration(fleet()[1:3]; dir=dir, save_visualization_scene=true)
        res = quiet(() -> run_constellation_ensemble(args; merge_results=true))
        paths = joinpath.(dir, ["constellation_members.feather", "constellation_members.csv", "constellation_summary.csv"])
        @test res.files == paths && all(isfile, paths)
        @test isequal(DataFrame(Arrow.Table(paths[1])), res.members)
        @test nrow(CSV.read(paths[2], DataFrame)) == nrow(res.members)
        summary = CSV.read(paths[3], DataFrame)
        @test summary.member_tag == ["sat_1_id_13", "sat_2_id_23", "sat_3_id_33"] && summary.rows == res.constellation.rows
        # The member folders are what the ensemble writes without the keyword, and
        # the viewer's discovery finds the same samples (labelled by identity from
        # the summary, see below).
        plain_dir = mktempdir()
        quiet(() -> run_constellation_ensemble(configuration(fleet()[1:3]; dir=plain_dir, save_visualization_scene=true)))
        for tag in ("sat_1_id_13", "sat_2_id_23", "sat_3_id_33")
            @test read(joinpath(dir, tag, "simulation_results.csv")) == read(joinpath(plain_dir, tag, "simulation_results.csv"))
        end
        discover = SM.SceneVisualization.discover_ensemble_samples
        @test [basename(s.directory) for s in first(discover(dir))] ==
            [basename(s.directory) for s in first(discover(plain_dir))] == ["sat_1_id_13", "sat_2_id_23", "sat_3_id_33"]
        # save_csv=false: no members CSV.
        no_csv_dir = mktempdir()
        no_csv = quiet(() -> run_constellation_ensemble(configuration(fleet()[1:2]; dir=no_csv_dir, save_csv=false);
            merge_results=true))
        @test !isfile(joinpath(no_csv_dir, "constellation_members.csv"))
        @test isfile(joinpath(no_csv_dir, "constellation_members.feather")) && length(no_csv.files) == 2
        # results=false: nothing written.
        @test isempty(quiet(() -> run_constellation_ensemble(configuration(fleet()[1:2]; results=false);
            merge_results=true)).files)
    end

    @testset "refusals before any member runs" begin
        started = Ref(0)
        counting(args) = (started[] += 1; OrdinaryDiffEq.DiscreteCallback((u, t, integrator) -> false, integrator -> nothing))
        args = configuration(fleet()[1:2])
        @test_throws "cannot return solver metadata" run_constellation_ensemble(args; merge_results=true,
            return_solver_metadata=true, extra_callback_factories=(counting,))
        clash = [SM.SaveField(:pos_again, (u, t, integrator) -> t; per_satellite=false, column_prefix="pos")]
        @test_throws "has the column prefix \"pos\"" run_constellation_ensemble(args; merge_results=true,
            save_fields=vcat(SM.default_save_fields(member(args, 1)), clash), extra_callback_factories=(counting,))
        @test_throws "results_layout must be :wide or :long" run_constellation_ensemble(
            configuration(fleet()[1:2]; results_layout=:tall); extra_callback_factories=(counting,))
        bound = configuration(fleet()[1:2]; effectors=(SM.InverseSquaredJ2GravityModel(), Bound(2)))
        @test_throws "binds it to member 2" run_constellation_ensemble(bound; merge_results=true,
            extra_callback_factories=(counting,))
        @test started[] == 0
        # Without the keyword, the metadata pair runs as before.
        @test quiet(() -> run_constellation_ensemble(args; return_solver_metadata=true)) isa MonteCarloResult
    end

    @testset "run-wide and per-spacecraft names" begin
        args = configuration(fleet()[1:2])
        # Save fields for the members, which run one spacecraft each.
        defaults = SM.default_save_fields(member(args, 1))
        # A run-wide column that meets a per-spacecraft one only once flattened.
        late = [SM.SaveField(:late, (u, t, integrator) -> t; per_satellite=false, column_prefix="pos_1")]
        @test_throws "run-wide column \"pos_1\"" quiet(() -> run_constellation_ensemble(args; merge_results=true,
            save_fields=vcat(defaults, late)))
        # A run-wide field of its own name joins the members table under it.
        clock = [SM.SaveField(:clock, (u, t, integrator) -> t^2)]
        res = quiet(() -> run_constellation_ensemble(args; merge_results=true, save_fields=vcat(defaults, clock)))
        @test all(res.constellation.success)
        @test res.members.clock == res.members.time .^ 2
        @test findfirst(==("clock"), names(res.members)) > findfirst(==("heat_load"), names(res.members))
    end

    @testset "members merged from disk" begin
        dir = mktempdir()
        res = quiet(() -> run_constellation_ensemble(configuration(fleet(); dir=dir, save_visualization_scene=true);
            merge_results=true))
        in_memory = quiet(() -> run_constellation_ensemble(configuration(fleet(); results=false,
            save_visualization_scene=true); merge_results=true))
        # Written member by member and read in place, the table built in memory.
        @test !(res.members.pos_1 isa Vector) && res.members.pos_1 isa AbstractVector
        @test isequal(res.members, in_memory.members)
        @test isequal(select(res.constellation, Not(:elapsed_s)), select(in_memory.constellation, Not(:elapsed_s)))
        @test in_memory.members.pos_1 isa Vector
        # The merged files the ensemble wrote are not replaced unless asked.
        before = read(joinpath(dir, "constellation_summary.csv"))
        @test_throws "already holds merged constellation results" read_constellation_results(dir)
        @test read(joinpath(dir, "constellation_summary.csv")) == before
        # The same tables from the member folders, without running anything; the
        # members' status comes from the summary the ensemble wrote.
        merged = read_constellation_results(dir; overwrite=true)
        @test isequal(merged.members, in_memory.members)
        @test isequal(select(merged.constellation, Not([:elapsed_s, :read_status])), select(res.constellation, Not(:elapsed_s)))
        @test merged.constellation.elapsed_s ≈ res.constellation.elapsed_s
        @test all(==("read"), merged.constellation.read_status)
        @test merged.files == joinpath.(dir, ["constellation_members.feather", "constellation_summary.csv"])
        @test isequal(DataFrame(Arrow.Table(merged.files[1])), in_memory.members)
        # Members written in the long layout give the same table.
        long_dir = mktempdir()
        quiet(() -> run_constellation_ensemble(configuration(fleet(); dir=long_dir, save_visualization_scene=true,
            results_layout=:long)))
        @test isequal(read_constellation_results(long_dir).members, in_memory.members)
        # Without a summary, wide members have no plane or slot, and their status
        # is not known.
        rm(joinpath(dir, "constellation_summary.csv"))
        bare = read_constellation_results(dir; overwrite=true)
        @test all(ismissing, bare.members.plane) && all(ismissing, bare.members.slot)
        @test isequal(select(bare.members, Not([:plane, :slot])), select(in_memory.members, Not([:plane, :slot])))
        @test all(ismissing, bare.constellation.success) && all(ismissing, bare.constellation.error)
        @test all(ismissing, bare.constellation.elapsed_s) && all(==("read"), bare.constellation.read_status)
        # A member the summary records as failed has no rows, although it left
        # partial results.
        failed_dir = mktempdir()
        spacecraft = [three_body(k) for k in 1:3]
        quiet(() -> run_constellation_ensemble(configuration(spacecraft; dir=failed_dir,
            effectors=(SM.InverseSquaredJ2GravityModel(), scoped(Fails(), BoundSpacecraft((2,))))); merge_results=true))
        @test isfile(joinpath(failed_dir, "sat_2_id_2", "simulation_results.feather"))
        partial = read_constellation_results(failed_dir; overwrite=true)
        @test sort(unique(partial.members.spacecraft)) == Int32[1, 3] && partial.constellation.rows[2] == 0
        @test !partial.constellation.success[2] && occursin("this member's model failed", partial.constellation.error[2])
        @test partial.constellation.read_status == ["read", "member failed", "read"]
        @test_throws "No constellation member folders" read_constellation_results(mktempdir())
    end

    @testset "members merged from disk keep every member, and read only what they can" begin
        dir = mktempdir()
        res = quiet(() -> run_constellation_ensemble(configuration(fleet()[1:4]; dir=dir); merge_results=true))
        tags = res.constellation.member_tag
        # A member whose folder is gone keeps its row; a member whose results file
        # is gone is not read. Both are named in a warning.
        rm(joinpath(dir, tags[2]); recursive=true)
        rm(joinpath(dir, tags[3], "simulation_results.feather"))
        merged = @test_logs (:warn, Regex("$(tags[2]) \\(folder missing\\).*$(tags[3]) \\(files differ from manifest\\)")) match_mode=:any read_constellation_results(dir; overwrite=true)
        @test merged.constellation.member_tag == tags && merged.constellation.rows == [res.constellation.rows[1], 0, 0, res.constellation.rows[4]]
        @test merged.constellation.read_status == ["read", "folder missing", "files differ from manifest", "read"]
        @test merged.constellation.success == res.constellation.success
        @test merged.constellation.spacecraft_id == res.constellation.spacecraft_id
        @test sort(unique(merged.members.spacecraft)) == Int32[1, 4]
        # The summary written keeps the member whose folder is gone.
        @test nrow(CSV.read(joinpath(dir, "constellation_summary.csv"), DataFrame)) == 4
        # A bundle of another version or layout is refused, naming its folder.
        manifest_of(tag) = joinpath(dir, tag, "simulation_results.manifest.toml")
        edited(tag, changes) = (m = TOML.parsefile(manifest_of(tag)); merge!(m, changes);
            open(io -> TOML.print(io, m), manifest_of(tag), "w"))
        original = read(manifest_of(tags[4]))
        edited(tags[4], Dict("schema_version" => "9", "layout" => "columnar"))
        @test_throws Regex("$(tags[4]).*schema version \"9\"") read_constellation_results(dir; overwrite=true)
        # So is a wide bundle of more than one spacecraft.
        write(manifest_of(tags[4]), original)
        edited(tags[4], Dict("spacecraft_count" => 2))
        @test_throws "wide results bundle of 2 spacecraft" read_constellation_results(dir; overwrite=true)
        # A results file whose size is not the one its manifest records is not read.
        write(manifest_of(tags[4]), original)
        feather_path = joinpath(dir, tags[4], "simulation_results.feather")
        feather_bytes = read(feather_path)
        write(feather_path, vcat(feather_bytes, UInt8[0]))
        grown = quiet(() -> read_constellation_results(dir; overwrite=true))
        @test grown.constellation.read_status[4] == "files differ from manifest" && grown.constellation.rows[4] == 0
        write(feather_path, feather_bytes)
        # A run that streamed its results and did not finish has no rows.
        edited(tags[4], Dict("streaming" => true))
        unfinished = quiet(() -> read_constellation_results(dir; overwrite=true))
        @test unfinished.constellation.read_status[4] == "unfinished run" && unfinished.constellation.rows[4] == 0
        # Members without a results bundle are not read, with a warning.
        bundle_off = mktempdir()
        withenv("SPACEAGORA_SAVE_BUNDLE" => "0") do
            quiet(() -> run_constellation_ensemble(configuration(fleet()[1:2]; dir=bundle_off)))
        end
        none = @test_logs (:warn, r"no results bundle") match_mode=:any read_constellation_results(bundle_off)
        @test none.constellation.read_status == ["no results bundle", "no results bundle"] && nrow(none.members) == 0
    end

    @testset "the viewer labels members by identity" begin
        dir = mktempdir()
        quiet(() -> run_constellation_ensemble(configuration(fleet(); dir=dir, save_visualization_scene=true);
            merge_results=true))
        samples, scalar_name = SM.SceneVisualization.discover_ensemble_samples(dir)
        @test scalar_name == "plane"
        @test [s.label for s in samples] ==
            ["id 13, plane 1, slot 1", "id 23, plane 1, slot 2", "id 33, plane 2, slot 1", "id 43, plane 2, slot 2",
             "id 53, plane 3, slot 1", "id 63, plane 3, slot 2", "id 99"]
        @test isequal([s.scalar for s in samples], [1.0, 1.0, 2.0, 2.0, 3.0, 3.0, NaN])
        # Without the summary, as before: folder names and no scalar.
        rm(joinpath(dir, "constellation_summary.csv"))
        samples, scalar_name = SM.SceneVisualization.discover_ensemble_samples(dir)
        @test scalar_name == "" && [s.label for s in samples] == sort(filter(n -> startswith(n, "sat_"), readdir(dir)))
        @test all(s -> isnan(s.scalar), samples)
    end

    @testset "a resumed ensemble" begin
        dir = mktempdir()
        spacecraft = [three_body(k) for k in 1:3]
        first_leg = configuration(spacecraft; dir=dir, mission_time=200.0, checkpoint_enabled=true, checkpoint_interval_s=100.0)
        quiet(() -> run_constellation_ensemble(first_leg))
        resumed = configuration(spacecraft; dir=dir, mission_time=400.0, checkpoint_enabled=true,
            checkpoint_interval_s=100.0, resume_from_checkpoint=true)
        res = quiet(() -> run_constellation_ensemble(resumed; merge_results=true))
        # Each member's rows are its whole run, as its own results file holds it.
        for (k, tag) in enumerate(res.constellation.member_tag)
            written = CSV.read(joinpath(dir, tag, "simulation_results.csv"), DataFrame)
            rows = res.members[res.members.spacecraft .== k, :]
            @test rows.time == written.time && rows.pos_1 ≈ written.sc1_pos_1
            @test first(rows.time) == 0.0 && last(rows.time) == 400.0
        end
    end
end

end # module ConstellationResultsTests
