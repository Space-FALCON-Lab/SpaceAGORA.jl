module LongLayoutTests
# `results_layout=:long`: the spacecraft and run tables hold the wide table's
# values rearranged, the files and manifest, the viewer reading a long run, and
# the refusals.
using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel
using DataFrames
using Arrow
using CSV
using JSON
using TOML
using Logging
using StaticArrays
import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

const SM = SpaceAGORA.SimulationModel
const SV = SM.SceneVisualization
const PLANET = make_no_gram_planet(:earth)
const IDENTITY = ["spacecraft", "spacecraft_id", "plane", "slot", "time"]

# A Walker constellation of `n` spacecraft in 3 planes, or `n` spacecraft without
# a constellation slot (`walker=false`), written to `dir` in `layout`.
function config(; dir, layout=:wide, attitude=false, n=6, walker=true, mission_time=300.0, scene=false)
    spacecraft(index) = make_three_body_spacecraft(bus_dims=(1.0, 1.0, 1.0), panel_dims=(0.01, 1.0, 0.5),
        bus_mass=200.0, panel_mass_each=5.0, panel_offset_y=1.0,
        ic=InitialCondition(ra=7.0e6, rp=7.0e6, i=53.0, ω=0.0, Ω=0.0, ν=30.0 * index), prop_mass=0.0, id=index)
    fleet = walker ? walker_constellation((plane, slot, index) -> spacecraft(index);
        total=n, planes=3, phasing=1, altitude=6.0e5, planet=PLANET, i=53) : [spacecraft(i) for i in 1:n]
    base = make_example_config(planet=PLANET, spacecraft=fleet[1], mission_time=mission_time,
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(SM.InverseSquaredJ2GravityModel(),), density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=attitude, keplerian=true, EI_km=120.0,
        verbose=false, results=true, results_directory=dir)
    settings = SM.SimConfig._with_fields(base.simulation_settings; results_layout=layout, save_visualization_scene=scene)
    return SM.SimConfig._with_configuration(base; simulation_settings=settings,
        dynamics_model=SM.DynamicsModel(fleet, base.dynamics_model.dynamic_effectors))
end

# Save fields beyond the defaults: a per-spacecraft NamedTuple, a per-spacecraft
# Dict, a run-wide scalar and a run-wide NamedTuple.
function extra_fields(args)
    return vcat(SM.default_save_fields(args), [
        SM.SaveField(:pair, (u, t, integrator) -> [(a=Float64(i), b=t + i) for i in eachindex(u.sc)]; per_satellite=true),
        SM.SaveField(:keyed, (u, t, integrator) -> [Dict("x" => i, "y" => 2 * i) for i in eachindex(u.sc)];
            per_satellite=true, column_prefix="keyed"),
        SM.SaveField(:clock, (u, t, integrator) -> t^2),
        SM.SaveField(:summary, (u, t, integrator) -> (count=length(u.sc), half=t / 2)),
    ])
end

quiet(f) = with_logger(f, ConsoleLogger(stderr, Logging.Error))

# A control model that holds no command and reports thruster levels for one spacecraft.
struct ReportedThrusters <: SM.AbstractControlEffectorModel
    firing::Int
end
SpaceAGORA.update_control!(::ReportedThrusters, u, p, t::Float64, i::Int64) = nothing
SpaceAGORA.control_wrench(::ReportedThrusters, x, p, t) = (SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0))
SM.ControlHooks.control_thruster_levels(c::ReportedThrusters, i::Int) = i == c.firing ? [0.25, 0.75] : nothing

# The wide table of a long spacecraft table and run table, in the wide table's
# column order: the per-spacecraft column `sc<i>_<c>` is spacecraft `i`'s rows of `c`.
function pivot(long::DataFrame, run, wide_names)
    S = maximum(long.spacecraft)
    out = DataFrame(time=long.time[long.spacecraft .== 1])
    for name in wide_names
        name == "time" && continue
        m = match(r"^sc(\d+)_(.*)$", name)
        out[!, name] = m === nothing ? run[!, name] : long[long.spacecraft .== parse(Int, m[1]), m[2]]
    end
    return out
end

@testset "the long tables hold the wide table's values" begin
    for attitude in (false, true)
        wide_args, long_args = config(dir=mktempdir(); attitude), config(dir=mktempdir(), layout=:long, attitude=attitude)
        fields = extra_fields(wide_args)
        wide = quiet(() -> run_simulation(wide_args; return_results=true, save_fields=fields))
        long = quiet(() -> run_simulation(long_args; return_results=true, save_fields=fields))
        S = length(wide_args.dynamics_model.spacecraft)
        @test S == 6
        # The spacecraft table: identity first, time-major rows, no `sc<i>_` names.
        @test names(long.table)[1:5] == IDENTITY
        @test nrow(long.table) == S * nrow(wide.table)
        @test long.table.spacecraft == repeat(Int32.(1:S), nrow(wide.table))
        @test long.table.time == repeat(wide.table.time; inner=S)
        @test !any(name -> startswith(name, "sc"), names(long.table)[6:end])
        # The run table: time and the run-wide fields only.
        @test names(run_table(long)) == ["time"; [n for n in names(wide.table) if n != "time" && !startswith(n, "sc")]]
        @test "clock" in names(run_table(long)) && "summary_half" in names(run_table(long))
        @test run_table(wide) === nothing
        # Pivoted back, the wide table, values and element types alike.
        back = pivot(long.table, run_table(long), names(wide.table))
        @test isequal(back, wide.table)
        @test all(name -> eltype(back[!, name]) == eltype(wide.table[!, name]), names(wide.table))
        # The NamedTuple and Dict fields flatten as in the wide table.
        @test all(name -> name in names(long.table), ["pair_a", "pair_b", "keyed_x", "keyed_y"])
        @test "sc3_pair_a" in names(wide.table) && "sc3_keyed_y" in names(wide.table)
        # long_table rearranges a wide run into the same spacecraft table.
        @test isequal(long_table(wide), long.table)
        @test long_table(long) === long.table
    end
end

@testset "fields whose shape differs between spacecraft" begin
    rows(table, k) = table[table.spacecraft .== k, :]
    # A user field with no entries for one spacecraft and two for the other, in
    # both orders, and a NamedTuple field whose keys differ between them.
    for short in (1, 2)
        other = 3 - short
        ragged = SM.SaveField(:rag, (u, t, integrator) -> [i == short ? Float64[] : [t, 2t] for i in 1:2];
            per_satellite=true, column_prefix="rag")
        keyed = SM.SaveField(:named, (u, t, integrator) -> [i == short ? (a=t,) : (a=t, b=2t) for i in 1:2];
            per_satellite=true, column_prefix="named")
        wide_args = config(dir=mktempdir(), n=2, walker=false, scene=true)
        long_args = config(dir=mktempdir(), layout=:long, n=2, walker=false, scene=true)
        gap = SM.SaveField(:gap, (u, t, integrator) -> [i == short ? missing : t for i in 1:2];
            per_satellite=true, column_prefix="gap")
        fields = vcat(SM.default_save_fields(wide_args), [ragged, keyed, gap])
        wide = quiet(() -> run_simulation(wide_args; return_results=true, save_fields=fields))
        long = quiet(() -> run_simulation(long_args; return_results=true, save_fields=fields))
        @test "sc$(other)_rag_2" in names(wide.table) && !("sc$(short)_rag_1" in names(wide.table))
        @test "sc$(other)_named_b" in names(wide.table) && !("sc$(short)_named_b" in names(wide.table))
        for name in ("rag_1", "rag_2", "named_b")
            @test rows(long.table, other)[!, name] == wide.table[!, "sc$(other)_$(name)"]
            @test all(ismissing, rows(long.table, short)[!, name])
        end
        @test rows(long.table, short).named_a == wide.table[!, "sc$(short)_named_a"]
        @test isequal(long_table(wide), long.table)
        # The viewer offers a channel only when every spacecraft has its column. The
        # long run's wide view leaves out the columns a spacecraft lacks, as its wide
        # table has none, so the same channels are offered.
        channels = [(column="rag_2",), (column="named_b",), (column="named_a",), (column="gap",)]
        payload(args) = JSON.parsefile(export_viewer_payload(args; textures=false, channels=channels))
        wide_json, long_json = payload(wide_args), payload(long_args)
        @test [channel["name"] for channel in long_json["frames"]["channels"]] == ["named_a", "gap"]
        @test isequal(wide_json["frames"], long_json["frames"])
        # In-memory long results and a converted wide table retain the same schema,
        # including present all-missing fields and absent ragged components.
        for table in (long.table, long_table(wide))
            restored = SV._wide_view(table, run_table(long), 2)
            @test Set(names(restored)) == Set(names(wide.table))
            @test all(name -> isequal(restored[!, name], wide.table[!, name]), names(wide.table))
        end
        legacy = copy(long.table)
        emptymetadata!(legacy)
        @test !("sc$(short)_gap" in names(SV._wide_view(legacy, run_table(long), 2)))
        damaged = copy(long.table)
        metadata!(damaged, "spaceagora.long_column_spacecraft.v1", "{}"; style=:note)
        @test_throws ArgumentError SV._wide_view(damaged, run_table(long), 2)
    end
    # The thruster levels of a spacecraft with thrusters beside one without, in
    # both orders; the viewer reads the long run as the wide one.
    for with_thrusters in (1, 2)
        without = 3 - with_thrusters
        dirs = (mktempdir(), mktempdir())
        args = (config(dir=dirs[1], n=2, walker=false, scene=true),
                config(dir=dirs[2], layout=:long, n=2, walker=false, scene=true))
        for a in args, location in ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
            push!(a.dynamics_model.spacecraft[with_thrusters].root.thrusters, SM.Thruster(max_thrust=2.0,
                location=MVector{3, Float64}(location...), direction=MVector{3, Float64}(1.0, 0.0, 0.0)))
        end
        control = SM.ControlModel(control_effectors=(ReportedThrusters(with_thrusters),), control_rates=[10.0])
        args = map(a -> SM.SimConfig._with_configuration(a; control_model=control), args)
        wide = quiet(() -> run_simulation(args[1]; return_results=true))
        long = quiet(() -> run_simulation(args[2]; return_results=true))
        @test "sc$(with_thrusters)_thruster_level_2" in names(wide.table)
        @test !("sc$(without)_thruster_level_1" in names(wide.table))
        @test rows(long.table, with_thrusters).thruster_level_2 == wide.table[!, "sc$(with_thrusters)_thruster_level_2"]
        @test all(==(0.75), rows(long.table, with_thrusters).thruster_level_2[2:end])
        @test all(ismissing, rows(long.table, without).thruster_level_1)
        @test isequal(long_table(wide), long.table)
        payload(dir) = JSON.parsefile(export_viewer_payload(joinpath(dir, "simulation_results"); textures=false))
        wide_json, long_json = payload(dirs[1]), payload(dirs[2])
        @test wide_json["frames"]["thruster_counts"] == long_json["frames"]["thruster_counts"]
        @test isequal(wide_json["frames"], long_json["frames"])
    end
end

@testset "identity columns" begin
    dir = mktempdir()
    walker = quiet(() -> run_simulation(config(dir=dir, layout=:long, n=6); return_results=true)).table
    rows = walker[walker.time .== 0.0, :]
    @test rows.spacecraft == Int32.(1:6) && rows.spacecraft_id == 1:6
    @test rows.plane == Int16[1, 1, 2, 2, 3, 3] && rows.slot == Int16[1, 2, 1, 2, 1, 2]
    @test eltype(walker.plane) == Union{Missing, Int16} && eltype(walker.spacecraft_id) == Int64
    free = quiet(() -> run_simulation(config(dir=mktempdir(), layout=:long, n=2, walker=false); return_results=true)).table
    @test all(ismissing, free.plane) && all(ismissing, free.slot) && eltype(free.slot) == Union{Missing, Int16}
    # A repeated id: the spacecraft column still tells the rows apart.
    args = config(dir=mktempdir(), layout=:long, n=2, walker=false)
    args.dynamics_model.spacecraft[2].id = 1
    same_id = quiet(() -> run_simulation(args; return_results=true)).table
    @test same_id.spacecraft_id == fill(1, nrow(same_id)) && same_id.spacecraft == repeat(Int32[1, 2], nrow(same_id) ÷ 2)
end

@testset "files and manifest" begin
    dir = mktempdir()
    args = config(dir=dir, layout=:long)
    result = quiet(() -> run_simulation(args; return_results=true, save_fields=extra_fields(args)))
    prefix = joinpath(dir, "simulation_results")
    for path in (prefix * "_spacecraft.csv", prefix * "_spacecraft.feather", prefix * "_run.csv", prefix * "_run.feather",
                 prefix * ".manifest.toml")
        @test isfile(path) && path in result.files
    end
    @test !isfile(prefix * ".csv") && !isfile(prefix * ".feather")
    manifest = TOML.parsefile(prefix * ".manifest.toml")
    @test manifest["schema_version"] == "2" && manifest["layout"] == "long"
    @test manifest["steps"] == nrow(run_table(result)) && manifest["spacecraft_count"] == 6
    @test sort(collect(keys(manifest["files"]))) == ["run_csv", "run_feather", "spacecraft_csv", "spacecraft_feather"]
    for (key, entry) in manifest["files"]
        @test entry["size_bytes"] == filesize(entry["path"])
        @test entry["sha256"] == SpaceAGORA.SimulationModel.IOSerialization._sha256_hex(entry["path"])
    end
    @test isequal(DataFrame(Arrow.Table(prefix * "_spacecraft.feather")), result.table)
    @test isequal(DataFrame(Arrow.Table(prefix * "_run.feather")), run_table(result))
    @test nrow(CSV.read(prefix * "_spacecraft.csv", DataFrame)) == nrow(result.table)
    # The CSVs follow save_csv.
    no_csv_dir = mktempdir()
    no_csv = config(dir=no_csv_dir, layout=:long)
    no_csv = SM.SimConfig._with_configuration(no_csv; simulation_settings=SM.SimConfig._with_fields(no_csv.simulation_settings; save_csv=false))
    quiet(() -> run_simulation(no_csv))
    @test isfile(joinpath(no_csv_dir, "simulation_results_spacecraft.feather"))
    @test !isfile(joinpath(no_csv_dir, "simulation_results_spacecraft.csv"))
    # Without a run-wide field: no run table, and an earlier one is removed.
    fields = filter(field -> field.per_satellite, SM.default_save_fields(args))
    only_spacecraft = quiet(() -> run_simulation(args; return_results=true, save_fields=fields))
    @test run_table(only_spacecraft) === nothing
    @test !isfile(prefix * "_run.feather") && !isfile(prefix * "_run.csv")
    @test !haskey(TOML.parsefile(prefix * ".manifest.toml")["files"], "run_feather")
end

@testset "the viewer reads a long run as the wide one" begin
    for n in (3, 66)
        wide_dir, long_dir = mktempdir(), mktempdir()
        quiet(() -> run_simulation(config(dir=wide_dir, n=n, scene=true)))
        quiet(() -> run_simulation(config(dir=long_dir, layout=:long, n=n, scene=true)))
        wide_scene = SV.read_visualization_scene(joinpath(wide_dir, "simulation_results_scene.json"))
        long_scene = SV.read_visualization_scene(joinpath(long_dir, "simulation_results_scene.json"))
        @test wide_scene.results_layout === :wide && long_scene.results_layout === :long
        @test long_scene.results_feather == "simulation_results_spacecraft.feather"
        @test !occursin("\"layout\":\"long\"", read(joinpath(wide_dir, "simulation_results_scene.json"), String))
        wide_payload = export_viewer_payload(joinpath(wide_dir, "simulation_results"); textures=false)
        long_payload = export_viewer_payload(joinpath(long_dir, "simulation_results"); textures=false)
        wide_json, long_json = JSON.parsefile(wide_payload), JSON.parsefile(long_payload)
        @test isequal(wide_json["frames"], long_json["frames"])
        # The scene differs only where it names the results file and its layout.
        for json in (wide_json, long_json)
            delete!(json["scene"]["results"], "feather")
            delete!(json["scene"]["results"], "layout")
        end
        @test isequal(wide_json, long_json)
    end
    # The wide view reads the rows without copying them.
    dir = mktempdir()
    quiet(() -> run_simulation(config(dir=dir, layout=:long, n=3, scene=true)))
    scene = SV.read_visualization_scene(joinpath(dir, "simulation_results_scene.json"))
    view_df = SV._saved_results_table(scene, joinpath(dir, scene.results_feather))
    @test view_df[!, "sc2_pos_1"] isa SubArray && nrow(view_df) == nrow(DataFrame(Arrow.Table(joinpath(dir, "simulation_results_run.feather"))))
    # Rows in another order are read by their key.
    long = DataFrame(Arrow.Table(joinpath(dir, scene.results_feather)))
    shuffled = long[sortperm(collect(zip(long.spacecraft, .-long.time))), :]
    @test isequal(SV._wide_view(shuffled, nothing, 3)[!, "sc3_pos_2"], SV._wide_view(long, nothing, 3)[!, "sc3_pos_2"])
    @test_throws ArgumentError SV._wide_view(long[2:end, :], nothing, 3)
end

@testset "refusals" begin
    args = config(dir=mktempdir())
    bad = SM.SimConfig._with_configuration(args; simulation_settings=SM.SimConfig._with_fields(args.simulation_settings; results_layout=:tall))
    @test_throws ArgumentError run_simulation(bad)
    campaign = mktempdir()
    # The sample fails with the message, whether the campaign rethrows it or
    # records it on the sample.
    refusal = try
        result, _ = quiet(() -> run_monte_carlo_visualization(seed -> config(dir=mktempdir(), layout=:long, n=2, walker=false),
            [1], campaign; threads=1, fail_fast=true, export_page=false))
        sample = only(result.samples)
        sample.success ? "" : (sample.error isa Exception ? sprint(showerror, sample.error) : string(sample.error))
    catch err
        sprint(showerror, err)
    end
    @test occursin("needs results_layout=:wide", refusal)
end

@testset "long output refuses fields that replace row identities" begin
    args = (dynamics_model=(spacecraft=[(id=11, constellation_slot=nothing), (id=22, constellation_slot=nothing)],),
            simulation_settings=(results_layout=:wide,))
    times = [0.0, 10.0]
    cases = Any[(name, [[100.0, 200.0], [110.0, 210.0]]) for name in
             ("spacecraft", "spacecraft_id", "plane", "slot", "time")]
    push!(cases, ("spacecraft", [[(id=100.0,), (id=200.0,)], [(id=110.0,), (id=210.0,)]]))
    push!(cases, ("spacecraft", [[Dict("id" => 100.0), Dict("id" => 200.0)],
                                [Dict("id" => 110.0), Dict("id" => 210.0)]]))
    for (prefix, values) in cases
        field = SM.SaveField(:custom, (u, t, i) -> nothing; per_satellite=true, column_prefix=prefix)
        rows = [Dict(:custom => value) for value in values]
        wide = SM.IOOutputs._build_results_dataframe(times, rows, [field], args)
        before = deepcopy(wide)
        @test wide.time == times
        @test_throws ArgumentError SM.IOOutputs._build_long_results_tables(times, rows, [field], args)
        @test_throws ArgumentError long_table(SimulationResults(wide, args, String[], nothing))
        @test isequal(wide, before)
    end
    # The flattened name, not the prefix, determines whether a collision exists.
    field = SM.SaveField(:custom, (u, t, i) -> nothing; per_satellite=true, column_prefix="spacecraft")
    rows = [Dict(:custom => [(temperature=2.0,), (temperature=3.0,)]) for _ in times]
    safe = SM.IOOutputs._build_long_results_tables(times, rows, [field], args).spacecraft
    @test safe.spacecraft == Int32[1, 2, 1, 2]
    @test safe.spacecraft_temperature == [2.0, 3.0, 2.0, 3.0]
    run_field = SM.SaveField(:custom, (u, t, i) -> nothing; column_prefix="time")
    @test_throws ArgumentError SM.IOOutputs._build_long_results_tables(times,
        [Dict(:custom => 100.0), Dict(:custom => 110.0)], [run_field], args)
end

end # module LongLayoutTests
