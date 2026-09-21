using Test, Arrow
include("replay_feather_animation.jl")
const Replay = FeatherAnimationReplay

@testset "Replay startup routing" begin
    path = joinpath(@__DIR__, "replay_feather_animation.jl")
    stub = :(module FeatherAnimationReplay
        const calls = Any[]
        main(arguments=Main.ARGS) = push!(calls, arguments)
    end)
    for (interactive, program, expected) in ((true, "", 1), (false, "", 0),
                                            (false, "another_script.jl", 0), (false, path, 1))
        sandbox = Module(gensym(:ReplayStartup))
        Core.eval(sandbox, :(const PROGRAM_FILE = $program))
        Core.eval(sandbox, :(isinteractive() = $interactive))
        Base.include(expression -> expression isa Expr && expression.head == :module ?
            deepcopy(stub) : expression, sandbox, path)
        calls = Base.invokelatest(() -> getproperty(getproperty(sandbox, :FeatherAnimationReplay), :calls))
        @test length(calls) == expected
        interactive && @test only(calls) == String[]
    end
    @test Replay.DEFAULT_BUNDLE == Replay.OPTIONS.bundle_directory
    @test Replay.OPTIONS.duration_seconds == 100.0
    @test Replay.OPTIONS.animation_fps == 30
end

@testset "Feather animation replay" begin
    trajectory = Replay.ReplayTrajectory([0.0, 10.0],
        [[0.0, 0, 0, 0, 0, 0], [100.0, 0, 0, 20, 0, 0]])
    @test trajectory(0.0) == first(trajectory.u)
    @test trajectory(10.0) == last(trajectory.u)
    @test trajectory(2.5) ≈ [6.25, 0, 0, 5, 0, 0]
    @test_throws ArgumentError trajectory(-1.0)
    @test_throws ArgumentError trajectory(11.0)
    @test !isdefined(Replay, :laser_forces)

    for target_id in (1, 3)
        mktempdir() do directory
            times = [0.0, 10.0]
            columns = Dict{Symbol,AbstractVector}(:time=>times)
            for satellite in 1:3, (field, offset) in (("pos", 0), ("vel", 3)), component in 1:3
                state = [Replay.R_EARTH + (satellite == target_id ? 1050e3 : 1000e3),
                    satellite*1000.0, 0.0, 0.0, 7000.0, 0.0]
                columns[Symbol("sc$(satellite)_$(field)_$(component)")] = fill(state[offset+component], 2)
            end
            metadata = Dict("source"=>"Kuang prototype", "scenario"=>"prototype_N20_h1000km_t1150km_ih0_it0_e0_nu0_T60s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt1s")
            target_id == 3 && (metadata["target_id"] = "3")
            Arrow.write(joinpath(directory, "trajectory.feather"), columns; metadata)
            helper_id = target_id == 1 ? 2 : 1
            pair = minmax(helper_id, target_id)
            geometry = (encounter_id=[1], sc_a=[pair[1]], sc_b=[pair[2]],
                start_time_s=[1.0], end_time_s=[10.0], end_clipped=[true])
            lasers = (encounter_id=[1, 1], sc_a=fill(pair[1], 2), sc_b=fill(pair[2], 2),
                start_time_s=[3.0, 8.0], end_time_s=[5.0, 10.0], end_clipped=[false, true])
            Arrow.write(joinpath(directory, "geometry_encounters.feather"), geometry; metadata)
            Arrow.write(joinpath(directory, "laser_on.feather"), lasers; metadata)
            bundle = Replay.load_replay_bundle(directory)
            @test bundle.helper_num == 2
            @test bundle.scenario_caption == "N20_h1000km_t1150km_ih0_it0_e0_nu0_gve_sma"
            @test last(bundle.stored_order) == target_id
            @test bundle.sol.u[1][13] == Replay.R_EARTH + 1050e3
            @test bundle.sol.u[1][1] == Replay.R_EARTH + 1000e3
            for (time, encounter, laser) in ((0.0, false, false), (1.0, true, false),
                    (3.0, true, true), (5.0, true, false), (8.0, true, true), (10.0, true, true))
                links = Replay.recorded_link_states(bundle, time)
                @test links.encounters[1, 3] == encounter
                @test links.active[1, 3] == Float32(laser)
                @test links.kinds[1, 3] == (laser ? 2 : 0)
                @test links.active == transpose(links.active)
                @test all((links.active .== 0) .| links.encounters)
            end
            invalid_lasers = merge(lasers, (start_time_s=[0.0, 8.0],))
            if target_id == 1
                result = Replay.replay_feather_animation(directory; duration_seconds=1.0, animation_fps=11,
                    show_earth=false, output_file=joinpath(directory, "replay.mp4"))
                makie = Replay.GLMakie
                @test filesize(result.output_file) > 0
                axis = makie.content(result.fig[1, 1])
                @test first(split(axis.title[], '\n')) == bundle.scenario_caption
                @test startswith(last(split(axis.title[], '\n')), "t = ")
                @test axis.titlefont[] == :regular
                heading = makie.content(result.fig[0, 1])
                @test !heading.tellwidth[]
                @test heading.text[] == "Multi-Satellite Orbital Animation"
                bounds = axis.layoutobservables.computedbbox[]
                @test bounds.origin ≈ [46.0, 46.0] atol=0.1
                @test bounds.widths ≈ [1108.0, 671.36] atol=0.1
                encounter_plot = only(filter(plot -> haskey(plot.attributes, :label) && plot.label[] == "Encounter detected", axis.scene.plots))
                laser_plot = only(filter(plot -> haskey(plot.attributes, :label) && plot.label[] == "Open-cavity pair", axis.scene.plots))
                legend = only(filter(block -> block isa makie.Legend, result.fig.content))
                encounter_entry = only(filter(entry -> entry.label[] == "Encounter detected", legend.entrygroups[][1][2]))
                laser_entry = only(filter(entry -> entry.label[] == "Laser active", legend.entrygroups[][1][2]))
                for (frame, encounter, laser) in ((1, false, false), (2, true, false), (4, true, true), (6, true, false), (11, true, true))
                    result.controls.current_time[] = frame
                    @test encounter_plot.visible[] == encounter
                    @test laser_plot.visible[] == laser
                    @test encounter_entry.labelcolor[] == (encounter ? :black : :white)
                    @test laser_entry.labelcolor[] == (laser ? :black : :white)
                    @test encounter_plot[1][] == laser_plot[1][]
                end
                @test encounter_plot.linewidth[] == laser_plot.linewidth[] == 2
            end
            Arrow.write(joinpath(directory, "laser_on.feather"), invalid_lasers; metadata)
            @test_throws ArgumentError Replay.load_replay_bundle(directory)
            empty_lasers = (; (name=>values[1:0] for (name, values) in pairs(lasers))...)
            Arrow.write(joinpath(directory, "laser_on.feather"), empty_lasers; metadata)
            without_lasers = Replay.load_replay_bundle(directory)
            @test Replay.recorded_link_states(without_lasers, 3.0).encounters[1, 3]
            @test iszero(sum(Replay.recorded_link_states(without_lasers, 3.0).active))
        end
    end
end