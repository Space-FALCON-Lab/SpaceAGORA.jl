import GLMakie
using GeometryBasics, FileIO
include("test16_feather.jl")
include("functions/10_Animation_ver2.jl")

struct CircularAnimationSolution
    t::Vector{Float64}
    u::Vector{Vector{Float64}}
end

function (sol::CircularAnimationSolution)(time)
    angle = 2pi * time / last(sol.t)
    state = copy(first(sol.u))
    for satellite in 1:2
        radius = norm(first(sol.u)[6*(satellite-1)+1:6*(satellite-1)+3])
        state[6*(satellite-1)+1:6*(satellite-1)+3] = radius .* [cos(angle), sin(angle)/sqrt(2), sin(angle)/sqrt(2)]
        state[6*(satellite-1)+4:6*(satellite-1)+6] = (2pi * radius / last(sol.t)) .* [-sin(angle), cos(angle)/sqrt(2), cos(angle)/sqrt(2)]
    end
    state
end

struct PairHandoffSolution
    t::Vector{Float64}
    u::Vector{Vector{Float64}}
end

(sol::PairHandoffSolution)(time) = sol.u[round(Int, time) + 1]

@testset "Helper-referenced radial scaling with real altitude labels" begin
    orbits = [(a_m=R_EARTH + altitude, e=0.0, i_deg=0.0, Ω_deg=0.0,
        ω_deg=0.0, ν_deg=0.0) for altitude in (1000e3, 1050e3)]
    sol, params, _, _ = mktempdir() do directory
        run_open_cavity_multi(orbits; T_seconds=60.0, helper_num=1,
            Pm=zeros(2, 2), max_range=200e3,
            cavity=Dict{Tuple{Int,Int},Dict{Symbol,Any}}((1, 2)=>Dict(:B=>100.0, :Pin=>10000.0)),
            result_plots=false, verbose=false, IMG_DIR=directory)
    end
    original_states = deepcopy(sol.u)
    reference_radius = R_EARTH + 1000e3
    expected_radius(radius, scale) = scale == 1 || radius <= reference_radius ? radius :
        reference_radius + scale * (radius - reference_radius)
    display_limit = 12000e3
    physical_axis_limit(scale) = (reference_radius + (display_limit - reference_radius) / scale) / 1e3
    axes = []
    earth_bounds = []
    mktempdir() do directory
        for scale in (1.0, 40.0)
            video = joinpath(directory, "scale_$(scale).mp4")
            fig, controls = animate_all_satellites_3d_smooth_helper_target(sol, params, 1;
                radial_exaggeration=scale, reference_altitude_km=1000.0, axis_limit_km=physical_axis_limit(scale),
                radial_ticks_km=[500, 1000, 1010, 1020, 1030, 1040, 1050],
                helper_trails=false, animation_fps=2, duration_seconds=1, output_file=video)
            controls.current_time[] = 1
            axis = GLMakie.content(fig[1, 1])
            markers = filter(plot -> plot isa GLMakie.Scatter, axis.scene.plots)
            @test norm(markers[1][1][][1]) ≈ expected_radius(R_EARTH + 1000e3, scale) atol=4
            @test norm(markers[2][1][][1]) ≈ expected_radius(R_EARTH + 1050e3, scale) atol=4
            link = only(filter(plot -> plot isa GLMakie.Lines && haskey(plot.attributes, :label) &&
                plot.label[] == "Open-cavity pair", axis.scene.plots))
            @test link.alpha[] == 1
            @test link.visible[]
            @test link.linewidth[] == 2
            @test link.color[] == GLMakie.Makie.to_color(:orange)
            @test occursin("ON", GLMakie.content(fig[1, 1]).title[])
            @test link[1][] == [markers[1][1][][1], markers[2][1][][1]]
            earth = only(filter(plot -> plot isa GLMakie.Mesh, axis.scene.plots))
            @test earth.visible[]
            push!(earth_bounds, GLMakie.boundingbox(earth))
            push!(axes, axis.finallimits[])
            if scale == 1
                @test axis.xtickformat[]([-1e6, 0, 1e6]) == ["-1000.0", "0.0", "1000.0"]
            else
                @test axis.xtickformat[] === GLMakie.Makie.automatic
            end
            @test axis.xlabel[] == (scale == 1 ? "X [km]" : "X / geocentric reference [km]")
            @test axis.ylabel[] == (scale == 1 ? "Y [km]" : "Y / geocentric reference [km]")
            @test axis.zlabel[] == (scale == 1 ? "Z [km]" : "Z / geocentric reference [km]")
            @test axis.xticklabelsvisible[]
            @test axis.yticklabelsvisible[]
            @test axis.zticklabelsvisible[]
            @test axis.xgridvisible[]
            @test axis.ygridvisible[]
            @test axis.zgridvisible[]
            @test sol.u == original_states
            @test controls.animation_task === nothing
            @test filesize(video) > 0
            controls.current_time[] = 2
            for satellite in 1:2
                actual_position = sol.u[end][6*(satellite-1)+1:6*(satellite-1)+3]
                expected_position = actual_position .* (expected_radius(norm(actual_position), scale) / norm(actual_position))
                @test collect(only(markers[satellite][1][])) ≈ expected_position atol=4
            end
            controls.current_time[] = 1
            if scale == 40
                @test !any(block -> block isa GLMakie.Axis, fig.content)
                for ticks in (axis.xticks[], axis.yticks[], axis.zticks[])
                    locations, labels = ticks
                    reference_values = [-7400, -7200, -3600, 0, 3600, 7200, 7400]
                    @test locations ≈ [sign(value) * expected_radius(abs(value) * 1e3, 40) for value in reference_values]
                    @test issorted(locations) && all(spacing -> spacing > 0, diff(locations))
                    @test maximum(abs, locations) <= 0.9 * display_limit
                    @test length(labels) == 7
                    @test parse.(Int, labels) == -reverse(parse.(Int, labels))
                    @test labels[5:end] == ["3600", "7200", "7400"]
                    @test locations[4] == 0
                    @test labels[4] == "0"
                    @test all(label -> !isempty(label), labels)
                    for (position, label) in zip(locations, labels)
                        physical_radius = abs(position) <= reference_radius ? abs(position) :
                            reference_radius + (abs(position) - reference_radius) / 40
                        @test parse(Int, label) % 100 == 0
                        @test parse(Float64, label) ≈ sign(position) * physical_radius / 1e3 atol=1e-9
                    end
                end
                @test !any(plot -> haskey(plot.attributes, :label) && occursin("altitude", string(plot.label[])), axis.scene.plots)
                output = joinpath(@__DIR__, "output", "smoke", "videos")
                mkpath(output)
                cp(video, joinpath(output, "radial_exaggeration_test.mp4"); force=true)
                GLMakie.save(joinpath(output, "radial_exaggeration_test.png"), fig)
            end
        end
        single_params = merge(params, Dict(:Pmatrix=>[0.0 10000.0; 0.0 0.0], :cavity=>empty(params[:cavity])))
        single_fig, single_controls = animate_all_satellites_3d_smooth_helper_target(sol, single_params, 1;
            radial_exaggeration=40.0, reference_altitude_km=1000.0,
            animation_fps=2, duration_seconds=1, output_file=joinpath(directory, "single.mp4"))
        single_controls.current_time[] = 1
        single_axis = GLMakie.content(single_fig[1, 1])
        single_link = only(filter(plot -> plot isa GLMakie.Lines && haskey(plot.attributes, :label) &&
            plot.label[] == "Single-pass pair", single_axis.scene.plots))
        @test single_link.color[] == GLMakie.Makie.to_color(:green)
        @test single_link.linewidth[] == 2
        @test single_link.visible[]
        @test single_link.alpha[] == 1
        encounter_states = [vcat([R_EARTH + 1000e3, 100e3, 0, 0, 7000, 0],
            [R_EARTH + 1050e3, target_y, 0, 0, 7000, 0]) for target_y in (0.0, 200e3, 500e3)]
        encounter_sol = PairHandoffSolution([0.0, 1.0, 2.0], encounter_states)
        encounter_params = merge(params, Dict(:gve_schedule=>"gve_sma", :use_los=>false))
        encounter_fig, encounter_controls = animate_all_satellites_3d_smooth_helper_target(encounter_sol, encounter_params, 1;
            radial_exaggeration=40.0, reference_altitude_km=1000.0,
            animation_fps=3, duration_seconds=1, output_file=joinpath(directory, "encounter.mp4"))
        encounter_axis = GLMakie.content(encounter_fig[1, 1])
        encounter_line = only(filter(plot -> plot isa GLMakie.Lines && haskey(plot.attributes, :label) &&
            plot.label[] == "Encounter detected", encounter_axis.scene.plots))
        active_line = only(filter(plot -> plot isa GLMakie.Lines && haskey(plot.attributes, :label) &&
            plot.label[] == "Open-cavity pair", encounter_axis.scene.plots))
        @test encounter_line.color[] == GLMakie.Makie.to_color(:gray)
        @test encounter_line.linestyle[] == Float32[0, 3, 6]
        encounter_legend = only(filter(block -> block isa GLMakie.Legend, encounter_fig.content))
        encounter_swatch = only(filter(entry -> entry.label[] == "Encounter detected", encounter_legend.entrygroups[][1][2]))
        for frame in (1, 2, 3, 2, 1)
            encounter_controls.current_time[] = frame
            @test encounter_line.visible[] == (frame != 3)
            @test active_line.visible[] == (frame == 2)
            @test GLMakie.Makie.to_color(only(encounter_swatch.elements).attributes.linecolor[]) == GLMakie.Makie.to_color(frame != 3 ? :gray : :white)
            @test encounter_swatch.labelcolor[] == (frame != 3 ? :black : :white)
        end
        encounter_output = joinpath(@__DIR__, "output", "smoke", "videos")
        mkpath(encounter_output)
        GLMakie.save(joinpath(encounter_output, "inactive_encounter_test.png"), encounter_fig)
        cp(joinpath(directory, "encounter.mp4"), joinpath(encounter_output, "inactive_encounter_test.mp4"); force=true)
        handoff_states = [vcat([R_EARTH + 1000e3, 0, 0, 0, 7000, 0],
            [R_EARTH + 1000e3, 1000e3, 0, 0, 7000, 0],
            [R_EARTH + 1050e3, target_y, 0, 0, 7000, 0]) for target_y in (0.0, 500e3, 1000e3)]
        handoff_sol = PairHandoffSolution([0.0, 1.0, 2.0], handoff_states)
        handoff_params = merge(params, Dict(:N=>3, :Pmatrix=>zeros(3, 3), :use_los=>false,
            :cavity=>Dict((1, 3)=>copy(params[:cavity][(1, 2)]), (2, 3)=>copy(params[:cavity][(1, 2)]))))
        handoff_fig, handoff_controls = animate_all_satellites_3d_smooth_helper_target(handoff_sol, handoff_params, 2;
            radial_exaggeration=40.0, reference_altitude_km=1000.0,
            animation_fps=3, duration_seconds=1, output_file=joinpath(directory, "handoff.mp4"))
        handoff_axis = GLMakie.content(handoff_fig[1, 1])
        handoff_links = filter(plot -> plot isa GLMakie.Lines && haskey(plot.attributes, :label) &&
            plot.label[] == "Open-cavity pair", handoff_axis.scene.plots)
        handoff_legend = only(filter(block -> block isa GLMakie.Legend, handoff_fig.content))
        cavity_entry = only(filter(entry -> entry.label[] == "Laser active",
            handoff_legend.entrygroups[][1][2]))
        encounter_entry = only(filter(entry -> entry.label[] == "Encounter detected",
            handoff_legend.entrygroups[][1][2]))
        for frame in (1, 2, 3, 2, 1, 3)
            handoff_controls.current_time[] = frame
            @test [link.visible[] for link in handoff_links] == [frame == 1, frame == 3]
            @test cavity_entry.label[] == "Laser active"
            @test GLMakie.Makie.to_color(only(cavity_entry.elements).attributes.linecolor[]) == GLMakie.Makie.to_color(frame == 2 ? :white : :orange)
            @test cavity_entry.labelcolor[] == (frame == 2 ? :white : :black)
            @test encounter_entry.labelcolor[] == (frame == 2 ? :white : :black)
        end
        handoff_output = joinpath(@__DIR__, "output", "smoke", "videos")
        mkpath(handoff_output)
        GLMakie.save(joinpath(handoff_output, "legend_handoff_test.png"), handoff_fig)
        cp(joinpath(directory, "handoff.mp4"), joinpath(handoff_output, "legend_handoff_test.mp4"); force=true)
        off_params = merge(params, Dict(:max_range=>1.0))
        off_fig, off_controls = animate_all_satellites_3d_smooth_helper_target(sol, off_params, 1;
            radial_exaggeration=40.0, reference_altitude_km=1000.0, axis_limit_km=physical_axis_limit(40.0),
            animation_fps=2, duration_seconds=1, output_file=joinpath(directory, "off.mp4"))
        off_axis = GLMakie.content(off_fig[1, 1])
        off_link = only(filter(plot -> plot isa GLMakie.Lines && haskey(plot.attributes, :label) &&
            plot.label[] == "Open-cavity pair", off_axis.scene.plots))
        for frame in 1:2
            off_controls.current_time[] = frame
            @test !off_link.visible[]
            @test off_link.alpha[] == 0
            @test occursin("OFF", GLMakie.content(off_fig[1, 1]).title[])
        end
        @test !any(plot -> haskey(plot.attributes, :label) && endswith(string(plot.label[]), " projection"), off_axis.scene.plots)
        projection_fig, projection_controls = animate_all_satellites_3d_smooth_helper_target(sol, params, 1;
            show_projections=true, helper_trails=false, radial_exaggeration=40.0,
            reference_altitude_km=1000.0, axis_limit_km=physical_axis_limit(40.0),
            radial_ticks_km=[500, 1000, 1010, 1020, 1030, 1040, 1050],
            animation_fps=2, duration_seconds=1, output_file=joinpath(directory, "projections.mp4"))
        projection_axis = GLMakie.content(projection_fig[1, 1])
        projection_legend = only(filter(block -> block isa GLMakie.Legend, projection_fig.content))
        projection_legend_labels = [entry.label[] for entry in projection_legend.entrygroups[][1][2]]
        @test !any(label -> endswith(label, " projection"), projection_legend_labels)
        @test "Encounter detected" in projection_legend_labels
        @test "Laser active" in projection_legend_labels
        main_trails = filter(plot -> plot isa GLMakie.Lines && haskey(plot.attributes, :label) &&
            plot.label[] in ("Helper Sat", "Target Sat"), projection_axis.scene.plots)
        all_markers = filter(plot -> plot isa GLMakie.Scatter, projection_axis.scene.plots)
        @test length(all_markers) == 8
        for frame in 1:2
            projection_controls.current_time[] = frame
            for (plane_index, (plane, fixed_coordinate, plane_position)) in enumerate(
                    (("XY", 3, -0.999 * display_limit), ("XZ", 2, 0.999 * display_limit), ("YZ", 1, 0.999 * display_limit)))
                projected_trails = filter(plot -> plot isa GLMakie.Lines && haskey(plot.attributes, :label) &&
                    plot.label[] == "$plane projection", projection_axis.scene.plots)
                @test length(projected_trails) == 2
                @test [plot.visible[] for plot in projected_trails] == [false, true]
                for satellite in 1:2
                    source_points = main_trails[satellite][1][]
                    projected_points = projected_trails[satellite][1][]
                    @test length(projected_points) == length(source_points)
                    for component in 1:3
                        expected = component == fixed_coordinate ? fill(Float32(plane_position), length(source_points)) :
                            [point[component] for point in source_points]
                        @test [point[component] for point in projected_points] == expected
                    end
                    projected_marker = all_markers[(satellite-1)*4 + 1 + plane_index]
                    @test projected_marker.visible[]
                    @test only(projected_marker[1][]) == last(projected_points)
                end
            end
        end
        @test sol.u == original_states
    end
    @test earth_bounds[1] == earth_bounds[2]
    mktempdir() do directory
        for (helpers_enabled, target_enabled, master_enabled) in
                ((false, true, true), (true, false, true), (false, false, true), (true, true, false))
            fig, controls = animate_all_satellites_3d_smooth_helper_target(sol, params, 1;
                show_projections=master_enabled, helper_projections=helpers_enabled,
                target_projections=target_enabled, helper_trails=true, target_trail=true,
                radial_exaggeration=40.0, reference_altitude_km=1000.0, axis_limit_km=physical_axis_limit(40.0),
                animation_fps=2, duration_seconds=1, output_file=joinpath(directory, "selected.mp4"))
            axis = GLMakie.content(fig[1, 1])
            markers = filter(plot -> plot isa GLMakie.Scatter, axis.scene.plots)
            source_markers = filter(plot -> haskey(plot.attributes, :label) &&
                plot.label[] in ("Helper Sat", "Target Sat"), markers)
            selected_satellites = master_enabled ? findall([helpers_enabled, target_enabled]) : Int[]
            @test length(source_markers) == 2
            @test all(plot -> plot.visible[], source_markers)
            @test length(markers) == 2 + 3 * length(selected_satellites)
            for frame in 1:2
                controls.current_time[] = frame
                for (plane, coordinate, offset) in (("XY", 3, -0.999 * display_limit), ("XZ", 2, 0.999 * display_limit), ("YZ", 1, 0.999 * display_limit))
                    trails = filter(plot -> plot isa GLMakie.Lines && haskey(plot.attributes, :label) &&
                        plot.label[] == "$plane projection", axis.scene.plots)
                    @test length(trails) == length(selected_satellites)
                    for (trail, satellite) in zip(trails, selected_satellites)
                        position = only(source_markers[satellite][1][])
                        expected = GLMakie.Point3f(ntuple(component -> component == coordinate ? offset : position[component], 3))
                        @test trail.visible[]
                        @test last(trail[1][]) == expected
                        @test any(plot -> plot.visible[] && only(plot[1][]) == expected, markers)
                    end
                end
            end
        end
    end
    @test axes[1] == axes[2]
    circular_sol = CircularAnimationSolution([0.0, 6000.0], [copy(first(sol.u)), copy(first(sol.u))])
    output = joinpath(@__DIR__, "output", "smoke", "videos")
    mkpath(output)
    circle_fig, circle_controls = animate_all_satellites_3d_smooth_helper_target(circular_sol, params, 1;
        radial_exaggeration=40.0, reference_altitude_km=1000.0,
        radial_ticks_km=[0, 500, 1000, 1050], helper_trails=true,
        show_projections=true, helper_projections=false, target_projections=true,
        animation_fps=24, duration_seconds=3, output_file=joinpath(output, "helper_referenced_orbit_test.mp4"))
    circle_axis = GLMakie.content(circle_fig[1, 1])
    circle_trails = filter(plot -> plot isa GLMakie.Lines && haskey(plot.attributes, :label) &&
        plot.label[] in ("Helper Sat", "Target Sat"), circle_axis.scene.plots)
    circle_controls.current_time[] = 72
    for satellite in 1:2
        physical_radius = norm(first(sol.u)[6*(satellite-1)+1:6*(satellite-1)+3])
        rendered_points = circle_trails[satellite][1][]
        @test length(rendered_points) == 72
        for (point, angle) in zip(rendered_points, range(0, 2pi; length=72))
            @test norm(point) ≈ expected_radius(physical_radius, 40.0) atol=4
            @test collect(point) ./ norm(point) ≈ [cos(angle), sin(angle)/sqrt(2), sin(angle)/sqrt(2)] atol=2e-7
        end
    end
    for (helper_point, target_point) in zip(circle_trails[1][1][], circle_trails[2][1][])
        @test norm(helper_point) ≈ reference_radius atol=2
        @test norm(target_point) - norm(helper_point) ≈ 40 * 50e3 atol=2
    end
    GLMakie.save(joinpath(output, "helper_referenced_orbit_test.png"), circle_fig)
    @test_throws ArgumentError animate_all_satellites_3d_smooth_helper_target(sol, params, 1; radial_exaggeration=0)
    @test_throws ArgumentError animate_all_satellites_3d_smooth_helper_target(sol, params, 1; radial_exaggeration=NaN)
    @test_throws ArgumentError animate_all_satellites_3d_smooth_helper_target(sol, params, 1;
        radial_exaggeration=40, reference_altitude_km=-1)
    @test_throws ArgumentError animate_all_satellites_3d_smooth_helper_target(sol, params, 1; axis_limit_km=100)
    @test_throws ArgumentError animate_all_satellites_3d_smooth_helper_target(sol, params, 1;
        radial_exaggeration=40, axis_limit_km=7400)
    @test_throws ArgumentError animate_all_satellites_3d_smooth_helper_target(sol, params, 1; radial_ticks_km=[-1])
    @test_throws ArgumentError animate_all_satellites_3d_smooth_helper_target(sol, params, 1; radial_ticks_km=Float64[])
end