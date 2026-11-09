module BacklogAnimations

using CSV, DataFrames, Dates, FileIO, GeometryBasics, LinearAlgebra, Printf, SHA, Test, TOML
import GLMakie

const ROOT = @__DIR__
const DATA_ROOT = joinpath(ROOT, "results", "backlog")
const REFERENCE_ROOT = normpath(joinpath(ROOT, "..", "..", "Summary_Table_for_Other_Labs"))
const REFERENCE_ANIMATION = joinpath(REFERENCE_ROOT, "3_Kuang's Prototype Code_for_animation",
    "functions", "10_Animation_ver2.jl")
const EARTH_TEXTURE = joinpath(dirname(dirname(REFERENCE_ANIMATION)), "input", "8k_earth_daymap.jpg")
const CASES = ("producer_control", "production_competition", "busy_relay",
    "direct_route_control", "production_competition_refined", "busy_relay_refined")
const FPS = 30
const SIZE = (1600, 1080)
const INPUTS = ("spacecraft_before_run.csv", "settings.toml", "accepted_steps.csv",
    "candidate_scores.csv", "node_accounting.csv", "network_accounting.csv", "link_schedule.csv")
const TITLES = Dict(
    "producer_control" => "Empty producer: continuous data delivery",
    "production_competition" => "Two senders competing for one receiver",
    "busy_relay" => "Busy relay: incoming traffic can displace local production",
    "direct_route_control" => "Direct-route control: relay input disabled",
    "production_competition_refined" => "Competing senders: smaller maximum solver step",
    "busy_relay_refined" => "Busy relay: smaller maximum solver step")

edge_key(row) = (Int(row.tx), Int(row.tx_terminal), Int(row.rx), Int(row.rx_terminal))
edge_label(key) = "$(key[1]).$(key[2])->$(key[3]).$(key[4])"
matching_keys(text) = text == "empty" ? Set{String}() : Set(split(text, ';'))
matching_label(keys) = isempty(keys) ? "empty" : join(sort(collect(keys)), ';')
file_hash(path) = bytes2hex(open(sha256, path))

function write_toml(path, values)
    open(path, "w") do io
        TOML.print(io, values; sorted=true)
    end
end

function load_case(name)
    name in CASES || throw(ArgumentError("Unknown animation case: $name"))
    directory = joinpath(DATA_ROOT, name)
    hashes = Dict(file => file_hash(joinpath(directory, file)) for file in INPUTS)
    settings = TOML.parsefile(joinpath(directory, "settings.toml"))
    spacecraft = CSV.read(joinpath(directory, "spacecraft_before_run.csv"), DataFrame)
    history = CSV.read(joinpath(directory, "accepted_steps.csv"), DataFrame)
    candidates = CSV.read(joinpath(directory, "candidate_scores.csv"), DataFrame)
    nodes = CSV.read(joinpath(directory, "node_accounting.csv"), DataFrame)
    network = CSV.read(joinpath(directory, "network_accounting.csv"), DataFrame)
    schedule = CSV.read(joinpath(directory, "link_schedule.csv"), DataFrame)
    sort!(spacecraft, :spacecraft)
    sort!(nodes, [:time_s, :spacecraft])
    sort!(candidates, [:time_s, :tx, :tx_terminal, :rx, :rx_terminal])
    times = Float64.(history.time_s)
    n, samples = nrow(spacecraft), length(times)
    @test spacecraft.spacecraft == collect(1:n)
    @test n in (2, 3)
    @test samples >= 2 && first(times) == 0.0 && all(diff(times) .> 0)
    @test last(times) == settings["duration_s"]
    @test nrow(nodes) == n * samples
    @test nodes.time_s == repeat(times; inner=n)
    @test nodes.spacecraft == repeat(collect(1:n), samples)
    @test network.time_s == times
    destination = only(findall(spacecraft.destination))
    positions = Array{Float64}(undef, 3, n, samples)
    for (component, column) in enumerate((:x_m, :y_m, :z_m))
        positions[component, :, :] .= reshape(Float64.(nodes[!, column]), n, samples)
    end
    stored = reshape(Float64.(nodes.stored_bits), n, samples)
    @test all(isfinite, positions) && all(isfinite, stored)
    @test all(spacecraft.capacity_bits .> 0)
    @test all(stored .>= -1.0)
    @test all(stored .<= spacecraft.capacity_bits .+ 1.0)
    @test maximum(abs, nodes.balance_residual_bits) <= 1.0
    @test maximum(abs, network.network_residual_bits) <= 1.0
    received = reshape(Float64.(nodes.received_bits), n, samples)
    @test received[destination, :] ≈ network.delivered_bits atol=1e-6

    first_candidates = candidates[candidates.time_s .== first(times), :]
    edges = [(key=edge_key(row), label=edge_label(edge_key(row)), allowed=Bool(row.allowed),
        rate_bps=Float64(row.link_rate_bps)) for row in eachrow(first_candidates)]
    edge_count = length(edges)
    @test edge_count in 1:3
    @test candidates.time_s == repeat(times; inner=edge_count)
    @test edge_key.(eachrow(candidates)) == repeat(getproperty.(edges, :key), samples)
    @test candidates.allowed == repeat(getproperty.(edges, :allowed), samples)
    available = reshape(Bool.(candidates.available), edge_count, samples)
    selected = reshape(Bool.(candidates.selected), edge_count, samples)
    scores = reshape(Float64.(candidates.score_bps), edge_count, samples)
    @test all(.!selected .| available)
    @test all(isfinite, scores)
    for k in 1:samples
        expected = Set(edges[e].label for e in 1:edge_count if selected[e, k])
        @test expected == matching_keys(history.selected_matching[k])
        @test length(expected) == history.selected_count[k]
        @test sum(available[:, k]) == history.feasible_count[k]
        @test isapprox(sum(scores[selected[:, k], k]), history.selected_score_bps[k]; atol=1e-5, rtol=1e-12)
        endpoints = [(satellite, terminal) for e in 1:edge_count if selected[e, k]
            for (satellite, terminal) in ((edges[e].key[1], edges[e].key[2]), (edges[e].key[3], edges[e].key[4]))]
        @test allunique(endpoints)
        if k < samples
            for (kind, flags) in (("active", selected), ("contact", available))
                recorded = Set("$(row.a_spacecraft).$(row.a_terminal)->$(row.b_spacecraft).$(row.b_terminal)"
                    for row in eachrow(schedule) if row.window == kind && row.start_s <= times[k] < row.end_s)
                @test recorded == Set(edges[e].label for e in 1:edge_count if flags[e, k])
            end
        end
    end
    changes = [k for k in 2:samples if selected[:, k] != selected[:, k - 1]]
    @test length(changes) == count(history.changed)

    # Recover the planet radius from the recorded geometric clearance, not another simulation.
    edge = first(edges).key
    a, b = positions[:, edge[1], 1], positions[:, edge[3], 1]
    displacement = b - a
    closest = a + clamp(-dot(a, displacement) / dot(displacement, displacement), 0.0, 1.0) * displacement
    earth_radius_m = norm(closest) - first_candidates.planet_clearance_m[1]
    @test isapprox(earth_radius_m, 6_378_136.6; atol=1e-6, rtol=0)
    @test all(i -> isapprox(norm(positions[:, i, 1]) - earth_radius_m, 700e3; atol=1e-6, rtol=0), 1:n)
    return (; name, directory, hashes, settings, spacecraft, times, positions, stored,
        network, edges, available, selected, scores, changes, destination, earth_radius_m)
end

function playback_frames(bundle)
    times = bundle.times
    change_times = Set(times[bundle.changes])
    regular = range(first(times), last(times); length=round(Int, (last(times) - first(times)) * FPS) + 1)
    display_times = sort!(unique(vcat(collect(regular), collect(change_times))))
    frames = NamedTuple[]
    for t in display_times
        left = searchsortedlast(times, t)
        right = min(left + 1, length(times))
        fraction = left == right ? 0.0 : (t - times[left]) / (times[right] - times[left])
        reason = t == first(times) ? "INITIAL SELECTION" :
            t == last(times) ? "FINAL TOTALS" : t in change_times ? "SELECTION CHANGE" : ""
        repeat_count = t == first(times) ? round(Int, 1.5FPS) :
            t == last(times) ? 2FPS : t in change_times ? FPS ÷ 2 : 1
        frame = (; time_s=Float64(t), left, right, fraction, reason)
        append!(frames, fill(frame, repeat_count))
    end
    @test first(frames).time_s == first(times) && last(frames).time_s == last(times)
    @test issorted(getproperty.(frames, :time_s))
    @test all(t -> any(frame -> frame.time_s == t, frames), change_times)
    return frames
end

interpolate(values::AbstractVector, frame) =
    (1 - frame.fraction) * values[frame.left] + frame.fraction * values[frame.right]
interpolate(values::AbstractMatrix, frame) =
    (1 - frame.fraction) .* values[:, frame.left] .+ frame.fraction .* values[:, frame.right]
positions_at(bundle, frame) =
    (1 - frame.fraction) .* bundle.positions[:, :, frame.left] .+
    frame.fraction .* bundle.positions[:, :, frame.right]

function write_frame_map(bundle, frames, output)
    table = DataFrame(map(enumerate(frames)) do (index, frame)
        selected = Set(edge.label for (e, edge) in enumerate(bundle.edges) if bundle.selected[e, frame.left])
        available = Set(edge.label for (e, edge) in enumerate(bundle.edges) if bundle.available[e, frame.left])
        (; video_frame=index, video_time_s=(index - 1) / FPS, simulation_time_s=frame.time_s,
            recorded_sample=frame.left, recorded_sample_time_s=bundle.times[frame.left],
            interpolation_fraction=frame.fraction, review_hold=frame.reason,
            selected_matching=matching_label(selected), available_matching=matching_label(available),
            displayed_delivered_mbit=interpolate(bundle.network.delivered_bits, frame) / 1e6,
            displayed_production_mbit=interpolate(bundle.network.produced_bits, frame) / 1e6,
            displayed_overflow_mbit=interpolate(bundle.network.lost_bits, frame) / 1e6)
    end)
    CSV.write(joinpath(output, "frames.csv"), table)
    @test nrow(CSV.read(joinpath(output, "frames.csv"), DataFrame)) == length(frames)
    return table
end

function node_roles(bundle)
    n = nrow(bundle.spacecraft)
    roles = fill("Sender", n)
    roles[1] = "Producer"
    if occursin("relay", bundle.name) || bundle.name == "direct_route_control"
        roles[1] = "Source"
        roles[2] = "Relay + producer"
    elseif occursin("competition", bundle.name)
        roles[2] = "Stored-data sender"
    end
    roles[bundle.destination] = "Destination"
    return roles
end

function network_positions(bundle)
    nrow(bundle.spacecraft) == 2 && return Point2f[(-0.9, 0.3), (0.9, 0.3)]
    occursin("competition", bundle.name) &&
        return Point2f[(-0.95, 0.85), (-0.95, -0.25), (0.95, 0.3)]
    return Point2f[(-1.0, -0.15), (0.0, 0.9), (1.0, -0.15)]
end

function create_scene(bundle, frames, texture)
    n = nrow(bundle.spacecraft)
    current = GLMakie.Observable(1)
    state = GLMakie.lift(index -> frames[index], current)
    roles = node_roles(bundle)
    colors = [:royalblue, :purple, :seagreen][1:n]
    colors[bundle.destination] = :seagreen
    fig = GLMakie.Figure(size=SIZE, backgroundcolor=:white, fontsize=17)
    GLMakie.Label(fig[1, 1:2], TITLES[bundle.name]; fontsize=25, font=:bold, tellwidth=false)
    GLMakie.Label(fig[2, 1:2],
        "$(bundle.name) | 20 s saved simulation | dtmax = $(bundle.settings["dtmax_s"]) s | DataBacklogPolicy";
        fontsize=17, tellwidth=false)
    clock = GLMakie.lift(state) do frame
        @sprintf("Simulation t = %.6f s  |  latest scheduling sample t = %.6f s  |  %d / %d set changes",
            frame.time_s, bundle.times[frame.left], count(k -> k <= frame.left, bundle.changes), length(bundle.changes))
    end
    GLMakie.Label(fig[3, 1:2], clock; fontsize=19, tellwidth=false)

    orbit = GLMakie.Axis3(fig[4, 1]; title="Orbital context: inertial coordinates (km)",
        xlabel="X (km)", ylabel="Y (km)", zlabel="Z (km)", aspect=(1, 1, 1),
        azimuth=0.35, elevation=0.65, perspectiveness=0.4,
        xticks=[-8000, 0, 8000], yticks=[-8000, 0, 8000], zticks=[-8000, 0, 8000])
    radius = bundle.earth_radius_m / 1e3
    sphere = uv_normal_mesh(Tesselation(Sphere(Point3f(0, 0, 0), Float32(radius)), 48))
    GLMakie.mesh!(orbit, sphere; color=texture)
    extent = 1.13 * maximum(norm(bundle.positions[:, i, 1]) for i in 1:n) / 1e3
    GLMakie.limits!(orbit, (-extent, extent), (-extent, extent), (-extent, extent))
    orbital_radius = norm(bundle.positions[:, 1, 1]) / 1e3
    angles = range(0.0, 2pi; length=361)
    GLMakie.lines!(orbit, [Point3f(orbital_radius * cos(t), orbital_radius * sin(t), 0) for t in angles];
        color=(:gray, 0.5), linestyle=:dash, linewidth=1)
    orbital_points = GLMakie.lift(state) do frame
        positions = positions_at(bundle, frame) / 1e3
        [Point3f(positions[:, i]) for i in 1:n]
    end
    orbit_markers = GLMakie.scatter!(orbit, orbital_points; color=colors, markersize=13)
    for i in 1:n
        trail = GLMakie.lift(state) do frame
            points = [Point3f(bundle.positions[:, i, k] / 1e3) for k in 1:frame.left]
            push!(points, Point3f(positions_at(bundle, frame)[:, i] / 1e3))
        end
        GLMakie.lines!(orbit, trail; color=colors[i], linewidth=2)
    end

    diagram = GLMakie.Axis(fig[4, 2]; title="Directed links - enlarged schematic, not orbital geometry",
        limits=(-1.65, 1.65, -1.0, 1.6))
    GLMakie.hidedecorations!(diagram)
    GLMakie.hidespines!(diagram)
    points = network_positions(bundle)
    active_plots, contact_plots, excluded_plots, edge_labels = Any[], Any[], Any[], Any[]
    for (e, edge) in enumerate(bundle.edges)
        a, b = points[edge.key[1]], points[edge.key[3]]
        direction = normalize(b - a)
        normal = Point2f(-direction[2], direction[1])
        endpoints = [a + 0.12f0 * direction, b - 0.12f0 * direction]
        selected = GLMakie.lift(frame -> bundle.selected[e, frame.left], state)
        idle = GLMakie.lift(frame -> edge.allowed && bundle.available[e, frame.left] &&
            !bundle.selected[e, frame.left], state)
        unavailable = GLMakie.lift(frame -> edge.allowed && !bundle.available[e, frame.left], state)
        push!(contact_plots, GLMakie.lines!(diagram, endpoints;
            color=:gray, linestyle=:dash, linewidth=2, visible=idle))
        GLMakie.lines!(diagram, endpoints; color=:lightgray, linestyle=:dot, linewidth=1, visible=unavailable)
        push!(excluded_plots, GLMakie.lines!(diagram, endpoints;
            color=:lightgray, linestyle=:dot, linewidth=2, visible=!edge.allowed))
        push!(active_plots, GLMakie.lines!(diagram, endpoints;
            color=:darkorange, linewidth=5, visible=selected))
        tip = b - 0.16f0 * direction
        triangle = [tip, tip - 0.11f0 * direction + 0.05f0 * normal,
            tip - 0.11f0 * direction - 0.05f0 * normal]
        arrow_color = GLMakie.lift(frame -> bundle.selected[e, frame.left] ? :darkorange :
            edge.allowed && bundle.available[e, frame.left] ? :gray : :lightgray, state)
        GLMakie.poly!(diagram, triangle; color=arrow_color)
        label = GLMakie.lift(state) do frame
            edge.allowed || return "$(edge.label)\nEXCLUDED"
            @sprintf("%s | max %.1f Mbit/s\nscore %.3f Mbit/s",
                edge.label, edge.rate_bps / 1e6, bundle.scores[e, frame.left] / 1e6)
        end
        push!(edge_labels, label)
        GLMakie.text!(diagram, label; position=(a + b) / 2 + 0.26f0 * normal,
            align=(:center, :center), fontsize=13, color=edge.allowed ? :black : :gray)
    end
    GLMakie.scatter!(diagram, points; color=colors, markersize=47)
    for i in 1:n
        upper_label = n == 3 && !occursin("competition", bundle.name) && i == 2
        GLMakie.text!(diagram, string(i); position=points[i], align=(:center, :center),
            color=:white, fontsize=19, font=:bold)
        GLMakie.text!(diagram, roles[i]; position=points[i] + Point2f(0, upper_label ? 0.27 : -0.27),
            align=(:center, :center), fontsize=16, font=:bold)
        text = @sprintf("production %.1f Mbit/s\n%d terminal%s",
            bundle.spacecraft.production_bps[i] / 1e6, bundle.spacecraft.terminal_capacity[i],
            bundle.spacecraft.terminal_capacity[i] == 1 ? "" : "s")
        GLMakie.text!(diagram, text; position=points[i] + Point2f(0, upper_label ? 0.52 : -0.52),
            align=(:center, :center), fontsize=14)
    end

    storage_axis = GLMakie.Axis(fig[5, 1]; title="Onboard storage", ylabel="Capacity used (%)",
        xticks=(1:n, ["SC $i\n$(roles[i])" for i in 1:n]), limits=(0.4, n + 0.6, 0, 118))
    fill_percent = GLMakie.lift(frame -> 100 .* interpolate(bundle.stored, frame) ./ bundle.spacecraft.capacity_bits, state)
    bars = GLMakie.barplot!(storage_axis, collect(1:n), fill_percent; color=colors, width=0.55)
    for i in 1:n
        label = GLMakie.lift(frame -> @sprintf("%.3f / %.0f Mbit",
            interpolate(bundle.stored, frame)[i] / 1e6, bundle.spacecraft.capacity_bits[i] / 1e6), state)
        label_position = GLMakie.lift(values -> Point2f(i, values[i] + 5), fill_percent)
        GLMakie.text!(storage_axis, label; position=label_position, align=(:center, :center), fontsize=14)
    end

    max_data = maximum(vcat(bundle.network.delivered_bits, bundle.network.produced_bits, bundle.network.lost_bits)) / 1e6
    data_axis = GLMakie.Axis(fig[5, 2]; title="Delivery, new production and discarded production",
        xlabel="Simulation time (s)", ylabel="Mbit", limits=(0, last(bundle.times), 0, 1.12max_data))
    cumulative = Dict{Symbol, Any}()
    for (column, label, color, style, width) in ((:delivered_bits, "Delivered", :seagreen, nothing, 4),
        (:produced_bits, "Produced", :royalblue, :dash, 3), (:lost_bits, "Overflow", :firebrick, :dot, 3))
        history_points = GLMakie.lift(state) do frame
            points = [Point2f(bundle.times[k], bundle.network[k, column] / 1e6) for k in 1:frame.left]
            push!(points, Point2f(frame.time_s, interpolate(bundle.network[!, column], frame) / 1e6))
        end
        cumulative[column] = GLMakie.lines!(data_axis, history_points; color, label, linestyle=style, linewidth=width)
    end
    GLMakie.axislegend(data_axis; position=:lt, orientation=:horizontal, labelsize=13)
    cursor = GLMakie.lift(frame -> [frame.time_s], state)
    GLMakie.vlines!(data_axis, cursor; color=(:black, 0.35), linestyle=:dash)
    totals = GLMakie.lift(state) do frame
        @sprintf("DELIVERED %.3f Mbit    |    NEW PRODUCTION %.3f Mbit    |    OVERFLOW %.3f Mbit",
            interpolate(bundle.network.delivered_bits, frame) / 1e6,
            interpolate(bundle.network.produced_bits, frame) / 1e6,
            interpolate(bundle.network.lost_bits, frame) / 1e6)
    end
    GLMakie.Label(fig[6, 1:2], totals; fontsize=21, font=:bold, tellwidth=false)
    GLMakie.Label(fig[7, 1:2],
        "Solid orange: SELECTED   |   Dashed gray: AVAILABLE, IDLE   |   Dotted light gray: EXCLUDED / UNAVAILABLE";
        fontsize=16, tellwidth=false)
    pacing = GLMakie.lift(frame -> isempty(frame.reason) ?
        "Recorded-data replay; positions and counters are interpolated, selections are not." :
        "REVIEW HOLD: $(frame.reason) - simulation time is paused for readability.", state)
    GLMakie.Label(fig[8, 1:2], pacing; fontsize=16, tellwidth=false)
    GLMakie.Label(fig[9, 1:2],
        "Destination is a satellite, not a ground station. Delivered data can include initial backlog. No laser forces are shown.";
        fontsize=14, tellwidth=false)
    GLMakie.rowsize!(fig.layout, 4, 350)
    GLMakie.rowsize!(fig.layout, 5, 250)
    GLMakie.colsize!(fig.layout, 1, GLMakie.Relative(0.44))
    GLMakie.colgap!(fig.layout, 32)
    return (; fig, current, orbit_markers, active_plots, contact_plots, excluded_plots,
        edge_labels, bars, fill_percent, cumulative)
end

function check_rendered_state(scene, bundle, frame)
    expected_positions = positions_at(bundle, frame) / 1e3
    actual_positions = scene.orbit_markers[1][]
    @test length(actual_positions) == nrow(bundle.spacecraft)
    for i in eachindex(actual_positions)
        @test collect(actual_positions[i]) ≈ expected_positions[:, i] atol=1e-3 rtol=0
    end
    for (e, edge) in enumerate(bundle.edges)
        @test scene.active_plots[e].visible[] == bundle.selected[e, frame.left]
        @test scene.contact_plots[e].visible[] ==
            (edge.allowed && bundle.available[e, frame.left] && !bundle.selected[e, frame.left])
        @test scene.excluded_plots[e].visible[] == !edge.allowed
        @test occursin(edge.label, scene.edge_labels[e][])
    end
    @test scene.fill_percent[] ≈ 100 .* interpolate(bundle.stored, frame) ./ bundle.spacecraft.capacity_bits atol=1e-10
    for column in (:delivered_bits, :produced_bits, :lost_bits)
        displayed = last(scene.cumulative[column][1][])
        @test isapprox(displayed[1], frame.time_s; atol=2e-6, rtol=0)
        @test isapprox(displayed[2], interpolate(bundle.network[!, column], frame) / 1e6; atol=4e-6, rtol=0)
    end
end

function verify_video(path, expected_frames)
    ffmpeg = GLMakie.Makie.FFMPEG_jll.ffmpeg()
    ffprobe = GLMakie.Makie.FFMPEG_jll.ffprobe()
    fields = "codec_name,pix_fmt,width,height,r_frame_rate,nb_frames,duration"
    output = read(`$ffprobe -v error -select_streams v:0 -show_entries stream=$fields -of default=noprint_wrappers=1 $path`, String)
    info = Dict{String, String}()
    for line in split(strip(output), '\n')
        key, value = split(line, '='; limit=2)
        info[key] = value
    end
    @test info["codec_name"] == "h264"
    @test info["pix_fmt"] == "yuv420p"
    @test parse(Int, info["width"]) == SIZE[1]
    @test parse(Int, info["height"]) == SIZE[2]
    @test parse(Int, info["nb_frames"]) == expected_frames
    rate = parse.(Float64, split(info["r_frame_rate"], '/'))
    @test rate[1] / rate[2] == FPS
    @test isapprox(parse(Float64, info["duration"]), expected_frames / FPS; atol=1e-5, rtol=0)
    run(`$ffmpeg -v error -xerror -i $path -f null -`)
    @test filesize(path) > 0
    return info
end

function render_case(name, texture; record_video=true)
    result = nothing
    @testset failfast=true "$name animation" begin
        bundle = load_case(name)
        frames = playback_frames(bundle)
        output = joinpath(bundle.directory, "animation")
        mkpath(output)
        write_frame_map(bundle, frames, output)
        scene = create_scene(bundle, frames, texture)
        try
            for index in unique(vcat(1, [findfirst(f -> f.time_s == bundle.times[k], frames) for k in bundle.changes], length(frames)))
                scene.current[] = index
                check_rendered_state(scene, bundle, frames[index])
            end
            scene.current[] = 1
            GLMakie.save(joinpath(output, "preview.png"), scene.fig)
            for block in scene.fig.content
                block isa GLMakie.Label || continue
                bounds = block.layoutobservables.computedbbox[]
                @test all(bounds.origin .>= 0)
                @test all(bounds.origin .+ bounds.widths .<= collect(SIZE))
            end
            path = joinpath(output, "animation.mp4")
            println("Rendering $name: $(length(frames)) frames at $FPS fps.")
            if record_video
                GLMakie.record(scene.fig, path, eachindex(frames); framerate=FPS, visible=false) do index
                    scene.current[] = index
                    check_rendered_state(scene, bundle, frames[index])
                end
                video_info = verify_video(path, length(frames))
                scene.current[] = length(frames)
                GLMakie.save(joinpath(output, "final.png"), scene.fig)
                @test all(file_hash(joinpath(bundle.directory, file)) == hash for (file, hash) in bundle.hashes)
                metadata = Dict("case" => name, "status" => "passed", "generated_utc" => string(now(UTC)),
                    "fps" => FPS, "frames" => length(frames), "video_duration_s" => length(frames) / FPS,
                    "simulation_duration_s" => last(bundle.times), "recorded_samples" => length(bundle.times),
                    "selection_changes" => length(bundle.changes), "input_sha256" => bundle.hashes,
                    "renderer_sha256" => file_hash(@__FILE__), "reference_renderer_sha256" => file_hash(REFERENCE_ANIMATION),
                    "reference_renderer" => REFERENCE_ANIMATION, "active_project" => Base.active_project(),
                    "glmakie_version" => string(pkgversion(GLMakie)), "video_probe" => video_info,
                    "video_sha256" => file_hash(path),
                    "display_interpolation" => "linear positions/counters; left-held recorded selections/scores",
                    "review_holds" => "1.5 s at start; 0.5 s at every selection change; 2 s at final totals")
                write_toml(joinpath(output, "validation.toml"), metadata)
                result = (; case=name, video_path=path, frames=length(frames), fps=FPS,
                    video_duration_s=length(frames) / FPS, simulation_duration_s=last(bundle.times),
                    recorded_samples=length(bundle.times), selection_changes=length(bundle.changes),
                    delivered_mbit=last(bundle.network.delivered_bits) / 1e6,
                    overflow_mbit=last(bundle.network.lost_bits) / 1e6,
                    video_bytes=filesize(path), video_sha256=file_hash(path))
            else
                scene.current[] = length(frames)
                GLMakie.save(joinpath(output, "final.png"), scene.fig)
            end
        finally
            GLMakie.closeall()
        end
    end
    return result
end

function write_index(results)
    CSV.write(joinpath(DATA_ROOT, "animation_index.csv"), DataFrame(results))
    open(joinpath(ROOT, "ANIMATIONS.md"), "w") do io
        println(io, "# Animations of the six backlog cases\n\n",
            "[Open the video gallery](animations.html). These are the original **satellite-to-satellite** ",
            "experiments, not the later forty-spacecraft ground-station run. No simulation was rerun.\n")
        println(io, "| Case | MP4 | Preview | Playback (s) | Recorded set changes | Delivered (Mbit) | Overflow (Mbit) |")
        println(io, "|---|---|---|---:|---:|---:|---:|")
        for result in results
            base = "results/backlog/$(result.case)/animation"
            @printf(io, "| %s | [Video](%s/animation.mp4) | [Initial](%s/preview.png) / [Final](%s/final.png) | %.3f | %d | %.6f | %.6f |\n",
                result.case, base, base, base, result.video_duration_s,
                result.selection_changes, result.delivered_mbit, result.overflow_mbit)
        end
        println(io, "\n## Method and interpretation\n\n",
            "Adapted from [the reference GLMakie renderer](<",
            relpath(REFERENCE_ANIMATION, ROOT), ">) and its recorded-data replay approach: ",
            "reactive `Observable` frame updates, `Axis3` orbital context, `GLMakie.record` MP4 export, ",
            "PNG snapshots, and headless rendering under `xvfb-run`. The existing Earth texture ",
            "and rendering environment are reused read-only. The reference renderer itself is not modified.\n\n",
            "The original renderer is laser/helper/target-specific. This adaptation uses directed **data** ",
            "links, actual terminal IDs, storage bars, and delivery/production/overflow counters. ",
            "The right-hand network is explicitly schematic because the real spacecraft are only ",
            "10-20 km apart and overlap visually in a full Earth view. The dashed orbital circle ",
            "is a reference orbit; colored trails use the saved 20-second trajectories.\n\n",
            "- Every original selected-set change is represented, including the busy relay's ",
            "initial direct link lasting only about 0.000121 s.\n",
            "- A labelled review hold pauses simulation time for 1.5 s initially, 0.5 s at each ",
            "selection change, and 2 s at final totals. These are **presentation pauses**, not longer ",
            "physical contacts. Playback duration therefore differs between cases.\n",
            "- The base display grid is 30 frames per simulated second, with exact change times inserted. ",
            "Positions and cumulative counters use linear interpolation of saved values. The CSVs ",
            "do not contain velocities, so the reference's cubic Hermite interpolation is not claimed. ",
            "Selections, availability and scores are held from the latest accepted scheduler sample; ",
            "they are never interpolated or recomputed.\n",
            "- Orange means selected; dashed gray means available but idle; dotted light gray ",
            "means excluded or unavailable. Arrows indicate direction, not individual packets or forces. ",
            "The label's maximum rate and backlog score are **not** claims of actual carried rate.\n",
            "- All numerical source CSVs remain unchanged. The per-video `frames.csv` maps media frames ",
            "to simulation time, source samples, selections and displayed counters. ",
            "`validation.toml` includes input/code/video hashes and codec checks.\n\n",
            "## Validation\n\nInput states are checked against the accepted-step history and link windows. ",
            "Every encoded frame checks the rendered positions, connection visibility and data curves. ",
            "FFprobe verifies H.264, yuv420p, 1600x1080, 30 fps, frame count and duration; FFmpeg ",
            "decodes each complete video with errors treated as failures. ",
            "See [animation_index.csv](results/backlog/animation_index.csv) and each case's validation record.\n\n",
            "## Reproduce\n\nFrom the repository root, using the existing reference rendering environment:\n\n```sh\n",
            "xvfb-run -a julia --startup-file=no \\\n",
            "  --project=\"../Summary_Table_for_Other_Labs/2_SpaceAGORA.jl\" \\\n",
            "  kuang_later_tasks/animate_backlog.jl\n```\n\n",
            "Append one or more of the six case names to render only those cases. A full run refreshes ",
            "this index and the gallery; partial runs update only the requested case outputs.\n")
    end
    open(joinpath(ROOT, "animations.html"), "w") do io
        println(io, """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Six backlog-scheduling animations</title>
<style>body{font:16px system-ui;margin:32px auto;padding:0 20px;max-width:1500px;background:#f3f6fa;color:#17283b}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(450px,1fr));gap:24px}article{background:white;padding:18px;border-radius:12px}
video{width:100%;height:auto}h2{font-size:20px}p{line-height:1.5}a{color:#185ca3}@media(max-width:520px){.grid{grid-template-columns:1fr}}</style>
</head><body><h1>Six recorded backlog-scheduling experiments</h1>
<p>Satellite-to-satellite data transfer, not ground-station downlinks. Videos use the reference GLMakie animation method.
Review holds pause the simulation clock so brief selections can be inspected.</p><div class="grid">""")
        for result in results
            base = "results/backlog/$(result.case)/animation"
            println(io, "<article><h2>$(TITLES[result.case])</h2><p><code>$(result.case)</code></p>",
                "<video controls preload=\"metadata\" poster=\"$base/preview.png\" data-case=\"$(result.case)\">",
                "<source src=\"$base/animation.mp4\" type=\"video/mp4\">",
                "Your browser cannot play this video. <a href=\"$base/animation.mp4\">Download MP4</a>.</video>",
                "<p>Delivered: $(round(result.delivered_mbit; digits=6)) Mbit | ",
                "Overflow: $(round(result.overflow_mbit; digits=6)) Mbit | Set changes: $(result.selection_changes)</p>",
                "<p><a href=\"$base/animation.mp4\">MP4</a> | <a href=\"$base/final.png\">Final frame</a> | ",
                "<a href=\"$base/frames.csv\">Frame map</a></p></article>")
        end
        println(io, "</div><p><a href=\"ANIMATIONS.md\">Methods and verification</a> | ",
            "<a href=\"REPORT.md\">Original numerical report</a></p></body></html>")
    end
end

function main(names=collect(CASES))
    !isempty(names) && all(name -> name in CASES, names) && allunique(names) ||
        throw(ArgumentError("Pass unique case names from: $(join(CASES, ", "))"))
    GLMakie.activate!()
    texture = FileIO.load(EARTH_TEXTURE)[1:4:end, 1:4:end]
    status_path = joinpath(DATA_ROOT, "animations_status.toml")
    status = Dict{String, Any}("status" => "running", "started_utc" => string(now(UTC)), "cases" => names)
    write_toml(status_path, status)
    passed = false
    try
        results = map(name -> render_case(name, texture), names)
        if Set(names) == Set(CASES)
            write_index(results)
        end
        passed = true
    finally
        status["status"] = passed ? "passed" : "failed"
        status["completed_utc"] = string(now(UTC))
        write_toml(status_path, status)
    end
end

end # module BacklogAnimations

if abspath(PROGRAM_FILE) == @__FILE__
    BacklogAnimations.main(isempty(ARGS) ? collect(BacklogAnimations.CASES) : ARGS)
end
