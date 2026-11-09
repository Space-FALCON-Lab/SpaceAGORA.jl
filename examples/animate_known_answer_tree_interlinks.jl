module KnownAnswerTreeAnimation

using CSV
using DataFrames
using Printf
import Plots

const OUTPUT = normpath(joinpath(@__DIR__, "..", "output", "known_answer_tree_interlinks"))
const PAIRS = ["1-3", "2-3", "2-4"]
const ALL_PAIRS = ["1-2", "1-3", "1-4", "2-3", "2-4", "3-4"]
const GEOMETRY_COLUMNS = [Symbol("geometrically_feasible_", replace(pair, "-" => "_")) for pair in ALL_PAIRS]
const COLORS = ["#cbd5e1", "#16a34a", "#f97316"]
const EXCLUDED_COLOR = "#a8b3c2"
const NODE_POSITIONS = ((580.0, 170.0), (180.0, 170.0), (180.0, 570.0), (580.0, 570.0))
const LINK_LABEL_POSITIONS = ((380, 140), (490, 300), (625, 370), (135, 370), (270, 300), (380, 600))

function matching_pairs(value::AbstractString)
    value == "empty" && return Set{String}()
    pairs = split(value, ';')
    allunique(pairs) || throw(ArgumentError("Duplicate link in recorded matching: $value"))
    return Set(String.(pairs))
end

function read_recording(directory)
    filenames = ("spacecraft_before_run.csv", "candidates_before_run.csv", "accepted_step_comparison.csv")
    for filename in filenames
        isfile(joinpath(directory, filename)) || error(
            "Missing $filename in $directory. Run examples/known_answer_tree_interlinks.jl first.")
    end
    spacecraft, candidates, samples = map(filenames) do filename
        # Otherwise CSV interprets pair labels such as "1-3" as dates.
        CSV.read(joinpath(directory, filename), DataFrame;
            types=(index, name) -> name in (:pair, :role, :active, :available) ? String :
                name in GEOMETRY_COLUMNS ? Bool : nothing)
    end
    sort!(spacecraft, :index)
    spacecraft.index == collect(1:4) && all(==(1), spacecraft.terminal_capacity) ||
        throw(ArgumentError("This renderer requires the four-spacecraft, one-terminal tree check."))
    allunique(candidates.pair) && Set(candidates.pair) == Set(ALL_PAIRS) &&
        Set(candidates[candidates.allowed, :pair]) == Set(PAIRS) ||
        throw(ArgumentError("Expected all six pairs, with only the tree 1-3-2-4 allowed."))
    all(edge -> edge.pair == "$(edge.a_spacecraft)-$(edge.b_spacecraft)", eachrow(candidates)) ||
        throw(ArgumentError("Pair labels must match their spacecraft endpoints."))
    all(==(1), candidates.a_terminal) && all(==(1), candidates.b_terminal) ||
        throw(ArgumentError("The tree check uses terminal 1 at each endpoint."))
    excluded = only(eachrow(candidates[candidates.pair .== "1-4", :]))
    !excluded.allowed || throw(ArgumentError("Pair 1-4 must remain excluded from the allowed graph."))
    nrow(samples) >= 2 && all(isfinite, samples.time_s) && all(>(0), diff(samples.time_s)) ||
        throw(ArgumentError("Need at least two strictly increasing, finite scheduling times."))
    all(column -> column in propertynames(samples), GEOMETRY_COLUMNS) ||
        throw(ArgumentError("Recording lacks per-pair geometry. Rerun the simulation example " *
            "(known_answer_tree_interlinks.jl or changing_tree_interlinks.jl) before rendering."))
    all(column -> all(value -> value isa Bool, samples[!, column]), GEOMETRY_COLUMNS) ||
        throw(ArgumentError("Every pair needs a Boolean geometric-feasibility value at every sample."))
    allowed = candidates[candidates.allowed, :]
    for sample in eachrow(samples)
        active, available = matching_pairs(sample.active), matching_pairs(sample.available)
        issubset(active, available) && issubset(available, Set(PAIRS)) ||
            throw(ArgumentError("Selected/available links violate the allowed graph at $(sample.time_s) s."))
        feasible_allowed = Set(pair for (pair, column) in zip(ALL_PAIRS, GEOMETRY_COLUMNS)
            if pair in PAIRS && sample[column])
        available == feasible_allowed ||
            throw(ArgumentError("Recorded availability disagrees with geometry at $(sample.time_s) s."))
        use = zeros(Int, 4)
        for edge in eachrow(allowed)
            edge.pair in active || continue
            use[edge.a_spacecraft] += 1
            use[edge.b_spacecraft] += 1
        end
        all(use .<= spacecraft.terminal_capacity) && maximum(use) == sample.maximum_terminal_use ||
            throw(ArgumentError("Recorded terminal use is inconsistent at $(sample.time_s) s."))
        isfinite(sample.selected_total_m_per_s) && isfinite(sample.largest_edge_m_per_s) &&
            isapprox(sample.selected_total_m_per_s, sample.exhaustive_best_m_per_s; rtol=1e-12, atol=1e-14) ||
            throw(ArgumentError("Recorded score does not match the exhaustive check at $(sample.time_s) s."))
        !sample.excluded_1_4_active ||
            throw(ArgumentError("Excluded pair 1-4 was marked active at $(sample.time_s) s."))
        isfinite(sample.excluded_1_4_distance_m) && sample.excluded_1_4_distance_m > 0 ||
            throw(ArgumentError("Invalid excluded-pair distance at $(sample.time_s) s."))
        sample.excluded_1_4_geometrically_feasible == sample.geometrically_feasible_1_4 ==
            (sample.excluded_1_4_distance_m <= excluded.range_m) ||
            throw(ArgumentError("Excluded-pair geometry is inconsistent at $(sample.time_s) s."))
    end
    return (; spacecraft, candidates, samples, excluded_range_m=excluded.range_m)
end

function animation_frames(recording; frame_count::Int=121)
    frame_count >= 2 || throw(ArgumentError("frame_count must be at least two."))
    samples = recording.samples
    times = range(first(samples.time_s), last(samples.time_s); length=frame_count)
    # Hold the last recorded schedule; never interpolate link status or its score.
    indices = searchsortedlast.(Ref(samples.time_s), times)
    return map(times, indices) do time, index
        sample = samples[index, :]
        active, available = matching_pairs(sample.active), matching_pairs(sample.available)
        (; time_s=time, sample_time_s=sample.time_s, sample_index=index,
            status=[!sample[column] ? 0 : (pair in PAIRS ? (pair in active ? 2 : Int(pair in available)) : -1)
                for (pair, column) in zip(ALL_PAIRS, GEOMETRY_COLUMNS)],
            selected_total_m_per_s=sample.selected_total_m_per_s,
            largest_edge_m_per_s=sample.largest_edge_m_per_s,
            excluded_distance_m=sample.excluded_1_4_distance_m,
            excluded_feasible=sample.excluded_1_4_geometrically_feasible,
            maximum_terminal_use=sample.maximum_terminal_use)
    end
end

# Adapted from the ORACLE/forty-spacecraft example's Plots schematic GIF style.
function graph_frame(recording, frame)
    figure = Plots.plot(; size=(1200, 800), xlims=(0, 1200), ylims=(0, 800),
        yflip=true, axis=false, ticks=false, legend=false, grid=false,
        background_color=:white, foreground_color="#0f172a")
    label(x, y, text, size=12, color="#0f172a", align=:center) =
        Plots.annotate!(figure, x, y, Plots.text(text, size, color, align))
    for (pair, code, label_position) in zip(ALL_PAIRS, frame.status, LINK_LABEL_POSITIONS)
        edge = only(eachrow(recording.candidates[recording.candidates.pair .== pair, :]))
        a, b = NODE_POSITIONS[edge.a_spacecraft], NODE_POSITIONS[edge.b_spacecraft]
        color = code == -1 ? EXCLUDED_COLOR : COLORS[code + 1]
        Plots.plot!(figure, [a[1], b[1]], [a[2], b[2]]; color,
            linewidth=code == 2 ? 5 : 3, linestyle=code == 0 ? :dash : :solid)
        label(label_position..., pair, 13, code < 1 ? "#64748b" : color)
    end
    for spacecraft in eachrow(recording.spacecraft)
        x, y = NODE_POSITIONS[spacecraft.index]
        color = spacecraft.role == "helper" ? "#0369a1" : "#4338ca"
        Plots.scatter!(figure, [x], [y]; markersize=22, markerstrokewidth=0, color)
        label(x, y, string(spacecraft.index), 16, :white)
        role_y = spacecraft.index in (1, 2) ? y - 70 : y + 55
        label(x, role_y, uppercase(spacecraft.role), 12, color)
        label(x, role_y + 25, @sprintf("%.0f km | %.0f kg",
            spacecraft.altitude_m / 1000, spacecraft.mass_kg), 10, "#475569")
    end
    legend = [(-1, "Not allowed / feasible"), (1, "Available / idle"), (2, "Active")]
    any(==(0), frame.status) && push!(legend, (0, "Unavailable"))
    for (index, (code, text)) in enumerate(legend)
        y = 285 + 80(index - 1)
        color = code == -1 ? EXCLUDED_COLOR : COLORS[code + 1]
        Plots.plot!(figure, [790, 840], [y, y]; color,
            linewidth=3, linestyle=code == 0 ? :dash : :solid)
        label(865, y, "$text ($(count(==(code), frame.status)))", 13, "#0f172a", :left)
    end
    start_s, end_s = first(recording.samples.time_s), last(recording.samples.time_s)
    progress = (frame.time_s - start_s) / (end_s - start_s)
    cursor_x = 140 + 920progress
    label(600, 697, @sprintf("Simulation time: %.1f / %.1f s", frame.time_s, end_s), 13)
    Plots.plot!(figure, [140, 1060], [738, 738]; color="#e2e8f0", linewidth=5)
    Plots.plot!(figure, [140, cursor_x], [738, 738]; color="#0369a1", linewidth=5)
    Plots.scatter!(figure, [cursor_x], [738]; color="#0369a1", markersize=7, markerstrokewidth=0)
    for (fraction, time) in zip(range(0, 1; length=5), range(start_s, end_s; length=5))
        label(140 + 920fraction, 772, @sprintf("%.0f s", time), 10, "#475569")
    end
    return figure
end

function save_animation(; input=OUTPUT, output=input, frame_count::Int=121, fps::Int=12)
    fps > 0 || throw(ArgumentError("fps must be positive."))
    recording = read_recording(input)
    frames = animation_frames(recording; frame_count)
    mkpath(output)
    gif_path = joinpath(output, "interlinks.gif")
    preview_path = joinpath(output, "interlinks_preview.png")
    withenv("GKSwstype" => "100") do
        mktempdir() do directory
            animation = Plots.Animation(directory, String[])
            for (index, frame) in enumerate(frames)
                figure = graph_frame(recording, frame)
                index == 1 && Plots.savefig(figure, preview_path)
                Plots.frame(animation, figure)
            end
            Plots.gif(animation, gif_path; fps)
        end
    end
    println("Saved ", length(frames), " graph frames to ", gif_path)
    return (; gif_path, preview_path, frame_count=length(frames), fps)
end

end # module KnownAnswerTreeAnimation

if abspath(PROGRAM_FILE) == @__FILE__
    KnownAnswerTreeAnimation.save_animation()
end
