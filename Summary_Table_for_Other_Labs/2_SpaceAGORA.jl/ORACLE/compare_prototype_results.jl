using Arrow, CSV, DataFrames, LinearAlgebra

function compare_prototype_results(prototype_dir, spaceagora_dir, output_dir)
    prototype = DataFrame(Arrow.Table(joinpath(prototype_dir, "trajectory.feather")))
    spaceagora = DataFrame(Arrow.Table(joinpath(spaceagora_dir, "trajectory.feather")))
    prototype.time == spaceagora.time || error("Trajectory timestamps must match exactly")
    names(prototype) == names(spaceagora) || error("Trajectory schemas must match")
    spacecraft_ids = [parse(Int, matched.captures[1]) for name in names(prototype)
        for matched in (match(r"^sc(\d+)_pos_1$", name),) if matched !== nothing]
    states = DataFrame(spacecraft=Int[], quantity=String[], initial=Float64[],
        final=Float64[], maximum=Float64[], rms=Float64[])
    for spacecraft in spacecraft_ids, (field, quantity) in (("pos", "position_difference_m"), ("vel", "velocity_difference_mps"))
        columns = ["sc$(spacecraft)_$(field)_$(component)" for component in 1:3]
        differences = Matrix(spaceagora[:, columns]) - Matrix(prototype[:, columns])
        magnitudes = [norm(row) for row in eachrow(differences)]
        push!(states, (spacecraft, quantity, first(magnitudes), last(magnitudes),
            maximum(magnitudes), sqrt(sum(abs2, magnitudes)/length(magnitudes))))
    end
    intervals = DataFrame(kind=String[], sc_a=Int[], sc_b=Int[], occurrence=Int[],
        prototype_start_s=Float64[], spaceagora_start_s=Float64[], start_difference_s=Float64[],
        prototype_end_s=Float64[], spaceagora_end_s=Float64[], end_difference_s=Float64[],
        prototype_duration_s=Float64[], spaceagora_duration_s=Float64[], duration_difference_s=Float64[])
    for kind in ("geometry_encounters", "laser_on")
        baseline = DataFrame(Arrow.Table(joinpath(prototype_dir, kind*".feather")))
        candidate = DataFrame(Arrow.Table(joinpath(spaceagora_dir, kind*".feather")))
        pairs = unique(vcat(collect(zip(baseline.sc_a, baseline.sc_b)), collect(zip(candidate.sc_a, candidate.sc_b))))
        for pair in pairs
            left = sort(filter(row -> (row.sc_a, row.sc_b) == pair, baseline), :start_time_s)
            right = sort(filter(row -> (row.sc_a, row.sc_b) == pair, candidate), :start_time_s)
            nrow(left) == nrow(right) || error("Unmatched interval counts for $kind pair $pair")
            for occurrence in 1:nrow(left)
                original, current = left[occurrence, :], right[occurrence, :]
                (original.start_clipped, original.end_clipped) == (current.start_clipped, current.end_clipped) ||
                    error("Interval clipping differs for $kind pair $pair")
                push!(intervals, (kind, pair..., occurrence,
                    original.start_time_s, current.start_time_s, current.start_time_s-original.start_time_s,
                    original.end_time_s, current.end_time_s, current.end_time_s-original.end_time_s,
                    original.duration_s, current.duration_s, current.duration_s-original.duration_s))
            end
        end
    end
    diagnostics = DataFrame(quantity=String[], prototype=Float64[], spaceagora=Float64[], difference=Float64[])
    for axis in ("r", "t", "n")
        column = "dv_$(axis)_accumulated"
        push!(diagnostics, (column, prototype[end, column], spaceagora[end, column],
            spaceagora[end, column]-prototype[end, column]))
    end
    mkpath(output_dir)
    CSV.write(joinpath(output_dir, "state_differences.csv"), states)
    CSV.write(joinpath(output_dir, "interval_differences.csv"), intervals)
    CSV.write(joinpath(output_dir, "diagnostic_differences.csv"), diagnostics)
    println("Matched $(nrow(prototype)) timestamps; $(length(spacecraft_ids)) spacecraft")
    for table in (states, intervals, diagnostics)
        show(stdout, MIME("text/plain"), table; allrows=true, allcols=true)
        println()
    end
    return (; states, intervals, diagnostics)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 3 || error("Usage: compare_prototype_results.jl PROTOTYPE_FEATHER_DIR SPACEAGORA_FEATHER_DIR COMPARISON_DIR")
    compare_prototype_results(ARGS...)
end