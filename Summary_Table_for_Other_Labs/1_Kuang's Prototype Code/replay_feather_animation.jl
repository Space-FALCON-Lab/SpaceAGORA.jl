# xvfb-run -a julia --startup-file=no --project=2_SpaceAGORA.jl "1_Kuang's Prototype Code/replay_feather_animation.jl"

module FeatherAnimationReplay

###########
# Options #
###########
# Edit this folder path, not an individual Feather filename. It must contain
# trajectory.feather, geometry_encounters.feather, and laser_on.feather.
# Use an absolute path, or joinpath(@__DIR__, ...) for a path relative to this script.
const OPTIONS = (
    bundle_directory = joinpath(@__DIR__, "output", "feather",
        # "prototype_N20_h1000km_t1000km_ih0_it5_e0_nu0_T315355.97s_R200km_P10000W_B100_m227kg_gve_inc_J2true_dt10s"),
        # "prototype_N20_h1000km_t1000km_ih0_it5_e0_nu0_T315355.97s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s"),
        "prototype_N20_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s"),
        # "prototype_N20_h1000km_t1150km_ih0_it0_e0_nu0_T325021.624s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s"),
        # "prototype_N100_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s"),

    duration_seconds = 500.0,
    animation_fps    = 30,
)

##########################
# Dependencies and units #
##########################
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "2_SpaceAGORA.jl"); io=devnull)
using Arrow, DataFrames, LinearAlgebra, StaticArrays, GeometryBasics, FileIO, Printf
if Sys.islinux() && isempty(get(ENV, "DISPLAY", ""))
    let executable = Sys.which("Xvfb")
        executable === nothing && error("Headless replay requires Xvfb. Install xvfb or start Julia with a working DISPLAY.")
        display_process = open(`$executable -displayfd 1 -screen 0 1280x1024x24 -nolisten tcp`, "r")
        atexit() do
            process_running(display_process) && kill(display_process)
            wait(display_process)
            close(display_process)
        end
        display_number = readline(display_process)
        tryparse(Int, display_number) === nothing && error("Xvfb did not provide a display number.")
        ENV["DISPLAY"] = ":$display_number"
    end
end
import GLMakie

const R_EARTH = 6_378_137.0
const R_ATMDEF = R_EARTH + 100_000.0
idx(satellite, component) = 6*(satellite-1) + component
include("functions/6_OE_and_dv_in_RTN.jl")
include(joinpath(@__DIR__, "..", "3_Kuang's Prototype Code_for_animation", "functions", "10_Animation_ver2.jl"))

const DEFAULT_BUNDLE = OPTIONS.bundle_directory

#########################
# Trajectory smoothing #
#########################
# Replay saved states without solving the orbital dynamics again.
struct ReplayTrajectory
    t::Vector{Float64}
    u::Vector{Vector{Float64}}
end

function (trajectory::ReplayTrajectory)(time)
    first(trajectory.t) <= time <= last(trajectory.t) || throw(ArgumentError("Time is outside the saved trajectory"))
    left = searchsortedlast(trajectory.t, time)
    time == trajectory.t[left] && return copy(trajectory.u[left])
    step = trajectory.t[left+1] - trajectory.t[left]
    fraction = (time - trajectory.t[left]) / step
    state = similar(first(trajectory.u))
    # Cubic Hermite interpolation uses saved velocities to smooth positions;
    # differentiating the same curve gives the interpolated velocities.
    for satellite in 1:length(state)÷6, component in 1:3
        position_index = idx(satellite, component)
        velocity_index = position_index + 3
        position_left, position_right = trajectory.u[left][position_index], trajectory.u[left+1][position_index]
        velocity_left, velocity_right = trajectory.u[left][velocity_index], trajectory.u[left+1][velocity_index]
        state[position_index] = (2*fraction^3 - 3*fraction^2 + 1)*position_left +
            (fraction^3 - 2*fraction^2 + fraction)*step*velocity_left +
            (-2*fraction^3 + 3*fraction^2)*position_right + (fraction^3 - fraction^2)*step*velocity_right
        state[velocity_index] = (6*fraction^2 - 6*fraction)/step*position_left +
            (3*fraction^2 - 4*fraction + 1)*velocity_left +
            (-6*fraction^2 + 6*fraction)/step*position_right + (3*fraction^2 - 2*fraction)*velocity_right
    end
    return state
end

######################
# Load saved results #
######################
function load_replay_bundle(directory)
    ### 1. Validate metadata and identify the target ###
    trajectory_table = Arrow.Table(joinpath(directory, "trajectory.feather"))
    metadata = Arrow.getmetadata(trajectory_table)
    metadata === nothing && throw(ArgumentError("Trajectory metadata is required"))
    source = get(metadata, "source", "")
    source in ("Kuang prototype", "SpaceAGORA") ||
        throw(ArgumentError("This replay loader expects a Kuang prototype or SpaceAGORA bundle"))
    spacecraft_ids = sort!([parse(Int, matched.captures[1]) for column in propertynames(trajectory_table)
        for matched in (match(r"^sc(\d+)_pos_1$", string(column)),) if matched !== nothing])
    count = length(spacecraft_ids)
    count >= 2 && spacecraft_ids == collect(1:count) || throw(ArgumentError("Expected contiguous spacecraft IDs"))
    target_id = parse(Int, source == "SpaceAGORA" ? metadata["target_id"] : get(metadata, "target_id", "1"))
    target_id in spacecraft_ids || throw(ArgumentError("Invalid target ID"))

    ### 2. Reorder saved states: helpers first, target last for the renderer ###
    order = vcat(filter(!=(target_id), spacecraft_ids), target_id)
    inverse_order = invperm(order)
    times = Float64.(trajectory_table.time)
    length(times) >= 2 && all(isfinite, times) && all(diff(times) .> 0) ||
        throw(ArgumentError("Expected at least two strictly increasing finite times"))
    states = [zeros(6*count) for time in times]
    for (render_id, stored_id) in enumerate(order), (field, offset) in (("pos", 0), ("vel", 3)), component in 1:3
        values = getproperty(trajectory_table, Symbol("sc$(stored_id)_$(field)_$(component)"))
        for sample in eachindex(times)
            states[sample][idx(render_id, offset+component)] = values[sample]
        end
    end
    all(state -> all(isfinite, state), states) || throw(ArgumentError("Nonfinite trajectory state"))

    ### 3. Load recorded encounter and laser intervals in the same spacecraft order ###
    function read_intervals(filename)
        table = Arrow.Table(joinpath(directory, filename))
        interval_metadata = Arrow.getmetadata(table)
        interval_metadata !== nothing && get(interval_metadata, "scenario", nothing) == metadata["scenario"] ||
            throw(ArgumentError("Mismatched scenario in $filename"))
        intervals = [(pair=minmax(inverse_order[row.sc_a], inverse_order[row.sc_b]),
            encounter_id=Int(row.encounter_id), start=Float64(row.start_time_s),
            stop=Float64(row.end_time_s), end_clipped=Bool(row.end_clipped)) for row in eachrow(DataFrame(table))]
        all(interval -> interval.pair[1] != interval.pair[2] &&
            first(times) <= interval.start <= interval.stop <= last(times), intervals) ||
            throw(ArgumentError("Invalid interval in $filename"))
        return intervals
    end
    geometry = read_intervals("geometry_encounters.feather")
    lasers = read_intervals("laser_on.feather")
    all(laser -> any(encounter -> laser.encounter_id == encounter.encounter_id && laser.pair == encounter.pair &&
        encounter.start <= laser.start <= laser.stop <= encounter.stop, geometry), lasers) ||
        throw(ArgumentError("Laser interval has no containing encounter"))
    ### 4. Build display parameters and the scenario caption ###
    pairs = unique(interval.pair for interval in vcat(geometry, lasers))
    params = Dict{Symbol,Any}(:N=>count, :Pmatrix=>zeros(count, count),
        :cavity=>Dict(pair=>nothing for pair in pairs))
    scenario = metadata["scenario"]
    caption_parts = match(r"^(?:prototype|spaceagora)_(N\d+_h[^_]+_t[^_]+_ih[^_]+_it[^_]+_e[^_]+_nu[^_]+)_T.*?_m[^_]+kg_(.+)_J2", scenario)
    scenario_caption = caption_parts === nothing ? scenario : join(caption_parts.captures, "_")
    return (; sol=ReplayTrajectory(times, states), params, geometry, lasers,
        helper_num=count-1, stored_order=order, target_id, scenario_caption,
        earth_radius=parse(Float64, get(metadata, "earth_radius_m", string(source == "SpaceAGORA" ? 6_378_136.6 : R_EARTH))))
end

#########################
# Recorded link display #
#########################
function recorded_link_states(bundle, time)
    count = bundle.params[:N]
    active = zeros(Float32, count, count)
    encounters = falses(count, count)
    kinds = zeros(Int, count, count)
    contains(interval) = interval.start <= time < interval.stop ||
        (interval.end_clipped && time == interval.stop)
    for interval in bundle.geometry
        if contains(interval)
            first_id, second_id = interval.pair
            encounters[first_id, second_id] = encounters[second_id, first_id] = true
        end
    end
    # A firing laser also remains a detected encounter; the indicators are not exclusive.
    for interval in bundle.lasers
        if contains(interval)
            first_id, second_id = interval.pair
            active[first_id, second_id] = active[second_id, first_id] = 1f0
            kinds[first_id, second_id] = kinds[second_id, first_id] = 2
            encounters[first_id, second_id] = encounters[second_id, first_id] = true
        end
    end
    return (; active, encounters, kinds)
end

###################
# Render animation #
###################
# By default, save beside the Feather output tree under videos/<scenario>/.
function replay_feather_animation(directory=DEFAULT_BUNDLE;
        duration_seconds=OPTIONS.duration_seconds, animation_fps=OPTIONS.animation_fps,
    output_file=joinpath(dirname(dirname(abspath(directory))), "videos", basename(normpath(directory)), basename(normpath(directory)) * ".mp4"),
        render_options...)
    ### 1. Load the saved bundle and configure the view ###
    bundle = load_replay_bundle(directory)
    GLMakie.activate!()
    defaults = (show_earth=true, helper_trails=false, target_trail=true,
        scenario_caption=bundle.scenario_caption,
        show_projections=true, helper_projections=false, target_projections=true,
        radial_exaggeration=40.0, earth_radius=bundle.earth_radius,
        reference_altitude_km=(norm(first(bundle.sol.u)[1:3])-bundle.earth_radius)/1e3,
        radial_ticks_km=[500.0, 1000.0, 1010.0, 1020.0, 1030.0, 1040.0, 1050.0])
    ### 2. Render frames using recorded events rather than recomputed laser forces ###
    fig, controls = animate_all_satellites_3d_smooth_helper_target(bundle.sol, bundle.params, bundle.helper_num;
        output_file, duration_seconds, animation_fps,
        recorded_links=time -> recorded_link_states(bundle, time), merge(defaults, (; render_options...))...)
    return (; fig, controls, output_file, bundle)
end

######################
# Command-line entry #
######################
# With no arguments, use OPTIONS above. Positional arguments override each setting:
# replay_feather_animation.jl [BUNDLE_DIRECTORY [PLAYBACK_SECONDS [FPS]]]
function main(arguments=ARGS)
    length(arguments) <= 3 || error("Usage: replay_feather_animation.jl [BUNDLE_DIRECTORY [PLAYBACK_SECONDS [FPS]]]")
    directory = isempty(arguments) ? DEFAULT_BUNDLE : abspath(arguments[1])
    duration_seconds = length(arguments) >= 2 ? parse(Float64, arguments[2]) : OPTIONS.duration_seconds
    animation_fps = length(arguments) >= 3 ? parse(Int, arguments[3]) : OPTIONS.animation_fps
    return replay_feather_animation(directory; duration_seconds, animation_fps)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    FeatherAnimationReplay.main()
elseif isinteractive() && isempty(PROGRAM_FILE)
    FeatherAnimationReplay.main(String[])
end