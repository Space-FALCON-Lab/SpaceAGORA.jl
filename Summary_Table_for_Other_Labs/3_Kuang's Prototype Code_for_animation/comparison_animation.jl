import GLMakie
using GeometryBasics, FileIO
GLMakie.activate!()
include("functions/10_Animation_ver2.jl")

function run_animated_comparison(settings, smoke)
    video_options = (
        show_earth = true,
        helper_trails = false,
        target_trail = true,
        show_projections = true,
        helper_projections = false,
        target_projections = true,
        radial_exaggeration = 40.0,
        reference_altitude_km = settings.helper_altitude_km,
        axis_limit_km = nothing,
        radial_ticks_km = [500.0, 1000.0, 1010.0, 1020.0, 1030.0, 1040.0, 1050.0],
        animation_fps = smoke ? 2 : 30,
    )
    return run_prototype(settings, smoke; video=true,
        video_duration_seconds=smoke ? 2.0 : 500.0, video_options)
end