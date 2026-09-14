# Shared helpers for the mission-recreation viewer demos (Magellan at Venus,
# Cassini at Titan, Apollo 11 at the Moon): mission SPK download, SPICE
# state sampling relative to the central body in J2000 (the frame the
# integrator and the viewer use), apoapsis / closest-approach searches for
# the initial epoch, and the "ghost" reference table for
# `export_visualization(...; references=...)`.
#
# Everything here is dev tooling: it needs the GRAM/SPICE asset tier and
# network access the first time a mission kernel is fetched.
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(REPO_ROOT, "examples", "common.jl"))

using SPICE
using StaticArrays
using LinearAlgebra
using Downloads
using JSON
using Arrow
using DataFrames

const DEMO_OUT_ROOT = get(ENV, "SPACEAGORA_VIEWER_DEMO_OUT", joinpath(REPO_ROOT, "output", "viewer_demos"))
const MISSION_SPK_DIR = joinpath(SPICE_PATH, "spk", "missions")
const MODELS_DIR = joinpath(REPO_ROOT, "data", "models")
const HARMONICS_DIR = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data")

demo_outdir(name::AbstractString) = (d = joinpath(DEMO_OUT_ROOT, String(name)); mkpath(d); d)

"""
    ensure_mission_kernel(filename, url) -> path

The mission SPK under `SPICE/spk/missions/`, downloaded from NAIF when absent.
"""
function ensure_mission_kernel(filename::AbstractString, url::AbstractString)::String
    path = joinpath(MISSION_SPK_DIR, String(filename))
    if !isfile(path)
        mkpath(MISSION_SPK_DIR)
        println("Downloading ", url, " -> ", path)
        tmp = path * ".part"
        Downloads.download(String(url), tmp)
        mv(tmp, path; force=true)
    end
    return path
end

"""
A mission trajectory in SPICE: the spacecraft `target`, the body `center`
whose centre the states are taken relative to (the simulation's primary),
and the kernels that carry them.
"""
struct MissionSpice
    name::String
    target::String
    center::String
    kernels::Vector{String}
end

function furnish!(m::MissionSpice)
    lock(RuntimeServices.SPICE_LOCK) do
        foreach(furnsh, m.kernels)
    end
    return m
end

"Elapsed-time ET of a `SM.InitialTime`."
function et_of(t::SM.InitialTime)::Float64
    stamp = string(lpad(t.year, 4, '0'), "-", lpad(t.month, 2, '0'), "-", lpad(t.day, 2, '0'), "T",
        lpad(t.hour, 2, '0'), ":", lpad(t.minute, 2, '0'), ":", lpad(floor(Int, t.second), 2, '0'))
    frac = Float64(t.second) - floor(Float64(t.second))
    return lock(RuntimeServices.SPICE_LOCK) do
        str2et(stamp) + frac
    end
end

et_of(utc::AbstractString)::Float64 = lock(RuntimeServices.SPICE_LOCK) do
    str2et(String(utc))
end

"`SM.InitialTime` of an ET (microsecond resolution)."
function initial_time_of(et::Float64)::SM.InitialTime
    utc = lock(RuntimeServices.SPICE_LOCK) do
        et2utc(et, "ISOC", 6)
    end
    d, tm = split(utc, 'T')
    y, mo, dd = parse.(Int, split(d, '-'))
    hh, mm, ss = split(tm, ':')
    return SM.InitialTime(year=y, month=mo, day=dd, hour=parse(Int, hh), minute=parse(Int, mm), second=parse(Float64, ss))
end

utc_of(et::Float64)::String = lock(RuntimeServices.SPICE_LOCK) do
    et2utc(et, "ISOC", 3)
end

"Spacecraft position (m) and velocity (m/s) relative to the centre in J2000 at `et`."
function spice_state_m(m::MissionSpice, et::Float64)
    st = lock(RuntimeServices.SPICE_LOCK) do
        spkezr(m.target, et, "J2000", "NONE", m.center)[1]
    end
    return SVector{3, Float64}(st[1], st[2], st[3]) .* 1e3, SVector{3, Float64}(st[4], st[5], st[6]) .* 1e3
end

function radial_velocity(m::MissionSpice, et::Float64)::Float64
    r, v = spice_state_m(m, et)
    return dot(r, v) / norm(r)
end

distance_m(m::MissionSpice, et::Float64)::Float64 = norm(spice_state_m(m, et)[1])

# Golden-section minimum of f on [a, b].
function _golden_min(f, a::Float64, b::Float64; tol::Float64=1e-3)
    φ = (sqrt(5.0) - 1.0) / 2.0
    c = b - φ * (b - a); d = a + φ * (b - a)
    fc = f(c); fd = f(d)
    while abs(b - a) > tol
        if fc < fd
            b, d, fd = d, c, fc
            c = b - φ * (b - a); fc = f(c)
        else
            a, c, fc = c, d, fd
            d = a + φ * (b - a); fd = f(d)
        end
    end
    return (a + b) / 2
end

# Bisection root of g on [a, b] with g(a) > 0 > g(b).
function _bisect(g, a::Float64, b::Float64; tol::Float64=1e-4)
    ga = g(a)
    while b - a > tol
        c = (a + b) / 2
        gc = g(c)
        if (gc > 0) == (ga > 0)
            a, ga = c, gc
        else
            b = c
        end
    end
    return (a + b) / 2
end

"""
    closest_approach_et(m, et_guess; window_s=6h, step_s=60) -> et

ET of the minimum distance to the centre within `et_guess ± window_s`.
"""
function closest_approach_et(m::MissionSpice, et_guess::Float64; window_s::Float64=6 * 3600.0, step_s::Float64=60.0)::Float64
    grid = collect((et_guess - window_s):step_s:(et_guess + window_s))
    d = [distance_m(m, et) for et in grid]
    k = argmin(d)
    lo = grid[max(1, k - 1)]; hi = grid[min(length(grid), k + 1)]
    return _golden_min(et -> distance_m(m, et), lo, hi)
end

"""
    first_apoapsis_et_after(m, et_from; window_s=1 day, step_s=120) -> et

First apoapsis (radial velocity crossing from positive to negative) after `et_from`.
"""
function first_apoapsis_et_after(m::MissionSpice, et_from::Float64; window_s::Float64=86_400.0, step_s::Float64=120.0)::Float64
    grid = collect(et_from:step_s:(et_from + window_s))
    rv = [radial_velocity(m, et) for et in grid]
    for i in 1:(length(grid) - 1)
        if rv[i] >= 0.0 && rv[i + 1] < 0.0
            return _bisect(et -> radial_velocity(m, et), grid[i], grid[i + 1])
        end
    end
    throw(ArgumentError("$(m.name): no apoapsis within $(window_s) s after ET $(et_from)."))
end

"Cartesian initial condition from the mission SPK at `et`."
function cartesian_ic_at(m::MissionSpice, et::Float64; kwargs...)::SM.CartesianInitialCondition
    r, v = spice_state_m(m, et)
    return SM.CartesianInitialCondition(r, v; kwargs...)
end

"""
    spice_reference(m, prefix; name, target=1, color, stride=1) -> NamedTuple

The ghost table for `export_visualization(prefix; references=[...])`: the
mission SPK sampled at the run's saved times (from `<prefix>.feather`),
relative to the centre in J2000, positions and velocities in SI units.
Times outside the kernel's coverage are dropped.
"""
function spice_reference(m::MissionSpice, prefix::AbstractString; name::AbstractString=m.name, target::Int=1,
                         color::AbstractString="#ff8c69", opacity::Real=0.45, stride::Int=1)
    scene = JSON.parsefile(String(prefix) * "_scene.json")
    et0 = Float64(scene["epoch"]["et_start_s"])
    df = DataFrame(Arrow.Table(String(prefix) * ".feather"))
    times = Float64.(df.time)[1:stride:end]
    keep = Int[]
    pos = zeros(3, length(times)); vel = zeros(3, length(times))
    for (k, t) in enumerate(times)
        try
            r, v = spice_state_m(m, et0 + t)
            pos[:, k] .= r; vel[:, k] .= v
            push!(keep, k)
        catch err
            occursin("SPICE(SPKINSUFFDATA)", sprint(showerror, err)) || rethrow()
        end
    end
    isempty(keep) && throw(ArgumentError("$(m.name): the kernel covers none of the run's saved times."))
    length(keep) == length(times) || println("  reference '", name, "': kernel covers ", length(keep), " of ", length(times), " saved times")
    return (name=String(name), t_s=times[keep], pos_m=pos[:, keep], vel_mps=vel[:, keep], target=target, color=String(color), opacity=Float64(opacity))
end

"""
    run_or_reuse!(args, outdir) -> prefix

Run the simulation with the visualization sidecar unless `outdir` already
holds one (set `SPACEAGORA_DEMO_FORCE=1` to rerun).
"""
function run_or_reuse!(args, outdir::AbstractString)::String
    prefix = joinpath(outdir, "simulation_results")
    if isfile(prefix * "_scene.json") && get(ENV, "SPACEAGORA_DEMO_FORCE", "0") != "1"
        println("results present in ", outdir, "; skipping the simulation")
    else
        t = @elapsed run_simulation(args; visualization=true)
        println("simulation: ", round(t; digits=1), " s")
    end
    return prefix
end

"Summary of the saved trajectory: rows, span, altitude range, periapsis passes below `alt_m`."
function summarize_run(prefix::AbstractString, planet; alt_m::Float64=250e3)
    df = DataFrame(Arrow.Table(String(prefix) * ".feather"))
    r = sqrt.(df.sc1_pos_1 .^ 2 .+ df.sc1_pos_2 .^ 2 .+ df.sc1_pos_3 .^ 2)
    alt = r .- planet.Rp_e
    passes = count(i -> alt[i] < alt_m && alt[i - 1] >= alt_m, 2:nrow(df))
    println("rows=", nrow(df), " span=", round(df.time[end] / 3600; digits=2), " h  alt min/max km=",
        round(minimum(alt) / 1e3; digits=1), "/", round(maximum(alt) / 1e3; digits=1), "  passes below ", alt_m / 1e3, " km: ", passes)
    return (rows=nrow(df), alt_min_m=minimum(alt), alt_max_m=maximum(alt), passes=passes)
end

"Largest separation (m) between the run and a reference table over the shared times."
function reference_separation(prefix::AbstractString, ref)
    df = DataFrame(Arrow.Table(String(prefix) * ".feather"))
    sim = Dict(t => SVector(df.sc1_pos_1[i], df.sc1_pos_2[i], df.sc1_pos_3[i]) for (i, t) in enumerate(df.time))
    seps = [norm(sim[t] - SVector(ref.pos_m[1, k], ref.pos_m[2, k], ref.pos_m[3, k])) for (k, t) in enumerate(ref.t_s) if haskey(sim, t)]
    println("separation vs ", ref.name, ": start ", round(seps[1]; digits=1), " m, max ", round(maximum(seps) / 1e3; digits=2), " km, end ", round(seps[end] / 1e3; digits=2), " km")
    return seps
end

"Build the CDN (artifact) variant of an exported page with scripts/dev/viewer_demos/build_cdn_page.py."
function build_cdn_page(html::AbstractString, out::AbstractString, title::AbstractString, heading::AbstractString,
                        orbit::AbstractString, span::AbstractString, foot::AbstractString)
    script = joinpath(@__DIR__, "build_cdn_page.py")
    run(`python3 $script $html $out $title $heading $orbit $span $foot`)
    return String(out)
end
