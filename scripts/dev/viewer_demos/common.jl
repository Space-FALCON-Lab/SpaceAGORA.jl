# Shared helpers for the mission-recreation viewer demos (Magellan at Venus,
# Cassini at Titan, Apollo 11 at the Moon): mission SPK download, SPICE
# state sampling relative to the central body in J2000 (the frame the
# integrator and the viewer use), apoapsis / closest-approach searches for
# the initial epoch, and the "ghost" reference table for
# `export_visualization(...; references=...)`.
#
# Everything here is dev tooling: it needs the GRAM/SPICE asset tier and
# network access the first time a mission kernel is fetched.
const SPICE_PATH = abspath(get(ENV, "SPACEAGORA_SPICE_PATH", joinpath(@__DIR__, "..", "..", "..", "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE")))
include(joinpath(@__DIR__, "..", "..", "..", "examples", "common.jl"))
include(joinpath(@__DIR__, "demo_options.jl"))

using SPICE
using StaticArrays
using LinearAlgebra
using Downloads
using JSON
using Arrow
using DataFrames

const DEMO_OUT_ROOT = get(ENV, "SPACEAGORA_VIEWER_DEMO_OUT", joinpath(REPO_ROOT, "output", "viewer_demos"))
const MISSION_SPK_DIR = abspath(get(ENV, "SPACEAGORA_MISSION_SPK_DIR", joinpath(SPICE_PATH, "spk", "missions")))
const MODELS_DIR = joinpath(REPO_ROOT, "data", "models")
const HARMONICS_DIR = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data")

demo_outdir(name::AbstractString) = viewer_demo_options(name, 1.0; argv=String[]).output_dir

"""
    ensure_mission_kernel(filename, url) -> path

The mission SPK under `SPICE/spk/missions/`, downloaded from NAIF when absent.
"""
function ensure_mission_kernel(filename::AbstractString, url::AbstractString)::String
    basename(filename) == filename && !isempty(filename) && filename ∉ (".", "..") ||
        throw(ArgumentError("mission kernel filename must be a plain basename"))
    startswith(url, "https://naif.jpl.nasa.gov/") || throw(ArgumentError("mission kernels must come from the NAIF HTTPS archive"))
    path = joinpath(MISSION_SPK_DIR, String(filename))
    if !isfile(path)
        mkpath(MISSION_SPK_DIR)
        tmp, io = mktemp(MISSION_SPK_DIR); close(io)
        try
            println("Downloading ", url, " -> ", path)
            Downloads.download(String(url), tmp)
            filesize(tmp) > 0 || throw(ArgumentError("downloaded mission kernel is empty"))
            mv(tmp, path; force=false)
        finally
            isfile(tmp) && rm(tmp)
        end
    end
    filesize(path) > 0 || throw(ArgumentError("cached mission kernel is empty: $path"))
    return path
end

"""
A mission trajectory in SPICE: the spacecraft `target`, the body `center`
whose center the states are taken relative to (the simulation's primary),
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

"SPICE ET at the simulator's millisecond-resolved `SM.InitialTime`."
et_of(t::SM.InitialTime)::Float64 = SM.ephemerides_time_seconds(t, SM.SpiceEphemeridesModel())

et_of(utc::AbstractString)::Float64 = lock(RuntimeServices.SPICE_LOCK) do
    str2et(String(utc))
end

"""
Calendar fields for an ET, stored with `InitialTime`'s Float32 seconds.
Resolve these fields with `et_of(initial_time)` before sampling an initial state:
the engine and saved scene use millisecond resolution, not the search's raw ET.
"""
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

"Spacecraft position (m) and velocity (m/s) relative to the center in J2000 at `et`."
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

# Searches are local bracket refinements, not guarantees of a global extremum.
function _search_bounds(a, b, tol)
    all(isfinite, (a, b, tol)) && a < b && tol > 0 && isfinite(b - a) ||
        throw(ArgumentError("search bounds must be finite, increasing, with positive tolerance"))
    return nothing
end
_finite_objective(f, t) = (value = Float64(f(t)); isfinite(value) || throw(ArgumentError("nonfinite search objective at $t")); value)
function _search_grid(a::Float64, b::Float64, step_s::Float64)
    _search_bounds(a, b, step_s)
    count = ceil((b - a) / step_s)
    isfinite(count) && count <= 100_000 || throw(ArgumentError("search needs too many samples; widen step_s or narrow window_s"))
    grid = collect(range(a, b; length=max(2, Int(count) + 1)))
    all(diff(grid) .> 0) || throw(ArgumentError("search spacing is below floating-point time resolution"))
    return grid
end
function _golden_min(f, a::Float64, b::Float64; tol::Float64=1e-3)
    _search_bounds(a, b, tol)
    φ = (sqrt(5.0) - 1.0) / 2.0
    c = b - φ * (b - a); d = a + φ * (b - a)
    fc = _finite_objective(f, c); fd = _finite_objective(f, d)
    for _ in 1:256
        b - a <= tol && return a + (b - a) / 2
        a < c < d < b || throw(ArgumentError("minimum-search tolerance is below floating-point time resolution"))
        if fc < fd
            b, d, fd = d, c, fc
            c = b - φ * (b - a); fc = _finite_objective(f, c)
        else
            a, c, fc = c, d, fd
            d = a + φ * (b - a); fd = _finite_objective(f, d)
        end
    end
    throw(ArgumentError("minimum search exceeded its iteration limit"))
end
function _bisect(g, a::Float64, b::Float64; tol::Float64=1e-4)
    _search_bounds(a, b, tol)
    ga = _finite_objective(g, a); gb = _finite_objective(g, b)
    ga >= 0 && gb < 0 || throw(ArgumentError("apoapsis bracket requires nonnegative then negative radial velocity"))
    ga == 0 && return a
    for _ in 1:256
        b - a <= tol && return a + (b - a) / 2
        c = a + (b - a) / 2
        a < c < b || throw(ArgumentError("root-search tolerance is below floating-point time resolution"))
        gc = _finite_objective(g, c)
        gc == 0 && return c
        if gc > 0
            a = c
        else
            b = c
        end
    end
    throw(ArgumentError("apoapsis search exceeded its iteration limit"))
end
function _closest_approach_et(distance, et_guess::Float64; window_s::Float64=6 * 3600.0, step_s::Float64=60.0)
    isfinite(et_guess) && isfinite(window_s) && window_s > 0 || throw(ArgumentError("finite epoch and positive window_s required"))
    grid = _search_grid(et_guess - window_s, et_guess + window_s, step_s)
    d = [_finite_objective(distance, et) for et in grid]
    k = argmin(d)
    1 < k < length(grid) || throw(ArgumentError("closest approach is not bracketed inside the search window; widen or recenter it"))
    return _golden_min(distance, grid[k - 1], grid[k + 1])
end
"Find a bracketed local closest approach to the central body, near et_guess."
closest_approach_et(m::MissionSpice, et_guess::Float64; kwargs...) =
    _closest_approach_et(et -> distance_m(m, et), et_guess; kwargs...)
function _first_apoapsis_et_after(radial, et_from::Float64; window_s::Float64=86_400.0, step_s::Float64=120.0)
    isfinite(et_from) && isfinite(window_s) && window_s > 0 || throw(ArgumentError("finite epoch and positive window_s required"))
    grid = _search_grid(et_from, et_from + window_s, step_s)
    rv = [_finite_objective(radial, et) for et in grid]
    for i in 1:(length(grid) - 1)
        if rv[i] >= 0 && rv[i + 1] < 0
            return _bisect(radial, grid[i], grid[i + 1])
        end
    end
    throw(ArgumentError("no sampled positive-to-negative apoapsis crossing within the search window"))
end
"First sampled positive-to-negative radial-velocity crossing at or after et_from."
first_apoapsis_et_after(m::MissionSpice, et_from::Float64; kwargs...) =
    _first_apoapsis_et_after(et -> radial_velocity(m, et), et_from; kwargs...)

"Cartesian initial condition from the mission SPK at `et`."
function cartesian_ic_at(m::MissionSpice, et::Float64; kwargs...)::SM.CartesianInitialCondition
    r, v = spice_state_m(m, et)
    return SM.CartesianInitialCondition(r, v; kwargs...)
end

"""
    spice_reference(m, prefix; name, target=1, color, stride=1) -> NamedTuple

The ghost table for `export_visualization(prefix; references=[...])`: the
mission SPK sampled at the run's saved times (from `<prefix>.feather`),
relative to the center in J2000, positions and velocities in SI units.
Times outside the kernel's coverage are dropped.
"""
function spice_reference(m::MissionSpice, prefix::AbstractString; name::AbstractString=m.name, target::Int=1,
                         color::AbstractString="#ff8c69", opacity::Real=0.45, stride::Int=1)
    stride >= 1 || throw(ArgumentError("reference stride must be positive"))
    scene = JSON.parsefile(String(prefix) * "_scene.json")
    1 <= target <= length(scene["spacecraft"]) || throw(ArgumentError("reference target is a one-based spacecraft index"))
    isfinite(opacity) && 0 <= opacity <= 1 || throw(ArgumentError("reference opacity must lie in [0, 1]"))
    uppercase(strip(m.center)) == uppercase(String(scene["planet"]["name"])) || throw(ArgumentError("SPICE reference central body differs from the saved scene"))
    et0 = Float64(scene["epoch"]["et_start_s"])
    df = DataFrame(Arrow.Table(String(prefix) * ".feather"))
    all_times = Float64.(df.time)
    !isempty(all_times) && all(isfinite, all_times) && issorted(all_times) || throw(ArgumentError("saved reference times must be finite, nonempty and ordered"))
    isfinite(et0) || throw(ArgumentError("scene epoch must be finite"))
    times = all_times[1:stride:end]
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

Run with the visualization sidecar into a fresh directory. The legacy name is
retained for callers, but results are never silently reused or overwritten.
"""
function run_or_reuse!(args, outdir::AbstractString; kwargs...)::String
    output = abspath(outdir)
    output == abspath(args.simulation_settings.results_directory) || throw(ArgumentError("outdir must match the configuration results_directory"))
    args.simulation_settings.results || throw(ArgumentError("mission viewer demos require results=true"))
    args.simulation_settings.generate_filenames && throw(ArgumentError("mission demos require the standard simulation_results prefix"))
    ispath(output) && (!isdir(output) || !isempty(readdir(output))) && throw(ArgumentError("output must be absent or empty; choose a new --output-dir: $output"))
    mkpath(output)
    elapsed = @elapsed Base.invokelatest(run_simulation, args; visualization=true, kwargs...)
    println("simulation: ", round(elapsed; digits=1), " s")
    prefix = joinpath(output, "simulation_results")
    all(isfile, (prefix * ".feather", prefix * "_scene.json")) || error("simulation did not produce the required viewer results")
    return prefix
end

"Summary of sampled radius-minus-equatorial-radius altitude and downward threshold crossings, not geodetic altitude or exact periapsis."
function summarize_run(prefix::AbstractString, planet; alt_m::Float64=250e3)
    df = DataFrame(Arrow.Table(String(prefix) * ".feather"))
    r = sqrt.(df.sc1_pos_1 .^ 2 .+ df.sc1_pos_2 .^ 2 .+ df.sc1_pos_3 .^ 2)
    alt = r .- planet.Rp_e
    isempty(r) && throw(ArgumentError("cannot summarize an empty saved trajectory"))
    passes = count(i -> alt[i] < alt_m && alt[i - 1] >= alt_m, 2:nrow(df))
    println("rows=", nrow(df), " span=", round(df.time[end] / 3600; digits=2), " h  sampled spherical alt min/max km=",
        round(minimum(alt) / 1e3; digits=1), "/", round(maximum(alt) / 1e3; digits=1), "  passes below ", alt_m / 1e3, " km: ", passes)
    return (rows=nrow(df), alt_min_m=minimum(alt), alt_max_m=maximum(alt), passes=passes)
end

"Largest separation (m) between the run and a reference table over the shared times."
function reference_separation(prefix::AbstractString, ref)
    df = DataFrame(Arrow.Table(String(prefix) * ".feather"))
    target = Int(ref.target)
    target >= 1 || throw(ArgumentError("reference target must be a positive spacecraft index"))
    columns = [Symbol("sc$(target)_pos_$(axis)") for axis in 1:3]
    all(c -> c in propertynames(df), columns) || throw(ArgumentError("saved results lack position columns for reference target $target"))
    length(ref.t_s) == size(ref.pos_m, 2) && size(ref.pos_m, 1) == 3 || throw(ArgumentError("reference positions must be 3 by number of times"))
    all(isfinite, ref.t_s) && all(isfinite, ref.pos_m) || throw(ArgumentError("reference samples must be finite"))
    length(unique(df.time)) == nrow(df) || throw(ArgumentError("separation report requires unique saved times"))
    sim = Dict(t => SVector((df[i, c] for c in columns)...) for (i, t) in enumerate(df.time))
    seps = [norm(sim[t] - SVector(ref.pos_m[1, k], ref.pos_m[2, k], ref.pos_m[3, k])) for (k, t) in enumerate(ref.t_s) if haskey(sim, t)]
    isempty(seps) && throw(ArgumentError("reference and results share no exact saved times"))
    println("sampled separation vs ", ref.name, ": start ", round(seps[1]; digits=1), " m, max ", round(maximum(seps) / 1e3; digits=2), " km, end ", round(seps[end] / 1e3; digits=2), " km")
    return seps
end

"Explicit SPACEAGORA_DEMO_CDN=1 opt-in to the separately supplied CDN tool; default is standalone HTML only."
function build_cdn_page(html::AbstractString, out::AbstractString, title::AbstractString, heading::AbstractString,
                        orbit::AbstractString, span::AbstractString, foot::AbstractString)
    get(ENV, "SPACEAGORA_DEMO_CDN", "0") == "1" || return nothing
    script = joinpath(REPO_ROOT, "viewer", "build_cdn_page.py")
    isfile(script) || throw(ArgumentError("optional CDN export requires the separately installed build_cdn_page.py tool; standalone HTML is already complete"))
    run(`python3 $script $html $out $title $heading $orbit $span $foot`)
    return String(out)
end

"""
    demo_aero_effector(model_path, outdir; scale, rotation_deg, reference_area_m2, wall_temperature_k) -> effector

`AerodynamicCoefficientfM()` by default. Explicit mesh mode fits the current
mission geometry afresh with the separately supplied mesh-aerodynamics feature.
It never silently reuses a fitted model from another configuration.
"""
function demo_aero_effector(model_path::AbstractString, outdir::AbstractString; scale::Real=1.0, rotation_deg=(0.0, 0.0, 0.0),
                            reference_area_m2=nothing, wall_temperature_k::Real=300.0, degree::Int=10, articulations=())
    get(ENV, "SPACEAGORA_DEMO_MESH_AERO", "0") == "1" || return AerodynamicCoefficientfM()
    all(name -> isdefined(@__MODULE__, name), (:mesh_aero_panels, :fit_mesh_aero_surrogate, :AerodynamicCoefficientMeshSurrogate)) ||
        throw(ArgumentError("mesh mode requires the separately supplied mesh-aerodynamics feature"))
    panels = mesh_aero_panels(model_path; scale=scale, rotation_deg=rotation_deg, reference_area_m2=reference_area_m2, articulations=articulations)
    println("mesh aero: ", length(panels), " facets from ", basename(model_path))
    elapsed = @elapsed surrogate = fit_mesh_aero_surrogate(panels; degree=degree, n_directions=1200, verbose=true)
    println("mesh aero: fitted in ", round(elapsed; digits=1), " s")
    return AerodynamicCoefficientMeshSurrogate(surrogate; wall_temperature_k=wall_temperature_k)
end

"Results directory suffix so the mesh-aero variant of a case lands beside the box-model run."
function demo_case_name(name::AbstractString)
    mode = get(ENV, "SPACEAGORA_DEMO_MESH_AERO", "0")
    mode in ("0", "1") || throw(ArgumentError("SPACEAGORA_DEMO_MESH_AERO must be 0 or 1"))
    return mode == "1" ? String(name) * "_mesh_aero" : String(name)
end

function demo_model_path(filename::AbstractString)
    path = joinpath(MODELS_DIR, filename)
    isfile(path) && filesize(path) > 0 || throw(ArgumentError("mission model asset missing: $path; install the mission-assets package"))
    return path
end
