# The CYGNSS constellation flown together over the window its flight telemetry
# covers: 2025-06-06T00:00:00Z to 2025-06-09T24:00:00Z, 96 hours.
#
# Seven spacecraft, not the eight that were launched: CYGNSS FM06 (NORAD 41889)
# was no longer in orbit at this epoch. NASA lost contact with it in November
# 2022 and the satellite catalogue gives its decay date as 2024-06-13. The
# initial-conditions file's notes carry the citations.
#
# All seven are flown states. Every initial condition and every reference ghost
# on this page comes from the NASA CYGNSS Level 1 navigation solution for that
# spacecraft over this window; nothing is a catalogue propagation and nothing is
# a design orbit. There is therefore nothing to distinguish between and no
# provenance to colour-code: the page says once that the whole constellation is
# flown.
#
# The constellation has no propulsion, so by June 2025 drag had taken it well
# below its design orbit: 438 to 448 km, inclination 34.88 to 34.97 degrees,
# periods of 93.4 to 93.5 minutes, against the 510 km and 95 minutes of the
# design. Nothing here sanity-checks a flown state against the design values.
#
# Initial states come from `data/telemetry/CYGNSS/constellation_ics_20250606.json`
# and the reference ghosts from
# `data/telemetry/CYGNSS/cygnss_constellation_tracks_20250606_96hr.feather`
# (both gitignored; both built by `build_cygnss_ics.jl`). Samples the track
# table marks `arc_consistent = false` are skipped: about two in a thousand, the
# product's own bad navigation fixes, which a ghost would otherwise draw as a
# spike.
#
# Force model, following `docs/spaceagora_cygnss_reconstruction_record.md`
# section 3: EarthGGM05C gravity to degree and order 50, Sun and Moon third
# body, solar radiation pressure, and NRLMSISE-00 drag with real CelesTrak
# space-weather indices.
#
# Spacecraft geometry. The drawn and simulated vehicle is a GENERIC small
# satellite box built only from publicly published CYGNSS figures (body
# roughly 51 x 64 x 28 cm and a 1.67 m deployed array span, NASA/eoPortal;
# about 29 kg, the mass in the reconstruction record). It is deliberately not
# a to-scale reconstruction, and no restricted mission-configuration geometry
# is read, embedded or published here.
#
#   julia --project=. scripts/dev/viewer_demos/cygnss_constellation.jl
include(joinpath(@__DIR__, "common.jl"))
# The initial-conditions file's own loader and validator, written and owned by
# `build_cygnss_ics.jl` alongside it: it raises rather than handing back a
# state with no provenance, which is exactly what this page must never draw.
include(joinpath(@__DIR__, "cygnss_ics.jl"))
include(joinpath(@__DIR__, "cygnss_tracks.jl"))
using .CygnssICs: load_constellation_ics
using .CygnssTracks: load_constellation_tracks, track_names
using Dates
using Printf
using Statistics

import SpaceAGORA.TelemetryVerification as TV

const OUTDIR = demo_outdir("cygnss_constellation")
const TELEMETRY_DIR = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS")
const IC_PATH = get(ENV, "SPACEAGORA_DEMO_CYGNSS_ICS", joinpath(TELEMETRY_DIR, "constellation_ics_20250606.json"))
const EPOCH_UTC = "2025-06-06T00:00:00"
# SPACEAGORA_DEMO_CYGNSS_HOURS shortens the window for a smoke run; the page of
# record is the full 96 h.
const WINDOW_S = parse(Float64, get(ENV, "SPACEAGORA_DEMO_CYGNSS_HOURS", "96.0")) * 3600.0
# The flown position solution of every spacecraft over the window, at 1 Hz in
# J2000, written by `build_cygnss_ics.jl` from the NASA Level 1 product.
const TRACKS_PATH = get(ENV, "SPACEAGORA_DEMO_CYGNSS_TRACKS",
    joinpath(TELEMETRY_DIR, "cygnss_constellation_tracks_20250606_96hr.feather"))

# --- spacecraft, from published figures only -------------------------------
# Body 51 x 64 x 28 cm and a 1.67 m deployed array span are NASA/eoPortal
# CYGNSS mission-page figures; 29.00 kg is the mass stated in the SpaceAGORA
# CYGNSS reconstruction record, section 2. The array CHORD is not published
# in either source: 0.28 m (the body height) is an assumption, made so the
# wings are of the same order as the body they fold against, and it enters
# only through the drag and SRP reference area.
const CYGNSS_BUS_DIMS_M = (0.51, 0.64, 0.28)
const CYGNSS_ARRAY_SPAN_M = 1.67
const CYGNSS_ARRAY_CHORD_M = 0.28          # ASSUMPTION, not a published figure
const CYGNSS_BUS_MASS_KG = 29.0
const CYGNSS_ARRAY_MASS_EACH_KG = 0.0      # folded into the bus mass; the record's 29.00 kg is the whole vehicle
# Per-wing half-span: make_three_body_spacecraft gives EACH panel link the full
# dims[2]*dims[3] as its reference area, so this must be half the total span
# minus the body it grows out of.
const CYGNSS_ARRAY_HALF_SPAN_M = (CYGNSS_ARRAY_SPAN_M - CYGNSS_BUS_DIMS_M[2]) / 2
const CYGNSS_SRP_CR = 1.3                  # ASSUMPTION: a generic small-satellite reflectivity; the record states only "fixed coefficients"

# Starting effective drag scale. The reconstruction record calibrated about 0.3 against
# its own (restricted, non-public) geometry; that number does not transfer to
# the generic public box above, so it is NOT reused here. The default is 1.0 —
# an uncalibrated free-molecular coefficient on the published area — and the
# end-of-window separations the script prints are what they are.
const CYGNSS_DRAG_SCALE = parse(Float64, get(ENV, "SPACEAGORA_DEMO_CYGNSS_DRAG_SCALE", "1.0"))

# Fit the initial orbital energy and the effective drag scale of each spacecraft
# to that spacecraft's own flown track (see `fit_states!`).
# SPACEAGORA_DEMO_CYGNSS_FIT_SMA=0 propagates the initial-conditions file's
# states untouched, at the uncalibrated drag scale, instead.
const FIT_SMA = get(ENV, "SPACEAGORA_DEMO_CYGNSS_FIT_SMA", "1") != "0"


"""
    ScaledAero(model, scales)

`AerodynamicCoefficientfM` with its force and torque multiplied by a
PER-SPACECRAFT constant: `scales[i]` is the effective drag scale of spacecraft
`i`. It absorbs everything the generic box and the free-molecular coefficient
get wrong about the real vehicle's ballistic coefficient and about the
atmosphere the empirical model predicts, which is why the reconstruction record
carries one too (section 3, "calibrated effective drag scale"). The engine's
own `ScaledAerodynamicCoefficientfM` is internal to `TelemetryVerification` and
takes a single global scale, so the demo carries its own per-spacecraft copy.
"""
struct ScaledAero <: SM.AbstractForceTorqueModel
    model::SM.AerodynamicCoefficientfM
    scales::Vector{Float64}
end
@inline SM.environment_requirements(::ScaledAero) =
    SM.EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
function SM.calcForceTorque(m::ScaledAero, x::AbstractVector, param::SM.ODEParams, i::Int64)
    f, τ = SM.calcForceTorque(m.model, x, param, i)
    c = m.scales[min(i, length(m.scales))]
    return c .* f, c .* τ
end

# --- telemetry -------------------------------------------------------------

const _TRACKS = Ref{Any}(nothing)

"The track table, loaded and validated once per process."
function tracks()
    _TRACKS[] === nothing || return _TRACKS[]
    isfile(TRACKS_PATH) || throw(ArgumentError(
        "no flown-track table at $(TRACKS_PATH); build it with scripts/dev/viewer_demos/build_cygnss_ics.jl."))
    t = load_constellation_tracks(TRACKS_PATH; epoch_utc=EPOCH_UTC * "Z")
    _TRACKS[] = t
    return t
end

const _TRACK_CACHE = Dict{String, Any}()

"""
    track_series(name) -> (t_s, pos_m, vel_mps, note)

One spacecraft's flown arc out of the track table: elapsed seconds from the run
epoch, J2000 position and velocity, with the samples the table marks
`arc_consistent = false` removed.

That flag is the table's own check that a sample agrees with its immediate
neighbors to better than 40 m, several times the 9 m a one-second second
difference should show at this altitude. It is false on about two samples in a
thousand — the product's bad navigation fixes — and a ghost that connected them
would draw a visible spike.

A neighbor test cannot catch a RUN of bad fixes, because they agree with each
other. One such run survives it: five samples in the FM01 arc that place the
spacecraft up to 486 km off the orbit. So a second filter follows — drop samples
whose radius is more than 50 km from the arc's own median. A near-circular
440 km orbit varies by about five kilometers in radius over a revolution, so
that band keeps every real sample; across all seven arcs it removes those five
and nothing else.
"""
function track_series(name::AbstractString)
    get!(_TRACK_CACHE, String(name)) do
        df = tracks().table
        rows = findall(==(String(name)), df.name)
        isempty(rows) && throw(ArgumentError("no flown track for \"$(name)\" in $(TRACKS_PATH)"))
        flagged = [i for i in rows if df.arc_consistent[i]]
        radii = [hypot(df.pos_ii_1[i], df.pos_ii_2[i], df.pos_ii_3[i]) for i in flagged]
        r_med = median(radii)
        good = [flagged[k] for k in eachindex(flagged) if abs(radii[k] - r_med) <= 50.0e3]
        t = Float64.(df.time[good])
        pos = permutedims(hcat(Float64.(df.pos_ii_1[good]), Float64.(df.pos_ii_2[good]), Float64.(df.pos_ii_3[good])))
        vel = permutedims(hcat(Float64.(df.vel_ii_1[good]), Float64.(df.vel_ii_2[good]), Float64.(df.vel_ii_3[good])))
        (t_s=t, pos_m=pos, vel_mps=vel, label=String(name),
         note="NASA Level 1 navigation solution, $(length(good)) of $(length(rows)) samples " *
              "($(length(rows) - length(flagged)) flagged inconsistent, $(length(flagged) - length(good)) off the orbit)")
    end
end

"Whether the track table carries a flown arc for this spacecraft."
has_track(name::AbstractString) = String(name) in track_names(tracks())

# --- initial conditions ----------------------------------------------------

"""
    track_stub_states() -> Vector

The development stub used when `constellation_ics_20250606.json` is absent but
the flown track table is not, in that file's schema: every spacecraft's own
state at the common epoch, straight out of its track, all `telemetry`.

The track table carries a fix at exactly t = 0 for every spacecraft, so nothing
is propagated or extrapolated to reach the epoch. This is a convenience for
running the demo without the initial-conditions file; the file remains the
interface, since it is the thing that carries the provenance and the citations.
"""
function track_stub_states()
    out = Any[]
    for name in track_names(tracks())
        series = track_series(name)
        j = findfirst(t -> t >= 0.0, series.t_s)
        j === nothing && throw(ArgumentError("the track for $(name) does not reach the epoch"))
        r = SVector{3, Float64}(series.pos_m[:, j]); v = SVector{3, Float64}(series.vel_mps[:, j])
        push!(out, (name=name, norad_id=0, r_ii_m=collect(r), v_ii_m_s=collect(v),
            provenance="telemetry", ic=SM.CartesianInitialCondition(r, v),
            source="$(basename(TRACKS_PATH)) at t = $(series.t_s[j]) s"))
    end
    return out
end

"""
    load_states(planet) -> (states, stubbed::Bool)

The constellation initial conditions, from the JSON file when it exists and
from the track table when it does not.
"""
function load_states(planet)
    if !isfile(IC_PATH)
        println("!! ", IC_PATH, " is absent: taking every initial state from the flown track table instead.")
        return track_stub_states(), true
    end
    ics = load_constellation_ics(IC_PATH)
    startswith(ics.epoch_utc, EPOCH_UTC) || throw(ArgumentError(
        "$(IC_PATH) is at epoch $(ics.epoch_utc); this demo propagates from $(EPOCH_UTC)Z."))
    states = Any[(name=sc.name, norad_id=sc.norad_id === nothing ? 0 : sc.norad_id,
        r_ii_m=collect(sc.r_ii_m), v_ii_m_s=collect(sc.v_ii_m_s), provenance=sc.provenance,
        ic=SM.CartesianInitialCondition(sc.r_ii_m, sc.v_ii_m_s), source=sc.source)
        for sc in ics.spacecraft]
    return states, false
end

# --- configuration builders ------------------------------------------------

"One CYGNSS-class box on the given initial condition."
build_spacecraft(ic, id::Int) = make_three_body_spacecraft(
    bus_dims=CYGNSS_BUS_DIMS_M,
    panel_dims=(0.01, CYGNSS_ARRAY_HALF_SPAN_M, CYGNSS_ARRAY_CHORD_M),
    bus_mass=CYGNSS_BUS_MASS_KG,
    panel_mass_each=CYGNSS_ARRAY_MASS_EACH_KG,
    panel_offset_y=(CYGNSS_BUS_DIMS_M[2] + CYGNSS_ARRAY_HALF_SPAN_M) / 2,
    ic=ic, reflection_coefficient=CYGNSS_SRP_CR, prop_mass=0.0, id=id)

"""
    build_args(planet, spacecraft, initial_time, outdir; mission_time, num_steps_to_save)

The run configuration of record: EarthGGM05C 50x50, Sun and Moon third body,
solar radiation pressure and NRLMSISE-00 drag with real space-weather indices,
integrated with DP8. The reconstruction record (section 3) states that the
automatic stiff/nonstiff default costs 7.2 km on this scenario, so the solver
is pinned rather than left to the default.
"""
function build_args(planet, spacecraft::Vector, initial_time, outdir::AbstractString;
                    mission_time::Float64, num_steps_to_save::Int,
                    drag_scales::Vector{Float64}=fill(CYGNSS_DRAG_SCALE, length(spacecraft)))
    effectors = (
        GravitationalHarmonicsModel(50, 50, joinpath(HARMONICS_DIR, "EarthGGM05C.csv"), planet),
        NBodyGravityModel(body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet),
        SolarRadiationPressureModel(CYGNSS_SRP_CR, spacecraft[1].root.ref_area),
        ScaledAero(AerodynamicCoefficientfM(), copy(drag_scales)),
    )
    base = make_example_config(planet=planet, spacecraft=spacecraft[1], mission_time=mission_time,
        initial_time=initial_time, dynamic_effectors=effectors,
        density_model=NRLMSISE00AtmosphereModel(use_space_indices=true),
        orientation_sim=false, keplerian=true, EI_km=120.0, verbose=false, results=true,
        results_directory=String(outdir))
    return SM.SimulationConfiguration(
        file_paths=base.file_paths, simulation_settings=base.simulation_settings,
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true,
            number_of_orbits=1, mission_time=mission_time, orientation_sim=false,
            num_steps_to_save=num_steps_to_save, data_rate=5.0),
        environment_model=base.environment_model,
        dynamics_model=SM.DynamicsModel(spacecraft, effectors),
        guidance_model=base.guidance_model, navigation_model=base.navigation_model,
        control_model=base.control_model, initial_time=base.initial_time,
        integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-10,
            dt_max_orbit=30.0),
        solver_config=SM.SolverConfig(solver_mode=:dp8))
end

"""
    separation_rtn(df, k, series) -> (t_s, total, radial, along, cross, dx, dy, dz)

Separation of simulated spacecraft `k` from a flown track, scored on the run's
own saved times with the 1 Hz track interpolated onto them, resolved both in
the orbit frame (radial, along-track, cross-track) and in Cartesian components.
Saved times that fall inside a track gap longer than two seconds are dropped:
the product has outages of up to ten minutes, and interpolating across one
would report hundreds of kilometers of "error" that is entirely the gap.
"""
function separation_rtn(df::DataFrame, k::Int, series)
    saved_t = Float64.(df.time)
    sx = Float64.(df[!, "sc$(k)_pos_1"]); sy = Float64.(df[!, "sc$(k)_pos_2"]); sz = Float64.(df[!, "sc$(k)_pos_3"])
    vx = Float64.(df[!, "sc$(k)_vel_1"]); vy = Float64.(df[!, "sc$(k)_vel_2"]); vz = Float64.(df[!, "sc$(k)_vel_3"])
    t_tel = series.t_s; P = series.pos_m
    ts = Float64[]; tot = Float64[]; rad = Float64[]; alo = Float64[]; cro = Float64[]
    dx = Float64[]; dy = Float64[]; dz = Float64[]
    for i in eachindex(saved_t)
        t = saved_t[i]
        (t < t_tel[1] || t > t_tel[end]) && continue
        j = clamp(searchsortedfirst(t_tel, t), 2, length(t_tel))
        (t_tel[j] - t_tel[j - 1] > 2.0) && continue
        w = (t - t_tel[j - 1]) / (t_tel[j] - t_tel[j - 1])
        ref = SVector{3, Float64}(P[1, j - 1] + w * (P[1, j] - P[1, j - 1]),
                                  P[2, j - 1] + w * (P[2, j] - P[2, j - 1]),
                                  P[3, j - 1] + w * (P[3, j] - P[3, j - 1]))
        r = SVector{3, Float64}(sx[i], sy[i], sz[i]); v = SVector{3, Float64}(vx[i], vy[i], vz[i])
        d = r - ref
        R = normalize(r); N = normalize(cross(r, v)); T = cross(N, R)
        push!(ts, t); push!(tot, norm(d)); push!(rad, dot(d, R)); push!(alo, dot(d, T)); push!(cro, dot(d, N))
        push!(dx, d[1]); push!(dy, d[2]); push!(dz, d[3])
    end
    return (t_s=ts, total=tot, radial=rad, along=alo, cross=cro, dx=dx, dy=dy, dz=dz)
end

"Initial condition with the velocity magnitude rescaled to the vis-viva energy of `a_target`."
function ic_at_sma(r0::SVector{3, Float64}, v0::SVector{3, Float64}, a_target::Float64, planet)
    vmag = sqrt(planet.μ * (2.0 / norm(r0) - 1.0 / a_target))
    return SM.CartesianInitialCondition(r0, v0 .* (vmag / norm(v0)))
end

"""
    mean_sma_decay(t_s, a_m, period_s) -> (a_first, a_last, rate_m_per_s)

Secular decay of the orbit-averaged semimajor axis: the mean of the osculating
values over the first revolution, the mean over the last, and the slope between
them. Averaging over a whole revolution is what removes the short-period
variation; a partial revolution biases the mean instead.
"""
function mean_sma_decay(t_s::Vector{Float64}, a_m::Vector{Float64}, period_s::Float64)
    first_window = [a_m[k] for k in eachindex(t_s) if t_s[k] <= t_s[1] + period_s]
    last_window = [a_m[k] for k in eachindex(t_s) if t_s[k] >= t_s[end] - period_s]
    (length(first_window) < 2 || length(last_window) < 2) && return (NaN, NaN, NaN)
    a0 = mean(first_window); a1 = mean(last_window)
    span = (t_s[end] - period_s / 2) - (t_s[1] + period_s / 2)
    return (a0, a1, (a1 - a0) / span)
end

"Osculating semimajor axis along a telemetry series."
function series_sma(series, planet; stride::Int=10)
    idx = 1:stride:length(series.t_s)
    a = [1.0 / (2.0 / norm(SVector{3, Float64}(series.pos_m[:, k])) -
                dot(SVector{3, Float64}(series.vel_mps[:, k]), SVector{3, Float64}(series.vel_mps[:, k])) / planet.μ)
         for k in idx]
    return collect(Float64, series.t_s[idx]), a
end

"Osculating semimajor axis of simulated spacecraft `k` over the saved rows."
function run_sma(df::DataFrame, k::Int, planet; stride::Int=1)
    t = Float64.(df.time)[1:stride:end]
    x = Float64.(df[!, "sc$(k)_pos_1"])[1:stride:end]; y = Float64.(df[!, "sc$(k)_pos_2"])[1:stride:end]; z = Float64.(df[!, "sc$(k)_pos_3"])[1:stride:end]
    vx = Float64.(df[!, "sc$(k)_vel_1"])[1:stride:end]; vy = Float64.(df[!, "sc$(k)_vel_2"])[1:stride:end]; vz = Float64.(df[!, "sc$(k)_vel_3"])[1:stride:end]
    a = [1.0 / (2.0 / norm(SVector{3, Float64}(x[j], y[j], z[j])) -
                (vx[j]^2 + vy[j]^2 + vz[j]^2) / planet.μ) for j in eachindex(t)]
    return t, a
end

"""
    fit_states!(states, planet, initial_time, et0) -> (states, drag_scales)

Fit two scalars per spacecraft to that spacecraft's own flown track over the
window: the magnitude of its initial velocity (equivalently its initial orbital
energy) and its effective drag scale. Position and velocity DIRECTION are left
exactly as the initial-conditions file gives them.

Every spacecraft on this page is flown, so every one of them is fitted against
its own telemetry and every separation this demo reports is an in-sample fit.
That is a change in what the number means, not in how it is produced: when only
two spacecraft had telemetry, the other five were shown unfitted and their
agreement meant nothing at all.

Why two. The Level 1 velocity is quantized at 1 m/s and the initial-conditions
file fits it from the position arc rather than reading it; whatever residual
error is left, the repository's own CYGNSS loader records that 0.6 m/s along
track is about a kilometer of semimajor axis and roughly 120 km of along-track
drift at 48 hours (`test/gmat_scenario_matrix.jl`,
`_build_cygnss_cyg04_96hr_inertial_reference`). That is the first scalar. The
second is drag: a free-molecular coefficient on a generic published box is not
this vehicle's ballistic coefficient, and NRLMSISE-00 is not the atmosphere
that actually flew, so the simulated orbit decays at the wrong rate. An initial
energy alone can null the along-track error at one instant, but a wrong decay
rate then bows the error by tens of kilometers in between. The reconstruction
record carries the same two quantities (section 3: "calibrated effective drag
scale about 0.3"; section 2.1: initial states "fitted from position over an arc
rather than taken from a single quantized velocity sample").

How. The drag scale is set by matching the secular decay of the orbit-averaged
semimajor axis to the telemetry's own, which is a measurement rather than a
search. The initial energy is then a Newton step through the two-body relation
`s = -1.5 r n t (da / a)`, solved on the MEAN along-track separation over the
window rather than on its last sample, so the error ends up centered on zero
instead of pinned to zero at the instant this demo reports.

This makes the reported separations an IN-SAMPLE fit over the scored window,
not a prediction, and they are reported as such.
"""
function fit_states!(states::Vector, planet, initial_time, et0::Float64;
                     iterations::Int=4, tolerance_m::Float64=1_000.0)
    n_sc = length(states)
    idx = [k for (k, s) in enumerate(states) if has_track(s.name)]
    isempty(idx) && return states, fill(CYGNSS_DRAG_SCALE, n_sc)
    fitdir = joinpath(OUTDIR, "fit"); mkpath(fitdir)
    cache_path = joinpath(fitdir, "fitted_parameters.json")

    r0 = Dict(k => SVector{3, Float64}(states[k].r_ii_m) for k in idx)
    v0 = Dict(k => SVector{3, Float64}(states[k].v_ii_m_s) for k in idx)

    function apply!(a::Dict{Int, Float64}, c::Dict{Int, Float64})
        for k in idx
            states[k] = merge(states[k], (ic=ic_at_sma(r0[k], v0[k], a[k], planet),
                fitted_sma_m=a[k], fitted_drag_scale=c[k],
                provenance=states[k].provenance * ", fitted"))
        end
        mean_scale = mean(c[k] for k in idx)
        return [k in idx ? c[k] : mean_scale for k in 1:n_sc]
    end

    # A rerun (re-exporting the page from results already on disk) reuses the
    # fitted values instead of repeating the fit propagations.
    if isfile(cache_path) && get(ENV, "SPACEAGORA_DEMO_FORCE", "0") != "1"
        cached = JSON.parsefile(cache_path)
        if all(k -> haskey(cached, states[k].name), idx)
            println("  reusing the fitted parameters in ", cache_path)
            a = Dict(k => Float64(cached[states[k].name]["sma_m"]) for k in idx)
            c = Dict(k => Float64(cached[states[k].name]["drag_scale"]) for k in idx)
            return states, apply!(a, c)
        end
    end

    a = Dict{Int, Float64}(k => TV.rvtoorbitalelement(r0[k], v0[k], planet)[1] for k in idx)
    c = Dict{Int, Float64}(k => CYGNSS_DRAG_SCALE for k in idx)
    for k in idx
        @printf("  %s: starting semimajor axis %.3f m (from the initial-conditions file), drag scale %.3f\n",
            states[k].name, a[k], c[k])
    end

    for iter in 1:iterations
        sc = [build_spacecraft(ic_at_sma(r0[k], v0[k], a[k], planet), j) for (j, k) in enumerate(idx)]
        args = build_args(planet, sc, initial_time, fitdir; mission_time=WINDOW_S,
            num_steps_to_save=4000, drag_scales=[c[k] for k in idx])
        rm(joinpath(fitdir, "simulation_results.feather"); force=true)
        run_simulation(args)
        fdf = DataFrame(Arrow.Table(joinpath(fitdir, "simulation_results.feather")))
        worst = 0.0
        for (j, k) in enumerate(idx)
            series = track_series(states[k].name)
            sep = separation_rtn(fdf, j, series)
            period_s = 2pi * sqrt(a[k]^3 / planet.μ)
            t_tel, a_tel = series_sma(series, planet)
            t_sim, a_sim = run_sma(fdf, j, planet; stride=4)
            _, _, rate_tel = mean_sma_decay(t_tel, a_tel, period_s)
            _, _, rate_sim = mean_sma_decay(t_sim, a_sim, period_s)
            # The decay measurement needs an arc long enough for the secular
            # trend to stand clear of the short-period variation of the
            # osculating semimajor axis. Over the 96 h window that is twenty
            # revolutions either side; over a shortened smoke window it is not,
            # and measuring anyway returns a decay of the wrong sign and drives
            # the scale into its clamp. Leave the drag alone in that case.
            long_enough = sep.t_s[end] >= 40 * period_s
            ratio = (long_enough && isfinite(rate_tel) && isfinite(rate_sim) && rate_sim < 0.0 && rate_tel < 0.0) ?
                clamp(rate_tel / rate_sim, 0.2, 5.0) : 1.0
            n_mean = sqrt(planet.μ / a[k]^3)
            # Center the along-track error over the window rather than null it
            # at the end. A semimajor-axis error drifts the along-track
            # position as `s = -1.5 r n t (da / a)`, whose mean over [0, T] is
            # half its end value; solving on the mean therefore leaves the
            # error symmetric about zero instead of hanging the whole excursion
            # in the middle of the window and pinning the last sample to zero.
            # The last sample is the number this demo reports, so it must not
            # be the thing the fit drives to zero by construction.
            da = 2.0 * mean(sep.along) * a[k] / (1.5 * norm(r0[k]) * n_mean * sep.t_s[end])
            rms = sqrt(mean(sep.total .^ 2))
            @printf("  %s pass %d%s: rms %.3f km, end %.3f km, mean along-track %+.3f km; SMA decay sim %+.4f vs telemetry %+.4f mm/s -> drag scale %.3f, da %+.2f m\n",
                states[k].name, iter, iter == 1 ? " (unfitted)" : "", rms / 1e3, sep.total[end] / 1e3,
                mean(sep.along) / 1e3, rate_sim * 1e3, rate_tel * 1e3,
                clamp(c[k] * ratio, 0.02, 50.0), da)
            c[k] = clamp(c[k] * ratio, 0.02, 50.0)
            a[k] += da
            worst = max(worst, rms)
        end
        worst <= tolerance_m && break
    end

    open(cache_path, "w") do io
        JSON.print(io, Dict(states[k].name => Dict("sma_m" => a[k], "drag_scale" => c[k]) for k in idx))
    end
    return states, apply!(a, c)
end

# --- run -------------------------------------------------------------------

planet = Earth("", SPICE_PATH)
initial_time = SM.InitialTime(year=2025, month=6, day=6, hour=0, minute=0, second=0.0)
et0 = et_of(EPOCH_UTC)

states, stubbed = load_states(planet)
println("initial conditions: ", stubbed ? "NOMINAL STUB" : IC_PATH)
for s in states
    r = SVector{3, Float64}(s.r_ii_m); v = SVector{3, Float64}(s.v_ii_m_s)
    oe = TV.rvtoorbitalelement(r, v, planet)
    @printf("  %-12s norad %5d  %-10s alt %6.1f km  i %5.2f deg  RAAN %6.2f deg  u %6.2f deg  T %5.1f min\n",
        s.name, s.norad_id, s.provenance, (oe[1] - planet.Rp_e) / 1e3, rad2deg(oe[3]), rad2deg(oe[4]),
        mod(rad2deg(oe[5] + oe[6]), 360.0), 2pi * sqrt(oe[1]^3 / planet.μ) / 60)
end

init_nrlmsise_space_indices!()

drag_scales = fill(CYGNSS_DRAG_SCALE, length(states))
if FIT_SMA
    println("fitting the initial energy and the drag scale of every spacecraft with a flown track:")
    states, drag_scales = fit_states!(states, planet, initial_time, et0)
    println("drag scales in the run: ", join([@sprintf("%s %.3f", states[k].name, drag_scales[k]) for k in eachindex(states)], ", "))
end

spacecraft = [build_spacecraft(s.ic, k) for (k, s) in enumerate(states)]
args = build_args(planet, spacecraft, initial_time, OUTDIR; mission_time=WINDOW_S,
    num_steps_to_save=6000, drag_scales=drag_scales)

prefix = run_or_reuse!(args, OUTDIR)

# --- the page's own copy of the run ----------------------------------------
# The page is exported from a trimmed copy of the results, not from the run
# itself, for two reasons.
#
# Provenance. The scene sidecar names spacecraft `sc<id>`; the viewer draws
# that name as the marker label and shows it in the selection panel. When the
# constellation is a mixture — some flown, some catalogue, some nominal — each
# spacecraft's own word goes in its name, because a reader must not have to
# open the source to find out which is which. When they all have the same
# provenance, as they do now that every state is flown, repeating the word
# seven times says nothing and only crowds the labels: the page footer says it
# once instead.
#
# Confidentiality. Spacecraft mass properties are not published here (see the
# header: the run's geometry is a generic box from public figures and no
# restricted mission-configuration file is read). Per-spacecraft mass is a
# default saved column that the viewer shows as a panel row, so the
# `sc<i>_mass` columns are dropped from the page's copy, and each link's
# `mass_kg` is zeroed in the scene, which nothing in the viewer reads. The link
# box DIMENSIONS stay, because the viewer needs them to draw the generic box at
# all; they are the published body envelope, not a reconstruction.
const mixed_provenance = length(unique(first(split(s.provenance, ',')) for s in states)) > 1
const PAGE_DIR = joinpath(OUTDIR, "page")
mkpath(PAGE_DIR)
page_prefix = joinpath(PAGE_DIR, "simulation_results")
let doc = JSON.parsefile(prefix * "_scene.json")
    for (k, s) in enumerate(states)
        k <= length(doc["spacecraft"]) || break
        doc["spacecraft"][k]["name"] = mixed_provenance ? "$(s.name) · $(s.provenance)" : s.name
        for link in doc["spacecraft"][k]["links"]
            link["mass_kg"] = 0.0
        end
    end
    open(page_prefix * "_scene.json", "w") do io
        JSON.print(io, doc)
    end
    full = DataFrame(Arrow.Table(prefix * ".feather"))
    dropped = [c for c in names(full) if occursin(r"^sc\d+_mass$", c)]
    Arrow.write(page_prefix * ".feather", select(full, Not(dropped)))
    println("page copy: dropped ", length(dropped), " mass columns")
end

# --- telemetry ghosts and the separation that scores the run ---------------
df = DataFrame(Arrow.Table(prefix * ".feather"))
saved_t = Float64.(df.time)
println("rows=", nrow(df), "  span=", round(saved_t[end] / 3600; digits=2), " h")

"""
    spacecraft_hex_color(index, count) -> String

The color the viewer gives spacecraft `index` of `count`, as a hex string, so a
ghost can be drawn in its own twin's color. Mirrors `spacecraftColor` in
`viewer/src/spacecraft.js`: hue spread evenly around the wheel, HSL saturation
0.85 and lightness 0.6.
"""
function spacecraft_hex_color(index::Int, count::Int)::String
    h = count <= 1 ? 0.12 : mod((index - 1) / count, 1.0)
    sat, light = 0.85, 0.6
    c = (1 - abs(2 * light - 1)) * sat
    x = c * (1 - abs(mod(h * 6, 2) - 1))
    m = light - c / 2
    r, g, b = h < 1/6 ? (c, x, 0.0) : h < 2/6 ? (x, c, 0.0) : h < 3/6 ? (0.0, c, x) :
              h < 4/6 ? (0.0, x, c) : h < 5/6 ? (x, 0.0, c) : (c, 0.0, x)
    to255(v) = clamp(round(Int, (v + m) * 255), 0, 255)
    return "#" * string(to255(r); base=16, pad=2) * string(to255(g); base=16, pad=2) * string(to255(b); base=16, pad=2)
end

references = Any[]
separations = Any[]
for (k, s_) in enumerate(states)
    has_track(s_.name) || continue
    series = track_series(s_.name)
    label = s_.name

    # Ghost: one sample every `stride` rows of the flown track. The record is
    # 345,600 samples per spacecraft and there are seven of them; a couple of
    # thousand each is indistinguishable at globe scale and keeps the page
    # inside its size budget. Each ghost takes its twin's own color, so the
    # translucent copy beside a spacecraft is unmistakably that spacecraft's.
    stride = max(1, cld(length(series.t_s), 2_500))
    idx = [j for j in 1:stride:length(series.t_s) if series.t_s[j] <= saved_t[end] + 1.0]
    println("ghost ", label, ": ", length(idx), " samples, ", series.note)
    push!(references, (name="$(label) flown track", t_s=series.t_s[idx], pos_m=series.pos_m[:, idx],
        vel_mps=series.vel_mps[:, idx], target=k,
        color=spacecraft_hex_color(k, length(states)), opacity=0.45))

    # Separation, scored by `separation_rtn` on the run's own saved times.
    sep = separation_rtn(df, k, series)
    at(hours) = (i = findlast(t -> t <= hours * 3600.0, sep.t_s); i === nothing ? NaN : sep.total[i])
    @printf("separation %s vs its own flown track: 0 h %.1f m | 24 h %.3f km | 48 h %.3f km | 72 h %.3f km | end %.3f km (max %.3f km, rms %.3f km)\n",
        label, sep.total[1], at(24) / 1e3, at(48) / 1e3, at(72) / 1e3, sep.total[end] / 1e3,
        maximum(sep.total) / 1e3, sqrt(mean(sep.total .^ 2)) / 1e3)
    @printf("  %s end-of-window components: radial %+.3f km  along-track %+.3f km  cross-track %+.3f km; first 48 h mean per-axis rms %.3f km\n",
        label, sep.radial[end] / 1e3, sep.along[end] / 1e3, sep.cross[end] / 1e3,
        # Mean per-axis RMSE over the first 48 hours, in Cartesian components:
        # the form the reconstruction record reports (0.967 km for FM4), so the
        # two can be read against each other directly.
        let h = findall(t -> t <= 48 * 3600.0, sep.t_s)
            isempty(h) ? NaN : mean([sqrt(mean(getproperty(sep, c)[h] .^ 2)) for c in (:dx, :dy, :dz)]) / 1e3
        end)
    push!(separations, (label=label, sep=sep))
end

# --- page ------------------------------------------------------------------
# Frame budget. 96 h of seven spacecraft is 69,121 saved rows; 3000 embedded
# frames is one sample every 115 s, about 50 per orbit, which the trails and
# the ground tracks both read well and which keeps the page inside its size
# budget alongside a 4k texture and seven 2500-sample flown-track ghosts. The
# planet-fixed frame opens on the picture the ground tracks belong to.
html = export_visualization(page_prefix; max_frames=3000, trail_orbits=1, texture_resolution="4k",
    frame=:planet_fixed, ground_tracks=true,
    title="AGORA CYGNSS · the constellation over its telemetry window",
    references=references)
println("html: ", html, " ", filesize(html))

# Elements of the flown constellation, for the page header: this is what drag
# has left of the design orbit, not the design orbit.
alts = Float64[]; incs = Float64[]; periods = Float64[]
for s in states
    oe = TV.rvtoorbitalelement(SVector{3, Float64}(s.r_ii_m), SVector{3, Float64}(s.v_ii_m_s), planet)
    push!(alts, (oe[1] - planet.Rp_e) / 1e3); push!(incs, rad2deg(oe[3]))
    push!(periods, 2pi * sqrt(oe[1]^3 / planet.μ) / 60)
end
orbit_note = @sprintf("%.0f-%.0f km, i %.2f-%.2f°, %.1f min", minimum(alts), maximum(alts),
    minimum(incs), maximum(incs), mean(periods)) * (stubbed ? " (development stub)" : "")

all_flown = all(s -> startswith(s.provenance, "telemetry"), states)
worst = isempty(separations) ? nothing : separations[argmax([x.sep.total[end] for x in separations])]
best = isempty(separations) ? nothing : separations[argmin([x.sep.total[end] for x in separations])]
separation_note = isempty(separations) ? "" : @sprintf(
    "At the end of the 96 hours the simulated spacecraft sit between %.2f km (%s) and %.2f km (%s) from their own flown tracks, and the root-mean-square separation over the whole window runs from %.2f to %.2f km. ",
    best.sep.total[end] / 1e3, best.label, worst.sep.total[end] / 1e3, worst.label,
    minimum([sqrt(mean(x.sep.total .^ 2)) for x in separations]) / 1e3,
    maximum([sqrt(mean(x.sep.total .^ 2)) for x in separations]) / 1e3)

cdn = build_cdn_page(html, joinpath(OUTDIR, "artifact.html"), "AGORA CYGNSS Constellation",
    "AGORA CYGNSS · $(length(states)) observatories, 2025-06-06 to 2025-06-09",
    orbit_note, "96 hours (≈62 orbits)",
    "$(length(states)) spacecraft, not the eight that launched: CYGNSS FM06 (NORAD 41889) is absent because it was no longer in orbit. " *
    "NASA lost contact with it in November 2022 and the satellite catalogue gives its decay date as 2024-06-13. " *
    (all_flown ?
        "All $(length(states)) are flown states. Every initial condition and every translucent ghost here comes from that spacecraft's own NASA Level 1 navigation solution over this window — nothing is a catalogue propagation and nothing is a design orbit. " :
        "Each spacecraft's label says where its initial state came from; only the telemetry-backed ones are flown states. ") *
    (FIT_SMA ? "Two scalars per spacecraft, the magnitude of its initial velocity and an effective drag scale, were fitted to that spacecraft's own track over this window, so the agreement between a spacecraft and its ghost is an in-sample fit and not a prediction. " : "") *
    "Click a solid spacecraft to read its separation from its own flown track. " *
    separation_note *
    "The spacecraft is drawn as a generic small-satellite box, not as a reconstruction of the flight geometry. " *
    "\"Ground tracks\" draws each sub-satellite point on the surface. Drag to orbit, wheel to zoom, Space to pause, F to follow.")
println("cdn: ", cdn, " ", filesize(cdn))
