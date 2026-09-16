# The CYGNSS constellation flown together over the window its flight telemetry
# covers: 2025-06-06T00:00:00Z to 2025-06-09T24:00:00Z, 96 hours.
#
# Seven spacecraft, not the eight that were launched: CYGNSS FM06 (NORAD 41889)
# was no longer in orbit at this epoch. NASA lost contact with it in November
# 2022 and the satellite catalogue gives its decay date as 2024-06-13. The
# initial-conditions file's notes carry the citations.
#
# The constellation has no propulsion, so by June 2025 drag had taken it well
# below its design orbit: 438 to 446 km, inclination 34.88 to 34.97 degrees,
# periods of 93.4 to 93.5 minutes, against the 510 km and 95 minutes of the
# design. Nothing here sanity-checks a flown state against the design values.
#
# Initial states come from `data/telemetry/CYGNSS/constellation_ics_20250606.json`
# (gitignored; built by the constellation initial-conditions script). Every
# spacecraft in that file carries a `provenance` field — `telemetry`,
# `catalogue` or `nominal` — and this demo puts that word in the spacecraft's
# own viewer label, because only the telemetry-backed ones are flown states
# and a reader must not have to open the source to find out which.
#
# When the file is absent the demo falls back to a development stub: FM1 and
# FM4 from their own telemetry and the rest of the launched eight on the
# published CYGNSS design orbit (510 km, 35 degrees, one plane, 45 degrees
# apart in argument of latitude — Spaceflight101 CYGNSS orbit design;
# NASA/eoPortal CYGNSS mission page), marked `nominal`. A stub run is not a
# reconstruction of the flown constellation and the page says so. The stub is
# a development convenience only: it does not know that FM06 is gone, and its
# design-orbit spacecraft sit about 70 km above where the real ones fly.
#
# Force model, following `docs/spaceagora_cygnss_reconstruction_record.md`
# section 3: EarthGGM05C gravity to degree and order 50, Sun and Moon third
# body, solar radiation pressure, and NRLMSISE-00 drag with real CelesTrak
# space-weather indices.
#
# Flight telemetry of FM1 and FM4 (`cyg01/cyg04_nasa_pvt_96hr.feather`) is
# attached as reference ghosts, so the simulated orbit can be read against
# the flight solution directly in the page.
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
using .CygnssICs: load_constellation_ics
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
# Unix seconds of EPOCH_UTC (2025-06-06T00:00:00Z), the time base the telemetry
# files' own `time` column and the run's elapsed time both count from.
const EPOCH_UNIX_S = 1_749_168_000.0

# Telemetry files, by the NORAD id of the spacecraft they belong to. Both span
# the window at about 1 Hz; FM4 carries precomputed J2000 columns, FM1 only the
# WGS84 Earth-fixed PVT solution (converted below through the same SPICE
# ITRF93 -> J2000 path the repository's own 96 h loader uses).
const TELEMETRY_FILES = Dict(
    41887 => (file="cyg01_nasa_pvt_96hr.feather", label="CYGNSS FM01"),
    41885 => (file="cyg04_nasa_pvt_96hr.feather", label="CYGNSS FM04"),
)

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

# Fit the initial orbital energy and the effective drag scale of the
# telemetry-backed spacecraft to their own telemetry (see `fit_states!`).
# SPACEAGORA_DEMO_CYGNSS_FIT_SMA=0 propagates the initial-conditions file's
# states untouched, at the uncalibrated drag scale, instead.
const FIT_SMA = get(ENV, "SPACEAGORA_DEMO_CYGNSS_FIT_SMA", "1") != "0"

# Published CYGNSS design orbit, used only by the nominal stub below.
const DESIGN_ALTITUDE_M = 510.0e3          # NASA/eoPortal, Spaceflight101
const DESIGN_INCLINATION_DEG = 35.0        # NASA/eoPortal, Spaceflight101
const DESIGN_SPACING_DEG = 45.0            # Spaceflight101 CYGNSS orbit design: eight spacecraft, one plane, 45 deg apart

const FM_NAMES = ["CYGNSS FM01", "CYGNSS FM02", "CYGNSS FM03", "CYGNSS FM04",
                  "CYGNSS FM05", "CYGNSS FM06", "CYGNSS FM07", "CYGNSS FM08"]
# Catalogue identities verified against the public catalogue by the
# initial-conditions agent; 41889 (FM06) returns no current element set.
const FM_NORAD = [41887, 41886, 41891, 41885, 41884, 41889, 41890, 41888]

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

"""
    drop_invalid_fixes(t, pos, vel; tolerance_m=50e3) -> (t, pos, vel, dropped)

Drop navigation fixes whose radius is more than `tolerance_m` from the record's
median radius. Not every PVT row is a valid fix: the reconstruction record
counts 2,479 flagged samples in the FM4 arc, 1.4 percent of it, including
zero-satellite solutions. A near-circular 440 km orbit varies by about five
kilometers in radius over a revolution, so a 50 km band keeps every real sample
and removes fixes that place the spacecraft hundreds of kilometers off the
orbit (the FM1 file carries a burst of these, the worst 490 km out).
"""
function drop_invalid_fixes(t::Vector{Float64}, pos::Matrix{Float64}, vel::Matrix{Float64};
                            tolerance_m::Float64=50.0e3)
    r = [norm(SVector{3, Float64}(pos[:, k])) for k in axes(pos, 2)]
    r_med = median(r)
    keep = findall(k -> abs(r[k] - r_med) <= tolerance_m, eachindex(r))
    return t[keep], pos[:, keep], vel[:, keep], length(r) - length(keep)
end

"""
    telemetry_series(file, et0) -> (t_s, pos_m, vel_mps)

The flight PVT record as J2000 Earth-centered state, on elapsed seconds from
the run epoch. `pos_ii_*`/`vel_ii_*` are used when the file carries them
(FM4); otherwise the WGS84 Earth-fixed PVT columns are rotated through SPICE,
exactly as `test/gmat_scenario_matrix.jl` does for the same product.
"""
function telemetry_series(file::AbstractString, et0::Float64; stride::Int=1)
    df = DataFrame(Arrow.Table(joinpath(TELEMETRY_DIR, String(file))))
    t_unix = Float64.(df[!, "pvt_unix_seconds"])
    perm = sortperm(t_unix)
    t_unix = t_unix[perm]
    # Elapsed seconds from the RUN epoch, not from the file's own first sample:
    # the first PVT row sits 1 s after 2025-06-06T00:00:00Z, and a 1 s error in
    # the Earth-fixed to J2000 rotation is about 500 m of position and 0.55 m/s
    # of velocity at this orbit. (Checked against the FM4 file's own pos_ii
    # columns: with this time base the two agree to a few centimeters.)
    t_rel = t_unix .- EPOCH_UNIX_S
    keep = 1:stride:length(t_rel)
    t_rel = t_rel[keep]

    if "pos_ii_1" in names(df)
        pos = permutedims(hcat(Float64.(df[!, "pos_ii_1"])[perm][keep],
                               Float64.(df[!, "pos_ii_2"])[perm][keep],
                               Float64.(df[!, "pos_ii_3"])[perm][keep]))
        vel = permutedims(hcat(Float64.(df[!, "vel_ii_1"])[perm][keep],
                               Float64.(df[!, "vel_ii_2"])[perm][keep],
                               Float64.(df[!, "vel_ii_3"])[perm][keep]))
        t_rel, pos, vel, dropped = drop_invalid_fixes(t_rel, pos, vel)
        return t_rel, pos, vel, "pos_ii/vel_ii (J2000 columns in the file), $(dropped) invalid fixes dropped"
    end

    TV._planet_from_name("earth")   # furnish the leap-second and orientation kernels
    xe = Float64.(df[!, "sc_pos_x_pvt_m"])[perm][keep]
    ye = Float64.(df[!, "sc_pos_y_pvt_m"])[perm][keep]
    ze = Float64.(df[!, "sc_pos_z_pvt_m"])[perm][keep]
    vxe = Float64.(df[!, "sc_vel_x_pvt_mps"])[perm][keep]
    vye = Float64.(df[!, "sc_vel_y_pvt_mps"])[perm][keep]
    vze = Float64.(df[!, "sc_vel_z_pvt_mps"])[perm][keep]
    n = length(t_rel)
    pos = zeros(3, n); vel = zeros(3, n)
    for k in 1:n
        r, v = TV._planet_fixed_to_j2000_state("earth", et0 + t_rel[k],
            SVector{3, Float64}(xe[k], ye[k], ze[k]), SVector{3, Float64}(vxe[k], vye[k], vze[k]))
        pos[:, k] .= r; vel[:, k] .= v
    end
    t_rel, pos, vel, dropped = drop_invalid_fixes(t_rel, pos, vel)
    return t_rel, pos, vel, "sc_pos/sc_vel_pvt rotated ITRF93 -> J2000 through SPICE, $(dropped) invalid fixes dropped"
end

const _TELEMETRY_CACHE = Dict{Int, Any}()

"Full-rate `telemetry_series` for a spacecraft, computed once per process."
function telemetry_for(norad::Int, et0::Float64)
    get!(_TELEMETRY_CACHE, norad) do
        entry = TELEMETRY_FILES[norad]
        t, pos, vel, note = telemetry_series(entry.file, et0; stride=1)
        (t_s=t, pos_m=pos, vel_mps=vel, note=note, label=entry.label)
    end
end

# --- initial conditions ----------------------------------------------------

"""
    nominal_stub_states(planet, et0) -> Vector

The development stub used until `constellation_ics_20250606.json` exists, in
that file's schema. FM1 and FM4 are taken from their own flight telemetry at
the common epoch (`provenance = "telemetry"`); the other six are placed on the
published CYGNSS design orbit — 510 km circular, 35 degrees, one plane, 45
degrees apart in argument of latitude (NASA/eoPortal CYGNSS mission page;
Spaceflight101 CYGNSS orbit design) — and marked `provenance = "nominal"`.
The nominal six are NOT flown states and the page labels them as such. The
plane of the nominal six is FM4's plane at the epoch, which is itself only one
of the constellation's planes: FM1 and FM4 are about 24 degrees apart in right
ascension of the ascending node in this window.
"""
function nominal_stub_states(planet, et0::Float64)
    telemetry = Dict{Int, Any}()
    for norad in keys(TELEMETRY_FILES)
        series = telemetry_for(norad, et0)
        t_rel, pos, vel, note = series.t_s, series.pos_m, series.vel_mps, series.note
        # The first PVT row is one second after the common epoch. Step it back
        # to the epoch with a two-body Taylor step; over one second the terms
        # this drops (J2 and below) are under a centimeter.
        dt = t_rel[1]
        r1 = SVector{3, Float64}(pos[:, 1]); v1 = SVector{3, Float64}(vel[:, 1])
        acc = -planet.μ * r1 / norm(r1)^3
        r0 = r1 - v1 * dt + 0.5 * acc * dt^2
        v0 = v1 - acc * dt
        telemetry[norad] = (r=r0, v=v0, dt=dt, note=note)
    end

    fm4 = telemetry[41885]
    oe = TV.rvtoorbitalelement(fm4.r, fm4.v, planet)
    raan_deg = rad2deg(oe[4])
    u0_deg = rad2deg(oe[5] + oe[6])          # argument of latitude of FM4 at the epoch
    a_m = planet.Rp_e + DESIGN_ALTITUDE_M

    out = Any[]
    for (k, name) in enumerate(FM_NAMES)
        norad = FM_NORAD[k]
        if haskey(telemetry, norad)
            tel = telemetry[norad]
            push!(out, (name=name, norad_id=norad, r_ii_m=collect(tel.r), v_ii_m_s=collect(tel.v),
                provenance="telemetry", ic=SM.CartesianInitialCondition(tel.r, tel.v),
                source="$(TELEMETRY_FILES[norad].file) row 1, $(tel.note), stepped back $(tel.dt) s to the epoch"))
        else
            u = mod(u0_deg + (k - 1) * DESIGN_SPACING_DEG, 360.0)
            ic = SM.InitialCondition(ra=a_m, rp=a_m, i=DESIGN_INCLINATION_DEG, ω=0.0, Ω=raan_deg, ν=u)
            r, v = TV.orbitalelemtorv(ic, planet)
            push!(out, (name=name, norad_id=norad, r_ii_m=collect(Float64, r), v_ii_m_s=collect(Float64, v),
                provenance="nominal", ic=ic,
                source="published design orbit: $(DESIGN_ALTITUDE_M / 1e3) km circular, $(DESIGN_INCLINATION_DEG) deg, $(DESIGN_SPACING_DEG) deg spacing in argument of latitude (NASA/eoPortal CYGNSS; Spaceflight101 CYGNSS orbit design); plane from the FM4 telemetry at the epoch"))
        end
    end
    return out
end

"""
    load_states(planet, et0) -> (states, stubbed::Bool)

The constellation initial conditions, from the JSON file when it exists and
from the nominal stub when it does not.
"""
function load_states(planet, et0::Float64)
    if !isfile(IC_PATH)
        println("!! ", IC_PATH, " is absent: falling back to the DEVELOPMENT STUB.")
        println("!! FM1 and FM4 come from their telemetry; the rest are the published design orbit, not flown states.")
        return nominal_stub_states(planet, et0), true
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
    separation_rtn(df, k, series) -> (t_s, total, radial, along, cross)

Separation of simulated spacecraft `k` from a telemetry series, scored on the
run's own saved times with the 1 Hz telemetry interpolated onto them. Saved
times that fall inside a telemetry gap longer than two seconds are dropped:
the record has outages of up to ten minutes, and interpolating across one
would report hundreds of kilometers of "error" that is entirely the gap.
"""
function separation_rtn(df::DataFrame, k::Int, series)
    saved_t = Float64.(df.time)
    sx = Float64.(df[!, "sc$(k)_pos_1"]); sy = Float64.(df[!, "sc$(k)_pos_2"]); sz = Float64.(df[!, "sc$(k)_pos_3"])
    vx = Float64.(df[!, "sc$(k)_vel_1"]); vy = Float64.(df[!, "sc$(k)_vel_2"]); vz = Float64.(df[!, "sc$(k)_vel_3"])
    t_tel = series.t_s; P = series.pos_m
    ts = Float64[]; tot = Float64[]; rad = Float64[]; alo = Float64[]; cro = Float64[]
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
    end
    return (t_s=ts, total=tot, radial=rad, along=alo, cross=cro)
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

Fit two scalars per telemetry-backed spacecraft to that spacecraft's own
telemetry over the window: the magnitude of its initial velocity (equivalently
its initial orbital energy) and its effective drag scale. Position and velocity
DIRECTION are left exactly as the initial-conditions file gives them.

Why two. The NASA PVT velocity is a GPS Doppler solution; the repository's own
CYGNSS loader records that an along-track error of about 0.6 m/s in it shifts
the semimajor axis by about a kilometer and produces roughly 120 km of
along-track drift at 48 hours (`test/gmat_scenario_matrix.jl`,
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
    idx = [k for (k, s) in enumerate(states) if haskey(TELEMETRY_FILES, s.norad_id)]
    isempty(idx) && return states, fill(CYGNSS_DRAG_SCALE, n_sc)
    fitdir = joinpath(OUTDIR, "fit"); mkpath(fitdir)
    cache_path = joinpath(fitdir, "fitted_parameters.json")

    r0 = Dict(k => SVector{3, Float64}(states[k].r_ii_m) for k in idx)
    v0 = Dict(k => SVector{3, Float64}(states[k].v_ii_m_s) for k in idx)

    function apply!(a::Dict{Int, Float64}, c::Dict{Int, Float64})
        for k in idx
            states[k] = merge(states[k], (ic=ic_at_sma(r0[k], v0[k], a[k], planet),
                fitted_sma_m=a[k], fitted_drag_scale=c[k],
                provenance=states[k].provenance * ", energy and drag fitted"))
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
            series = telemetry_for(states[k].norad_id, et0)
            sep = separation_rtn(fdf, j, series)
            period_s = 2pi * sqrt(a[k]^3 / planet.μ)
            t_tel, a_tel = series_sma(series, planet)
            t_sim, a_sim = run_sma(fdf, j, planet; stride=4)
            _, _, rate_tel = mean_sma_decay(t_tel, a_tel, period_s)
            _, _, rate_sim = mean_sma_decay(t_sim, a_sim, period_s)
            ratio = (isfinite(rate_tel) && isfinite(rate_sim) && rate_sim < 0.0 && rate_tel < 0.0) ?
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

states, stubbed = load_states(planet, et0)
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
    println("fitting the initial energy and the drag scale of the telemetry-backed spacecraft:")
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
# that name as the marker label and shows it in the selection panel. Putting
# the provenance word in the name is what makes it legible in the page without
# reading the source.
#
# Confidentiality. Spacecraft mass properties are not published here (see the
# header: the run's geometry is a generic box from public figures and no
# restricted mission-configuration file is read). Per-spacecraft mass is a
# default saved column that the viewer shows as a panel row, so the
# `sc<i>_mass` columns are dropped from the page's copy, and each link's
# `mass_kg` is zeroed in the scene, which nothing in the viewer reads. The link
# box DIMENSIONS stay, because the viewer needs them to draw the generic box at
# all; they are the published body envelope, not a reconstruction.
const PAGE_DIR = joinpath(OUTDIR, "page")
mkpath(PAGE_DIR)
page_prefix = joinpath(PAGE_DIR, "simulation_results")
let doc = JSON.parsefile(prefix * "_scene.json")
    for (k, s) in enumerate(states)
        k <= length(doc["spacecraft"]) || break
        doc["spacecraft"][k]["name"] = "$(s.name) · $(s.provenance)"
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

references = Any[]
separations = Any[]
for (k, s_) in enumerate(states)
    haskey(TELEMETRY_FILES, s_.norad_id) || continue
    series = telemetry_for(s_.norad_id, et0)
    label = series.label

    # Ghost: one sample every `stride` telemetry rows. The 1 Hz record is
    # 345,600 samples per spacecraft; a few thousand is indistinguishable at
    # globe scale and keeps the page inside its size budget.
    stride = max(1, cld(length(series.t_s), 4_000))
    idx = [j for j in 1:stride:length(series.t_s) if series.t_s[j] <= saved_t[end] + 1.0]
    println("ghost ", label, ": ", length(idx), " samples, ", series.note)
    push!(references, (name="$(label) flight PVT", t_s=series.t_s[idx], pos_m=series.pos_m[:, idx],
        vel_mps=series.vel_mps[:, idx], target=k,
        color=(s_.norad_id == 41887 ? "#ff8c69" : "#8cd7ff"), opacity=0.45))

    # Separation, scored by `separation_rtn` on the run's own saved times.
    sep = separation_rtn(df, k, series)
    at(hours) = (i = findlast(t -> t <= hours * 3600.0, sep.t_s); i === nothing ? NaN : sep.total[i])
    @printf("separation %s vs its own telemetry: 0 h %.1f m | 24 h %.3f km | 48 h %.3f km | 72 h %.3f km | end %.3f km (max %.3f km, rms %.3f km)\n",
        label, sep.total[1], at(24) / 1e3, at(48) / 1e3, at(72) / 1e3, sep.total[end] / 1e3,
        maximum(sep.total) / 1e3, sqrt(mean(sep.total .^ 2)) / 1e3)
    @printf("  end-of-window components: radial %+.3f km  along-track %+.3f km  cross-track %+.3f km\n",
        sep.radial[end] / 1e3, sep.along[end] / 1e3, sep.cross[end] / 1e3)
    push!(separations, (label=label, sep=sep))
end

# --- page ------------------------------------------------------------------
# Frame budget. 96 h of seven spacecraft is 69,121 saved rows; 3000 embedded
# frames is one sample every 115 s, about 50 per orbit, which the trails and
# the ground tracks both read well and which keeps the page inside its size
# budget alongside a 4k texture and two 4000-sample telemetry ghosts. The
# planet-fixed frame opens on the picture the ground tracks belong to.
html = export_visualization(page_prefix; max_frames=3000, trail_orbits=1, texture_resolution="4k",
    frame=:planet_fixed, ground_tracks=true,
    title="AGORA CYGNSS · the constellation over its telemetry window",
    references=references)
println("html: ", html, " ", filesize(html))

# Elements of the flown constellation, for the page header: this is what drag
# has left of the design orbit, not the design orbit.
alts = Float64[]; incs = Float64[]; raans = Float64[]; periods = Float64[]
for s in states
    oe = TV.rvtoorbitalelement(SVector{3, Float64}(s.r_ii_m), SVector{3, Float64}(s.v_ii_m_s), planet)
    push!(alts, (oe[1] - planet.Rp_e) / 1e3); push!(incs, rad2deg(oe[3]))
    push!(raans, mod(rad2deg(oe[4]), 360.0)); push!(periods, 2pi * sqrt(oe[1]^3 / planet.μ) / 60)
end
orbit_note = @sprintf("%.0f-%.0f km, i %.2f-%.2f°, %.1f min", minimum(alts), maximum(alts),
    minimum(incs), maximum(incs), mean(periods)) * (stubbed ? " (development stub)" : "")

provenance_counts = Dict{String, Int}()
for s in states
    # The label carries the fit as well; count by the source of the state itself.
    provenance_counts[first(split(s.provenance, ','))] = get(provenance_counts, first(split(s.provenance, ',')), 0) + 1
end
provenance_note = join(["$(v) $(k)" for (k, v) in sort(collect(provenance_counts); by=first)], " and ")
separation_note = isempty(separations) ? "" :
    "At the end of the 96 hours the simulated " *
    join(["$(x.label) is $(round(x.sep.total[end] / 1e3; digits=2)) km from its own telemetry" for x in separations], " and ") *
    @sprintf(" (root-mean-square over the window %s).",
        join([@sprintf("%.2f km", sqrt(mean(x.sep.total .^ 2)) / 1e3) for x in separations], " and "))

cdn = build_cdn_page(html, joinpath(OUTDIR, "artifact.html"), "AGORA CYGNSS Constellation",
    "AGORA CYGNSS · $(length(states)) observatories, 2025-06-06 to 2025-06-09",
    orbit_note, "96 hours (≈62 orbits)",
    "$(length(states)) spacecraft, not eight: CYGNSS FM06 (NORAD 41889) is absent because it was no longer in orbit. " *
    "NASA lost contact with it in November 2022 and the satellite catalogue gives its decay date as 2024-06-13. " *
    "Every remaining spacecraft's label says where its initial state came from — " * provenance_note * ". " *
    "Only the telemetry-backed ones are flown states. The catalogue ones are historical element sets propagated to the epoch: " *
    "their orbit planes are good to better than 0.01°, but their position along the orbit carries about 12 km of uncertainty, " *
    "measured by running the same element sets for the two spacecraft that do have flight states. " *
    (FIT_SMA ? "The telemetry-backed spacecraft carry \"energy and drag fitted\": two scalars each, the magnitude of the initial velocity and an effective drag scale, were fitted to their own telemetry over this window, so their agreement with the ghosts is an in-sample fit and not a prediction. " : "") *
    "The two translucent ghosts are the FM1 and FM4 flight position solutions; click a solid spacecraft to read its separation from its own telemetry. " *
    separation_note * " " *
    "The spacecraft is drawn as a generic small-satellite box, not as a reconstruction of the flight geometry. " *
    "\"Ground tracks\" draws each sub-satellite point on the surface. Drag to orbit, wheel to zoom, Space to pause, F to follow.")
println("cdn: ", cdn, " ", filesize(cdn))
