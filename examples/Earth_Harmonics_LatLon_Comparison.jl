include(joinpath(@__DIR__, "common.jl"))

using CSV
using DataFrames
using Plots
using LinearAlgebra

# Reference case for comparing SpaceAGORA against an external propagator (the
# MATLAB one) that seeds from a geodetic latitude/longitude and rotates the
# Earth underneath a spherical-harmonics gravity field. Both sides need the
# same three things to agree: the epoch, the geodetic start point, and the
# degree/order of the field. All three are in the block below.
#
# Repo-local: no GRAM, no SPICE, no licensed assets.

# ===========================================================================
# Case. Edit these.
# ===========================================================================

const INITIAL_LATITUDE_DEG = 0.0        # geodetic latitude of the start point
const INITIAL_LONGITUDE_DEG = 0.0       # planet-fixed longitude of the start point
const INITIAL_ALTITUDE_M = 500e3        # above the reference ellipsoid
const INCLINATION_DEG = 45.0            # orbit inclination through that point
const DESCENDING = false                # true selects the southbound crossing
const FLIGHT_PATH_ANGLE_DEG = 0.0       # 0 starts at an apsis
const ECCENTRICITY = 0.0                # 0 is circular; the start is periapsis

const HARMONICS_DEGREE = 40
const HARMONICS_ORDER = 40

const INITIAL_EPOCH = InitialTime(year=2026, month=9, day=8, hour=5, minute=0, second=0.0)
const MISSION_TIME_S = 86400.0          # one day
const DATA_RATE_S = 10.0                # seconds between saved rows

# Any name in `available_save_fields()`, plus `SaveField`s of your own. Every
# saved value is plotted: the ones named in PANEL_LABELS below get titled axes,
# anything else gets a generic panel per component.
const SAVED_VALUES = (:orbital_elements, :gravity_accel)

const OUTPUT_DIR = joinpath(REPO_ROOT, "output", "harmonics_latlon_comparison")

# ===========================================================================
# Plot labelling. One row per saved column: (title, y-axis label, unit scale).
# Add an entry to give a new saved value titled panels; leave it out and it
# still gets plotted, just with generic labels.
# ===========================================================================

const PANEL_LABELS = Dict{Symbol, Vector{Tuple{String, String, Float64}}}(
    :orbital_elements => [
        ("Semimajor Axis", "a (km)", 1e-3),
        ("Eccentricity", "e", 1.0),
        ("Inclination", "i (deg)", 1.0),
        ("RAAN", "Ω (deg)", 1.0),
        ("Argument of Periapsis", "ω (deg)", 1.0),
        ("True Anomaly", "ν (deg)", 1.0)
    ],
    :gravity_accel => [
        ("Gravity Acceleration X", "aₓ (m/s²)", 1.0),
        ("Gravity Acceleration Y", "a_y (m/s²)", 1.0),
        ("Gravity Acceleration Z", "a_z (m/s²)", 1.0)
    ],
    :aero_accel => [
        ("Aero Acceleration X", "aₓ (m/s²)", 1.0),
        ("Aero Acceleration Y", "a_y (m/s²)", 1.0),
        ("Aero Acceleration Z", "a_z (m/s²)", 1.0)
    ],
    :total_accel => [
        ("Total Acceleration X", "aₓ (m/s²)", 1.0),
        ("Total Acceleration Y", "a_y (m/s²)", 1.0),
        ("Total Acceleration Z", "a_z (m/s²)", 1.0)
    ]
)

# Columns a saved field actually wrote, in order: `sc1_<name>` for a scalar,
# `sc1_<name>_1..n` for a vector. Whichever shape the field has, this finds it.
function _field_columns(df::DataFrame, field::Symbol)::Vector{Symbol}
    present = propertynames(df)
    scalar = Symbol("sc1_$(field)")
    scalar in present && return [scalar]
    columns = Symbol[]
    k = 1
    while true
        column = Symbol("sc1_$(field)_$(k)")
        column in present || break
        push!(columns, column)
        k += 1
    end
    return columns
end

function _panel(df::DataFrame, t_hr::Vector{Float64}, field::Symbol, column::Symbol, index::Int)
    labels = get(PANEL_LABELS, field, Tuple{String, String, Float64}[])
    title, ylabel, scale = index <= length(labels) ?
        labels[index] : (String(column), String(column), 1.0)
    return plot(
        t_hr,
        df[!, column] .* scale;
        label=false,
        xlabel="Time since start (hr)",
        ylabel=ylabel,
        title=title,
        lw=2
    )
end

"""
    save_comparison_plots(csv_path, fields) -> Vector{String}

One figure per saved field, written to `<csv dir>/plots`. Every component of
every field is plotted, so adding a name to `SAVED_VALUES` is all it takes to
see it.
"""
function save_comparison_plots(csv_path, fields)
    csv_path === nothing && return String[]
    isfile(csv_path) || return String[]
    df = CSV.read(csv_path, DataFrame)
    nrow(df) == 0 && return String[]

    plots_dir = joinpath(dirname(csv_path), "plots")
    mkpath(plots_dir)
    t = Float64.(df[!, :time])
    t_hr = (t .- first(t)) ./ 3600.0

    saved_paths = String[]
    for field in fields
        field isa Symbol || continue
        columns = _field_columns(df, field)
        isempty(columns) && continue
        panels = [_panel(df, t_hr, field, column, k) for (k, column) in enumerate(columns)]
        rows = cld(length(panels), 2)
        figure = plot(
            panels...;
            layout=(rows, min(2, length(panels))),
            size=(1100, 300 * rows),
            left_margin=8 * Plots.mm,
            bottom_margin=5 * Plots.mm
        )
        path = joinpath(plots_dir, "$(field).png")
        savefig(figure, path)
        push!(saved_paths, path)
    end

    if !isempty(saved_paths)
        println("Saved plots:")
        foreach(path -> println("  ", abspath(path)), saved_paths)
    end
    return saved_paths
end

# ===========================================================================
# Run.
# ===========================================================================

planet = make_no_gram_planet(:earth)
ephemerides_model = SimpleEphemeridesModel()
earth_harmonics_file = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")

# Circular speed at the start radius, raised for the requested eccentricity: a
# horizontal pass is an apsis, and this speed puts periapsis here.
start_radius_m = norm(
    CartesianInitialCondition(
        planet;
        lat=INITIAL_LATITUDE_DEG,
        lon=INITIAL_LONGITUDE_DEG,
        alt=INITIAL_ALTITUDE_M,
        speed=0.0,
        azimuth=0.0,
        initial_time=INITIAL_EPOCH,
        ephemerides_model=ephemerides_model
    ).pos
)
initial_speed_mps = sqrt(planet.μ / start_radius_m * (1.0 + ECCENTRICITY))

ic = CartesianInitialCondition(
    planet;
    lat=INITIAL_LATITUDE_DEG,
    lon=INITIAL_LONGITUDE_DEG,
    alt=INITIAL_ALTITUDE_M,
    speed=initial_speed_mps,
    flight_path_angle=FLIGHT_PATH_ANGLE_DEG,
    inclination=INCLINATION_DEG,
    descending=DESCENDING,
    initial_time=INITIAL_EPOCH,
    ephemerides_model=ephemerides_model
)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.05, 2.05, 2.8),
    panel_dims=(0.01, 5.7 / 2.0, 1.0),
    bus_mass=620.0,
    panel_mass_each=10.0,
    panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
    ic=ic,
    prop_mass=200.0,
    id=1
)

args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=MISSION_TIME_S,
    initial_time=INITIAL_EPOCH,
    dynamic_effectors=(GravitationalHarmonicsModel(HARMONICS_DEGREE, HARMONICS_ORDER, earth_harmonics_file, planet),),
    density_model=NoAtmosphereModel(),
    ephemerides_model=ephemerides_model,
    orientation_sim=false,
    keplerian=true,
    data_rate=DATA_RATE_S,
    verbose=true,
    results_directory=OUTPUT_DIR
)

# The external propagator has to start from the same inertial state, so print
# what the geodetic point resolved to at this epoch.
println("Comparison case")
println("  epoch                 : ", INITIAL_EPOCH)
println("  geodetic start        : lat ", INITIAL_LATITUDE_DEG, " deg, lon ", INITIAL_LONGITUDE_DEG,
        " deg, alt ", INITIAL_ALTITUDE_M, " m")
println("  inclination           : ", INCLINATION_DEG, " deg (", DESCENDING ? "descending" : "ascending", ")")
println("  gravity field         : GGM05C ", HARMONICS_DEGREE, "x", HARMONICS_ORDER)
println("  inertial position (m) : ", ic.pos)
println("  inertial velocity (m/s): ", ic.vel)

csv_path = run_and_report(args; save_fields=default_save_fields(args; extra=SAVED_VALUES))
save_comparison_plots(csv_path, SAVED_VALUES)
