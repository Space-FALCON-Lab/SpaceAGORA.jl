#!/usr/bin/env julia
# Generates SpaceAGORA_Examples feather files matching the GMAT scenario matrix
# defined in data/telemetry/GMAT_Basilisk_STK_comp.py.
#
# Recreates the exact 24 scenarios (4 bodies × 3 gravity degrees × 2 third-body
# settings) and writes GMAT-compatible feather files to
# data/telemetry/SpaceAGORA_Examples/, overwriting any existing files.
#
# Usage:
#   julia --project=. scripts/generate_gmat_comparison_data.jl

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

import Pkg
if something(Base.active_project(), "") != joinpath(REPO_ROOT, "Project.toml")
    Pkg.activate(REPO_ROOT; io=devnull)
end

using SpaceAGORA
using Arrow
using CSV
using DataFrames
using StaticArrays
using LinearAlgebra

if !isdefined(@__MODULE__, :SM)
    const SM = SpaceAGORA.SimulationModel
end
if !isdefined(@__MODULE__, :TV)
    const TV = SpaceAGORA.TelemetryVerification
end
if !isdefined(@__MODULE__, :run_simulation)
    const run_simulation = SpaceAGORA.run_simulation
end
if !isdefined(@__MODULE__, :SPICE_PATH)
    const SPICE_PATH = joinpath(REPO_ROOT, "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE")
end

import SpaceAGORA.TelemetryVerification: make_three_body_spacecraft

using .SM

const OUTPUT_DIR = joinpath(REPO_ROOT, "data", "telemetry", "SpaceAGORA_Examples")
mkpath(OUTPUT_DIR)

# Harmonics coefficient files used in the GMAT comparison matrix
const HARMONICS = Dict{String, String}(
    "Earth" => joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "egm96.csv"),
    "Mars"  => joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "Mars50c.csv"),
    "Venus" => joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "MGNP180U.csv"),
    "Luna"  => joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "LP165P.csv"),
)

# Third-body perturbers per central body, matching GMAT Python script
const THIRD_BODIES = Dict{String, Tuple}(
    "Earth" => ("sun", "moon"),
    "Mars"  => ("sun",),
    "Luna"  => ("earth", "sun"),
    "Venus" => ("sun",),
)

const PRIMARY_NAMES = Dict{String, String}(
    "Earth" => "Earth",
    "Mars"  => "Mars",
    "Luna"  => "Moon",  # SpaceAGORA uses "Moon", not "Luna"
    "Venus" => "Venus",
)

# Initial Keplerian elements from the GMAT Python script — sma in km, angles in degrees
const INITIAL_CONDITIONS = Dict{String, NamedTuple}(
    "Earth" => (sma_km=7000.0, ecc=0.00001, inc=28.5, raan=45.0, aop=0.0, ta=0.0),
    "Mars"  => (sma_km=4000.0, ecc=0.002,   inc=45.0, raan=0.0,  aop=0.0, ta=90.0),
    "Luna"  => (sma_km=2000.0, ecc=0.002,   inc=30.0, raan=0.0,  aop=0.0, ta=0.0),
    "Venus" => (sma_km=7000.0, ecc=0.001,   inc=30.0, raan=0.0,  aop=0.0, ta=180.0),
)

# order=0 for J2 matches GMAT Order=0 (zonal-only, no tesseral terms)
const GRAVITY_CONFIGS = [
    (degree=0,  order=0,  tag="J0"),
    (degree=2,  order=0,  tag="J2"),
    (degree=50, order=50, tag="J50"),
]

# Epoch matches GMAT Python script: 01 Jan 2026 12:00:00.000
const INITIAL_TIME   = SM.InitialTime(year=2026, month=1, day=1, hour=12, minute=0, second=0.0)
const MISSION_TIME_S = 1_000_000.0  # matches GMAT propagation duration
const DATA_RATE_S    = 10.0         # save every 10 s — matches GMAT MaxStep and output cadence

# Lightweight spacecraft used in the GMAT comparison matrix
const SC_BUS_DIMS       = (0.205, 0.37, 0.08)
const SC_PANEL_DIMS     = (0.01, 0.0285, 0.0001)
const SC_BUS_MASS_KG    = 29.0
const SC_PANEL_MASS_KG  = 0.0
const SC_PANEL_OFFSET_M = 2.45
const SC_PROP_MASS_KG   = 0.0

function _make_planet(body_name::String)
    body_name == "Earth" && return SM.Earth("", SPICE_PATH)
    body_name == "Mars"  && return SM.Mars("", SPICE_PATH)
    body_name == "Luna"  && return SM.Moon("", SPICE_PATH)
    body_name == "Venus" && return SM.Venus("", SPICE_PATH)
    error("Unknown body: $body_name")
end

function _build_effectors(body_name::String, degree::Int, order::Int, tb::Bool, planet)
    eff = Any[SM.InverseSquaredGravityModel()]
    if degree > 0
        push!(eff, SM.GravitationalHarmonicsModel(degree, order, HARMONICS[body_name], planet))
    end
    if tb
        push!(eff, SM.NBodyGravityModel(
            body_names=THIRD_BODIES[body_name],
            primary_body_name=PRIMARY_NAMES[body_name],
            planet=planet
        ))
    end
    return Tuple(eff)
end

# Reformat standard SpaceAGORA feather output into the GMAT-compatible column set:
# ElapsedSecs (s), SMA (km), ECC, INC (deg), X/Y/Z (km), VX/VY/VZ (km/s)
function _format_output(df_raw::DataFrame, planet)::DataFrame
    t   = Float64.(df_raw[!, :time])
    x_m = Float64.(df_raw[!, :sc1_pos_1])
    y_m = Float64.(df_raw[!, :sc1_pos_2])
    z_m = Float64.(df_raw[!, :sc1_pos_3])
    vx  = Float64.(df_raw[!, :sc1_vel_1])
    vy  = Float64.(df_raw[!, :sc1_vel_2])
    vz  = Float64.(df_raw[!, :sc1_vel_3])

    n = length(t)
    sma_km  = Vector{Float64}(undef, n)
    ecc_v   = Vector{Float64}(undef, n)
    inc_deg = Vector{Float64}(undef, n)

    for i in 1:n
        r = SVector{3, Float64}(x_m[i], y_m[i], z_m[i])
        v = SVector{3, Float64}(vx[i], vy[i], vz[i])
        try
            oe = TV.rvtoorbitalelement(r, v, planet)
            sma_km[i]  = oe[1] * 1.0e-3
            ecc_v[i]   = oe[2]
            inc_deg[i] = rad2deg(oe[3])
        catch
            sma_km[i]  = NaN
            ecc_v[i]   = NaN
            inc_deg[i] = NaN
        end
    end

    return DataFrame(
        ElapsedSecs = t,
        SMA         = sma_km,
        ECC         = ecc_v,
        INC         = inc_deg,
        X           = x_m .* 1.0e-3,
        Y           = y_m .* 1.0e-3,
        Z           = z_m .* 1.0e-3,
        VX          = vx .* 1.0e-3,
        VY          = vy .* 1.0e-3,
        VZ          = vz .* 1.0e-3,
    )
end

for body_name in ("Earth", "Mars", "Luna", "Venus")
    for grav in GRAVITY_CONFIGS
        for tb in (false, true)
            tb_tag   = tb ? "TBTrue" : "TBFalse"
            sim_id   = "Sim_$(body_name)_$(grav.tag)_$(tb_tag)"
            dst_feather = joinpath(OUTPUT_DIR, "$(sim_id).feather")
            dst_csv     = joinpath(OUTPUT_DIR, "$(sim_id).csv")

            println("Running $(sim_id)...")

            planet = _make_planet(body_name)
            icp    = INITIAL_CONDITIONS[body_name]
            sma_m  = icp.sma_km * 1e3
            ic     = SM.InitialCondition(
                ra = sma_m * (1.0 + icp.ecc),
                rp = sma_m * (1.0 - icp.ecc),
                i  = icp.inc,
                ω  = icp.aop,
                Ω  = icp.raan,
                ν  = icp.ta,
            )

            spacecraft = make_three_body_spacecraft(
                bus_dims=SC_BUS_DIMS,
                panel_dims=SC_PANEL_DIMS,
                bus_mass=SC_BUS_MASS_KG,
                panel_mass_each=SC_PANEL_MASS_KG,
                panel_offset_y=SC_PANEL_OFFSET_M,
                ic=ic,
                prop_mass=SC_PROP_MASS_KG,
                id=1,
            )

            eff = _build_effectors(body_name, grav.degree, grav.order, tb, planet)

            mktempdir() do tmp
                config = SM.SimulationConfiguration(
                    simulation_settings=SM.SimulationSettings(
                        results=true,
                        verbose=false,
                        results_directory=tmp,
                        generate_plots=false,
                        generate_filenames=false,
                        normalize=false,
                        save_csv=false,
                    ),
                    mission_configuration=SM.MissionConfiguration(
                        mission_type=SM.MissionTime,
                        keplerian=false,
                        number_of_orbits=1,
                        mission_time=MISSION_TIME_S,
                        orientation_sim=false,
                        num_steps_to_save=5000,
                        data_rate=DATA_RATE_S,
                    ),
                    environment_model=SM.EnvironmentModel(
                        planet=planet,
                        EI=300.0,
                        density_model=SM.NoAtmosphereModel(),
                        ephemerides_model=SM.SpiceEphemeridesModel(),
                        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
                        topography=false,
                        wind=false,
                    ),
                    dynamics_model=SM.DynamicsModel([spacecraft], eff),
                    guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
                    navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
                    control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
                    initial_time=INITIAL_TIME,
                    integration_tolerances=SM.IntegrationTolerances(
                        reltol_orbit=1e-12,
                        abstol_orbit=1e-12,
                        dt_max_orbit=10.0,
                    ),
                )

                t_elapsed = @elapsed run_simulation(config)

                src = joinpath(tmp, "simulation_results.feather")
                if isfile(src)
                    df_raw = DataFrame(Arrow.Table(src))
                    df_out = _format_output(df_raw, planet)
                    Arrow.write(dst_feather, df_out)
                    CSV.write(dst_csv, df_out)
                    println("  Saved $(sim_id).{feather,csv}  rows=$(nrow(df_out))  $(round(t_elapsed; digits=1))s")
                else
                    @warn "No output file generated for $(sim_id)"
                end
            end
        end
    end
end

println("\nAll simulations complete. Output: $(OUTPUT_DIR)")
