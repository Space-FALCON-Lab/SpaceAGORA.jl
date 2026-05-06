const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")

const N_GRID = 3
const HARMONICS_DEGREE = 50
const HARMONICS_ORDER = 50
const OUTPUT_DIR = joinpath(REPO_ROOT, "output", "perturbation_aerobraking_study")

Base.@kwdef struct BodyStudyConfig
    planet_name::Symbol
    gram_planet_name::String
    harmonics_file::String
    nbody_names::Tuple{Vararg{String}}
    h_corridor_min_m::Float64
    h_corridor_max_m::Float64
    h_apo_initial_m::Float64
    EI_km::Float64
    mission_time_s::Float64
    inc_range_deg::Tuple{Float64,Float64}
    RAAN_range_deg::Tuple{Float64,Float64}
    omega_range_deg::Tuple{Float64,Float64}
end

h_corridor_mid(cfg::BodyStudyConfig) = (cfg.h_corridor_min_m + cfg.h_corridor_max_m) / 2.0

const BODY_CONFIGS = Dict{Symbol,BodyStudyConfig}(
    :earth => BodyStudyConfig(
        planet_name=:earth,
        gram_planet_name="earth",
        harmonics_file=joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv"),
        nbody_names=("Sun",),
        h_corridor_min_m=100_000.0,
        h_corridor_max_m=130_000.0,
        h_apo_initial_m=10_000_000.0,
        EI_km=200.0,
        mission_time_s=3600.0 * 300.0,
        inc_range_deg=(20.0, 90.0),
        RAAN_range_deg=(0.0, 300.0),
        omega_range_deg=(0.0, 300.0),
    ),
    :mars => BodyStudyConfig(
        planet_name=:mars,
        gram_planet_name="mars",
        harmonics_file=joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "GMM3.csv"),
        nbody_names=("Sun",),
        h_corridor_min_m=100_000.0,
        h_corridor_max_m=130_000.0,
        h_apo_initial_m=30_000_000.0,
        EI_km=160.0,
        mission_time_s=3600.0 * 800.0,
        inc_range_deg=(10.0, 100.0),
        RAAN_range_deg=(0.0, 300.0),
        omega_range_deg=(0.0, 300.0),
    ),
    :venus => BodyStudyConfig(
        planet_name=:venus,
        gram_planet_name="venus",
        harmonics_file=joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "MGNP180U.csv"),
        nbody_names=("Sun",),
        h_corridor_min_m=130_000.0,
        h_corridor_max_m=165_000.0,
        h_apo_initial_m=60_000_000.0,
        EI_km=250.0,
        mission_time_s=3600.0 * 24.0 * 50.0,
        inc_range_deg=(10.0, 90.0),
        RAAN_range_deg=(0.0, 300.0),
        omega_range_deg=(0.0, 300.0),
    ),
    :titan => BodyStudyConfig(
        planet_name=:titan,
        gram_planet_name="titan",
        harmonics_file=joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "titan5.csv"),
        nbody_names=("Saturn",),
        h_corridor_min_m=700_000.0,
        h_corridor_max_m=900_000.0,
        h_apo_initial_m=2_000_000.0,
        EI_km=1200.0,
        mission_time_s=3600.0 * 200.0,
        inc_range_deg=(10.0, 90.0),
        RAAN_range_deg=(0.0, 300.0),
        omega_range_deg=(0.0, 300.0),
    ),
)

const PERTURBATION_LEVELS = [:point_mass, :harmonics, :srp, :nbody, :full]

_ic_range(lo::Float64, hi::Float64, n::Int) = n == 1 ? [lo] : range(lo, hi; length=n)

function build_ic_grid(cfg::BodyStudyConfig; n::Int=N_GRID)
    h_peri_m = h_corridor_mid(cfg)
    inc_vals   = _ic_range(cfg.inc_range_deg[1],  cfg.inc_range_deg[2],  n)
    RAAN_vals  = _ic_range(cfg.RAAN_range_deg[1], cfg.RAAN_range_deg[2], n)
    omega_vals = _ic_range(cfg.omega_range_deg[1], cfg.omega_range_deg[2], n)

    ics = NamedTuple[]
    for inc_deg in inc_vals, RAAN_deg in RAAN_vals, omega_deg in omega_vals
        push!(ics, (
            h_peri_m=h_peri_m,
            h_apo_m=cfg.h_apo_initial_m,
            inc_deg=Float64(inc_deg),
            RAAN_deg=Float64(RAAN_deg),
            omega_deg=Float64(omega_deg),
        ))
    end
    return ics
end
