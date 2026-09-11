# Configuration struct and validation for the ORACLE open-cavity laser-link case.
# Placed in core/types per project convention for option/config structs.

using Printf

export OracleOptions

const ORACLE_PAPER_TARGET_ALTITUDES_KM       = (1150.0, 1050.0, 1000.0, 950.0, 850.0)
const ORACLE_PAPER_TARGET_INCLINATIONS_DEG   = (0.0, 0.5, 1.0)
const ORACLE_PAPER_HELPER_COUNTS             = (1, 50, 100, 150, 200, 250, 300)
const ORACLE_PAPER_FIXED_HELPER_ALTITUDE_KM  = 1000.0
const ORACLE_PAPER_FIXED_HELPER_INCLINATION_DEG = 0.0

"""
    OracleOptions

All tunable parameters for one ORACLE open-cavity laser-link scenario.
Pass to `build_oracle_case_config` and `run_oracle_open_cavity_case`.
"""
Base.@kwdef struct OracleOptions
    helpers::Int                   = 10
    helper_altitude_km::Float64    = 1050.0
    target_altitude_km::Float64    = 1000.0
    target_inclination_deg::Float64 = 0.0
    helper_inclination_deg::Float64 = 0.0
    target_nu_deg::Float64         = 0.0
    target_ecc::Float64            = 0.0
    orbits::Float64                = 10.0
    schedule::Symbol               = :naive_next_entering
    laser_range_km::Float64        = 200.0
    laser_power_w::Float64         = 10_000.0
    magnification::Float64         = 100.0
    beta::Float64                  = 1.0
    eta::Float64                   = 1.0
    mass_kg::Float64               = 227.0
    dt_max_s::Float64              = 10.0
    planet::Symbol                 = :earth  # planet symbol passed to make_no_gram_planet
    paper_grid::Bool               = false
    feather_only::Bool             = false
    output_dir::String             = "output"
    timeseries_points::Int         = 1001
    animate::Bool                  = false
end

function _validate_options(opts::OracleOptions)
    opts.helpers >= 1 || throw(ArgumentError("helpers must be >= 1."))
    opts.helper_altitude_km > 0.0 || throw(ArgumentError("helper_altitude_km must be positive."))
    opts.target_altitude_km > 0.0 || throw(ArgumentError("target_altitude_km must be positive."))
    (0.0 <= opts.target_ecc < 1.0) || throw(ArgumentError("target_ecc must be in [0, 1)."))
    opts.orbits > 0.0 || throw(ArgumentError("orbits must be positive."))
    opts.schedule in (:naive_next_entering, :positive_along_track,
                      :gve_sma, :gve_ecc, :gve_inc, :gve_raan, :gve_argp) ||
        throw(ArgumentError("schedule must be one of the supported ORACLE schedules."))
    opts.laser_range_km >= 0.0 || throw(ArgumentError("laser_range_km must be nonnegative."))
    opts.laser_power_w >= 0.0 || throw(ArgumentError("laser_power_w must be nonnegative."))
    opts.magnification >= 0.0 || throw(ArgumentError("magnification must be nonnegative."))
    opts.beta >= 0.0 || throw(ArgumentError("beta must be nonnegative."))
    opts.eta >= 0.0 || throw(ArgumentError("eta must be nonnegative."))
    opts.mass_kg > 0.0 || throw(ArgumentError("mass_kg must be positive."))
    opts.dt_max_s > 0.0 || throw(ArgumentError("dt_max_s must be positive."))
    opts.timeseries_points >= 2 || throw(ArgumentError("timeseries_points must be >= 2."))
    opts.planet in (:earth, :mars, :venus, :titan, :moon) ||
        throw(ArgumentError("planet must be a symbol supported by make_no_gram_planet (e.g. :earth, :mars)."))
    return nothing
end

# Return a copy of opts with keyword overrides applied.
_with(opts::OracleOptions; kwargs...) =
    OracleOptions(; (f => getfield(opts, f) for f in fieldnames(OracleOptions))..., kwargs...)
