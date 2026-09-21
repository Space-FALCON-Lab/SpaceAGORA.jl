using OrdinaryDiffEq, DiffEqCallbacks
using LinearAlgebra, StaticArrays
using Printf, Statistics
using DelimitedFiles
using Arrow, Test

#############
# Constants #
#############
const MU       = 3.986004418e14
const C        = 3.0e8
const R_EARTH  = 6_378_137.0
const R_ATMDEF = R_EARTH + 100_000.0
const Ẑ        = SVector(0.0, 0.0, 1.0)

@inline idx(satellite, offset) = 6*(satellite-1) + offset

#############
# Functions #
#############
include("functions/1_LOS_Metrics.jl")
include("functions/2_Laser_Forces_ver2.jl")
include("functions/3_Dynamics.jl")
include("functions/4_Diagnostics.jl")
include("functions/5_OE_Converters.jl")
include("functions/6_OE_and_dv_in_RTN.jl")
include("functions/9_Runners.jl")
include("functions/12_CSV_Write_Read.jl")
include("functions/13_Feather_Write.jl")
include("extract_feather_encounters.jl")

############
# Settings #
############
Base.@kwdef struct Test16FeatherOptions
    helpers::Int                    = 10
    helper_altitude_km::Float64      = 1050.0
    target_altitude_km::Float64      = 1000.0
    target_inclination_deg::Float64  = 0.0
    helper_inclination_deg::Float64  = 0.0
    target_nu_deg::Float64           = 0.0
    target_ecc::Float64              = 0.0
    laser_range_km::Float64          = 200.0
    laser_power_w::Float64          = 10_000.0
    magnification::Float64          = 100.0
    mass_kg::Float64                = 227.0
end

########
# Main #
########
function test16_feather(; smoke=false)
    opts = Test16FeatherOptions()

    ### 1. Define orbits ###
    helper_num = opts.helpers

    helper_oe = [
        (a_m   = R_EARTH + opts.helper_altitude_km*1e3,
         e     = 0.0,
         i_deg = opts.helper_inclination_deg,
         Ω_deg = 0.0,
         ω_deg = 0.0,
         ν_deg = (360.0 / helper_num) * (helper - 1)) for helper in 1:helper_num
    ]
    target_orbit = (
        a_m   = R_EARTH + opts.target_altitude_km*1e3,
        e     = opts.target_ecc,
        i_deg = opts.target_inclination_deg,
        Ω_deg = 0.0,
        ω_deg = 0.0,
        ν_deg = opts.target_nu_deg
    )
    orbits = vcat(helper_oe, [target_orbit])

    ### 2. Laser / cavity parameters ###
    power_matrix = zeros(length(orbits), length(orbits))
    cavity = Dict{Tuple{Int,Int},Dict{Symbol,Any}}(
        (helper, helper_num + 1) => Dict(
            :B   => opts.magnification,
            :Pin => opts.laser_power_w
        ) for helper in 1:helper_num
    )

    ### 3. Simulation duration and output directory ###
    duration = smoke ? 60.0 : 63071.0
    paths = scenario_paths(joinpath(@__DIR__, "output"), opts, duration;
        source="prototype", smoke, schedule=:none)

    ### 4. Run the simulation without plots ###
    println("START SIMULATION (test16-feather): $(duration) s")
    sol, params, _, _ = mktempdir() do image_dir
        run_open_cavity_multi(orbits;
            mass_kg      = opts.mass_kg,
            Pm           = power_matrix,
            cavity       = cavity,
            use_los      = true,
            min_range    = 0.0,
            max_range    = opts.laser_range_km * 1e3,
            stop_on_dv   = false,
            T_seconds    = duration,
            verbose      = false,
            result_plots = false,
            target_only  = true,
            IMG_DIR      = image_dir,
            helper_num   = helper_num,
            gve_schedule = :none,
            record_events = true
        )
    end

    ### 5. Save results to CSV and Feather ###
    csv_path = save_timeseries_csv(sol, params, helper_oe, target_orbit; csv_dir=paths.csv)
    feather_path = save_timeseries_feather(sol, params; feather_dir=paths.feather)

    ### 6. Read back and verify the saved results ###
    table = Arrow.Table(feather_path)
    satellite_order = vcat(params[:target_ids], params[:helper_ids])

    @testset "test16-feather native schema and round-trip" begin
        ### 6.1. Solver completion and saved times ###
        @test string(sol.retcode) == "Success"
        @test sol.t[end] == duration
        @test table.time == sol.t
        @test table.time == output_times(duration)
        @test isfile(joinpath(paths.feather, "geometry_encounters.feather"))
        @test isfile(joinpath(paths.feather, "laser_on.feather"))

        ### 6.2. SpaceAGORA column names and order ###
        vector_fields = ("pos", "vel", "wind", "drag", "lift", "cross")
        expected_names = [:time]
        for field in ("pos", "vel", "altitude", "latitude_deg", "longitude_deg", "mass",
                      "wind", "drag", "lift", "cross", "periapsis_altitude", "heat_rate", "heat_load")
            suffixes = field in vector_fields ? ("_1", "_2", "_3") : ("",)
            for output_id in 1:params[:N], suffix in suffixes
                push!(expected_names, Symbol("sc$(output_id)_$(field)$(suffix)"))
            end
        end
        append!(expected_names, [
            :dv_r_accumulated,
            :dv_t_accumulated,
            :dv_n_accumulated,
            :laser_active_helper
        ])
        @test collect(propertynames(table)) == expected_names
        @test length(expected_names) == 25*params[:N] + 5

        ### 6.3. Target-first spacecraft mapping, masses, and orbital elements ###
        for (output_id, prototype_id) in enumerate(satellite_order)
            for (field, offset) in (("pos", 0), ("vel", 3)), component in 1:3
                values = getproperty(table, Symbol("sc$(output_id)_$(field)_$(component)"))
                expected_values = [state[idx(prototype_id, offset + component)] for state in sol.u]
                @test values == expected_values
                @test eltype(values) == Float64
            end
            saved_masses = getproperty(table, Symbol("sc$(output_id)_mass"))
            saved_longitudes = getproperty(table, Symbol("sc$(output_id)_longitude_deg"))
            @test saved_masses == fill(params[:masses][prototype_id], length(sol.t))
            @test all(ismissing, saved_longitudes)
        end
        @test isapprox(table.sc1_periapsis_altitude[1], opts.target_altitude_km*1e3)
        @test isapprox(table.sc2_periapsis_altitude[1], opts.helper_altitude_km*1e3)

        ### 6.4. Accumulated laser delta-v and active helper IDs ###
        _, delta_v = delta_v_RTN_time_series(sol, params)
        for (component, axis) in enumerate(("r", "t", "n"))
            values = getproperty(table, Symbol("dv_$(axis)_accumulated"))
            @test values == delta_v[only(params[:target_ids])][component, :]
            @test first(values) == 0.0
            @test all(isfinite, values)
        end
        @test any(!iszero, table.dv_r_accumulated)
        @test table.laser_active_helper[1] == 2
        @test all(helper -> helper == 0 || helper in 2:params[:N], table.laser_active_helper)

        ### 6.5. Export metadata, CSV round-trip, and invalid target IDs ###
        @test Arrow.getmetadata(table)["prototype_satellite_order"] == join(satellite_order, ",")
        csv_sol, metadata, csv_params = load_timeseries_csv(csv_path)
        @test csv_sol.t == sol.t
        @test csv_sol.u == sol.u
        bad_params = merge(params, Dict(:target_ids=>Int[]))
        @test_throws ArgumentError save_timeseries_feather(sol, bad_params; feather_dir=paths.feather)
    end

    ### 7. Extract pairwise encounter statistics and states ###
    extract_feather_encounters(feather_path; maximum_range_m=opts.laser_range_km * 1e3, output_dir=paths.csv)

    println("Verified $(length(table.time)) rows and $(length(propertynames(table))) columns: $feather_path")
    return feather_path
end

######################
# Command-line entry #
######################
if abspath(PROGRAM_FILE) == @__FILE__
    all(argument -> argument == "--smoke", ARGS) || error("Usage: test16-feather.jl [--smoke]")
    test16_feather(; smoke="--smoke" in ARGS)
end