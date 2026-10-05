module SimulationResultsTests

using Test
using CSV
using DataFrames
using StaticArrays
using SpaceAGORA
using SpaceAGORA.SimulationModel
using SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

# Minimal no-GRAM, no-SPICE run: J2 gravity, no atmosphere, one hour.
function _results_test_config(; results::Bool, dir::String)
    planet = make_no_gram_planet(:earth)
    sc = make_three_body_spacecraft(
        bus_dims=(2.05, 2.05, 2.8), panel_dims=(0.01, 2.85, 1.0), bus_mass=620.0,
        panel_mass_each=10.0, panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
        ic=InitialCondition(ra=7.0e6, rp=7.0e6, i=45.0, ω=0.0, Ω=0.0, ν=0.0),
        prop_mass=200.0, id=1)
    return make_example_config(
        planet=planet, spacecraft=sc, mission_time=3600.0,
        initial_time=InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0),
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        density_model=NoAtmosphereModel(), ephemerides_model=SimpleEphemeridesModel(),
        orientation_sim=false, keplerian=true, EI_km=120.0, verbose=false,
        results=results, results_directory=dir)
end

# Same run with checkpointing on (1200 s segments over the 3600 s mission:
# three segments), via SimulationSettings fields.
function _checkpointed_config(; results::Bool, dir::String)
    a = _results_test_config(results=results, dir=dir)
    s = a.simulation_settings
    names = fieldnames(typeof(s))
    vals = NamedTuple{names}(map(n -> getfield(s, n), names))
    s2 = SpaceAGORA.SimulationModel.SimulationSettings(;
        merge(vals, (checkpoint_enabled=true, checkpoint_interval_s=1200.0))...)
    return SpaceAGORA.SimulationModel.SimConfig._with_configuration(a; simulation_settings=s2)
end

@testset "SimulationResults" begin
    csv_dir = mktempdir()
    run_simulation(_results_test_config(results=true, dir=csv_dir))
    csv_path = joinpath(csv_dir, "simulation_results.csv")
    @test isfile(csv_path)
    csv = CSV.read(csv_path, DataFrame)

    mem_dir = mktempdir()
    r = run_simulation(_results_test_config(results=false, dir=mem_dir); return_results=true)
    @test r isa SimulationResults
    @test isempty(readdir(mem_dir))                  # nothing written
    @test isempty(r.files)
    @test r.solution === nothing
    @test nrow(r.table) > 0
    @test names(r.table) == names(csv)
    @test nrow(r.table) == nrow(csv)
    for c in names(csv)
        eltype(csv[!, c]) <: Real || continue
        @test all(isapprox.(Float64.(r.table[!, c]), Float64.(csv[!, c]); rtol=1e-12, atol=1e-12))
    end

    # With results on, the returned table is the written one and the files are listed.
    csv_dir2 = mktempdir()
    r2 = run_simulation(_results_test_config(results=true, dir=csv_dir2); return_results=true)
    @test joinpath(csv_dir2, "simulation_results.csv") in r2.files
    @test names(r2.table) == names(csv)

    # Stale bundle from an earlier run in the same directory is not listed.
    r_stale = run_simulation(_results_test_config(results=false, dir=csv_dir); return_results=true)
    @test isfile(joinpath(csv_dir, "simulation_results.feather"))
    @test isempty(r_stale.files)

    # Checkpointed solve with results=false fills the table from saved_values.
    ck_a = run_simulation(_checkpointed_config(results=true, dir=mktempdir()); return_results=true)
    ck_b = run_simulation(_checkpointed_config(results=false, dir=mktempdir()); return_results=true)
    @test nrow(ck_b.table) > 0
    @test nrow(ck_b.table) == nrow(ck_a.table)
    @test names(ck_b.table) == names(ck_a.table)

    # Default return value is unchanged.
    @test run_simulation(_results_test_config(results=false, dir=mktempdir())) === nothing
    sol = run_simulation(_results_test_config(results=false, dir=mktempdir()); return_solution=true)
    @test length(sol.t) > 1

    # return_solution + return_results puts the solution on the struct.
    r3 = run_simulation(_results_test_config(results=false, dir=mktempdir());
                        return_results=true, return_solution=true)
    @test r3.solution !== nothing
    @test length(r3.solution.t) > 1
    @test_throws ArgumentError run_simulation(_results_test_config(results=false, dir=mktempdir());
                                              return_results=true, return_solver_metadata=true)
end

# RPO command log reachable through the configuration that ran.
module RPOResultsExample
using SpaceAGORA
using SpaceAGORA.SimulationModel
end
const _RPO_SPICE = joinpath(@__DIR__, "..", "..", "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
if !isdir(_RPO_SPICE)
    @info "SPICE kernels absent; skipping the RPO command-log results test" _RPO_SPICE
else
    # The RPO example needs the HYPR package; the unit runner has loaded it already.
    isdefined(Main, :SpaceAGORAHYPR) || Base.include(Main, joinpath(@__DIR__, "..", "helpers", "load_hypr.jl"))
    Base.include(RPOResultsExample, normpath(joinpath(@__DIR__, "..", "..", "examples", "Earth_RPO_CubeSat_MPC.jl")))
    @testset "SimulationResults RPO command log" begin
        dt = Base.invokelatest() do
            RPOResultsExample.build_rpo_cubesat_mpc_demo(;
            pso_n_particles=8, pso_n_iters=2, n_station_points=256, verbose=false,
            results_directory=mktempdir(), mission_time=10.0, data_rate_s=1.0,
            start_rtn=SVector(-8.0, -4.0, 2.0), goal_rtn=SVector(-5.5, -2.5, 1.0),
            record_control_commands=true)
        end
        r = run_simulation(dt.args; return_results=true)   # default isolate_state=true
        ran = only(c for c in r.configuration.control_model.control_effectors
                   if c isa SpaceAGORA.SimulationModel.RPOMPCControlModel)
        @test ran !== dt.control
        @test ran.command_log !== nothing
        @test !isempty(ran.command_log.t_s)
        @test isempty(dt.control.command_log.t_s)           # caller's object untouched
    end
end

end # module SimulationResultsTests
