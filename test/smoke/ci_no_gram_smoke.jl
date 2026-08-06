using Test
using StaticArrays
using SPICE
using SpaceAGORA

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SimulationEngine = SpaceAGORA.SimulationEngine
const SimulationModel = SpaceAGORA.SimulationModel
const RuntimeServices = SpaceAGORA.RuntimeServices
using .SimulationModel
const run_simulation = SpaceAGORA.run_simulation

import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

lock(RuntimeServices.SPICE_LOCK) do
    kclear()
end
# kclear() wipes CSPICE's kernel pool but not the furnished-kernel /
# planet-instance caches in src/environment/ephemerides/planets.jl; without
# this, a subsequent Earth(...)/Mars(...) call with a previously-seen key
# (e.g. from ci_clean_depot_smoke.jl / ci_threaded_smoke.jl running earlier
# in the same process via test/smoke/runtests.jl) silently skips
# re-furnishing and returns a planet built from kernels kclear() already
# wiped from CSPICE's pool.
SimulationModel.Planets._reset_furnished_kernels!()

spacecraft = make_three_body_spacecraft(
    bus_dims=(1.2, 1.1, 0.9),
    panel_dims=(0.01, 0.8, 0.4),
    bus_mass=150.0,
    panel_mass_each=2.0,
    panel_offset_y=0.7,
    ic=InitialCondition(
        ra=Mars().Rp_e + 220e3,
        rp=Mars().Rp_e + 150e3,
        i=28.0,
        ω=10.0,
        Ω=15.0,
        ν=165.0
    ),
    prop_mass=15.0,
    id=1
)

args = make_example_config(
    planet=Mars(),
    spacecraft=spacecraft,
    mission_time=90.0,
    initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredGravityModel(),),
    density_model=ExponentialAtmosphereModel(Mars()),
    ephemerides_model=SimpleEphemeridesModel(),
    orientation_sim=false,
    keplerian=true,
    EI_km=140.0,
    verbose=false,
    results=false,
    results_directory=joinpath(REPO_ROOT, "output", "ci_no_gram_smoke")
)

# The quickstart config is deliberately drag-free (exponential atmosphere for
# heating only, no aero effector), so the engine's density-without-aero
# diagnostic is expected here; pin it instead of asserting silence.
@test_logs (:warn, r"NO aerodynamic force will ever be applied") match_mode=:any run_simulation(args)
et0 = ephemerides_time_seconds(args.initial_time, args.environment_model.ephemerides_model)
l_pi0 = planet_frame_lpi(args.environment_model.planet, et0, args.environment_model.ephemerides_model)
@test !all(iszero, l_pi0)

println("No-GRAM smoke passed.")
