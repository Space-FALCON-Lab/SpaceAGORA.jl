# Regression guard for `@simd ivdep` on the flat harmonics batch kernel
# (`_harmonics_flat_batch_kernel!`, src/dynamics/coupled/perturbations.jl).
#
# `@simd ivdep` licenses LLVM to vectorise the three per-degree batch loops
# unconditionally instead of only above the batch size LLVM's own alias
# analysis clears on its own -- see the kernel's docstring and
# `git log -1 cd212833aa`. It does not license reassociation: each loop
# iteration writes only its own satellite's slot, so there is no reduction for
# the annotation to reassociate and the kernel's floating-point operations stay
# in the scalar kernel's order regardless of vectorisation. This file pins that
# with two complementary checks `test/unit/simulation/harmonics_batch_parity_tests.jl`
# does not already cover: batch sizes straddling common SIMD lane widths
# (2/4/8 doubles), and a full end-to-end trajectory compared between the flat
# (batch/SIMD) and per-satellite (scalar) RHS routes rather than just the raw
# kernel call.
module HarmonicsKernelIdentityTests

using Test
using StaticArrays
using LinearAlgebra
using Random
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SM_HKI = SpaceAGORA.SimulationModel
const SE_HKI = SpaceAGORA.SimulationEngine
const PERT_HKI = SM_HKI.DynamicEffectors.PerturbationEffectors

const HKI_GRAVITY_FILE = joinpath(
    normpath(joinpath(@__DIR__, "..", "..", "..")),
    "data", "Gravity_harmonics_data", "EarthGGM05C.csv",
)

const HKI_LPI = SMatrix{3, 3, Float64, 9}(
    0.9362934, -0.3513767, 0.0,
    0.3513767,  0.9362934, 0.0,
    0.0,        0.0,       1.0,
)

function hki_states(n::Int, seed::Int)
    rng = Random.MersenneTwister(seed)
    states = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        r = 6.771e6 + 2.0e6 * rand(rng)
        u = 2 * rand(rng) - 1
        lon = 2pi * rand(rng)
        s = sqrt(1 - u^2)
        v = zeros(Float64, 8)
        v[1] = r * s * cos(lon)
        v[2] = r * s * sin(lon)
        v[3] = r * u
        v[5] = 7.5e3
        v[7] = 500.0 + 100 * rand(rng)   # mass
        states[i] = v
    end
    return states
end

function hki_batch_forces(model, states::Vector{Vector{Float64}})
    B = length(states)
    work_items = collect(1:B)
    slots = zeros(Float64, 6, 1, B)
    pool = PERT_HKI._get_harmonics_batch_pool(model, 1, B)
    PERT_HKI._harmonics_flat_batch_kernel!(
        slots, 1, model, states, work_items, 1, B, HKI_LPI, pool[1],
    )
    return slots
end

function hki_scalar_force(model, state::Vector{Float64})
    workspace = PERT_HKI._make_harmonics_scratch_workspace(model)
    pos = SVector{3, Float64}(state[1], state[2], state[3])
    force, _ = PERT_HKI._harmonics_scalar_force_ii(model, workspace, pos, state[7], HKI_LPI)
    return force
end

@testset "the SIMD batch kernel is exact at batch sizes straddling common vector widths" begin
    @test isfile(HKI_GRAVITY_FILE)
    planet = SM_HKI.Earth()
    model = SM_HKI.GravitationalHarmonicsModel(30, 15, HKI_GRAVITY_FILE, planet)
    # 2, 4 and 8 doubles is every lane width AVX2/AVX-512 double-precision
    # vectorisation could pick; one below, at, and one above each is where a
    # SIMD remainder-loop bug would show up first.
    boundary_sizes = sort!(unique(vcat(collect(1:20), [31, 32, 33, 63, 64, 65, 127, 128, 129])))
    for B in boundary_sizes
        states = hki_states(B, 1000 + B)
        slots = hki_batch_forces(model, states)
        for i in eachindex(states)
            expected = hki_scalar_force(model, states[i])
            for k in 1:3
                @test slots[k, 1, i] === expected[k]
            end
        end
    end
end

# A second, independently-shaped field (order-truncated rather than square)
# at the same boundary sizes, since the recursion's column bound depends on
# both L and M.
@testset "the SIMD batch kernel is exact at boundary sizes on an order-truncated field" begin
    planet = SM_HKI.Earth()
    model = SM_HKI.GravitationalHarmonicsModel(50, 6, HKI_GRAVITY_FILE, planet)
    for B in (2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 63, 64, 65)
        states = hki_states(B, 2000 + B)
        slots = hki_batch_forces(model, states)
        for i in eachindex(states)
            expected = hki_scalar_force(model, states[i])
            for k in 1:3
                @test slots[k, 1, i] === expected[k]
            end
        end
    end
end

# ── Simulation-level identity: flat (SIMD batch) route vs. per-satellite
# (scalar) route, full trajectory ──────────────────────────────────────────

function hki_spacecraft(planet, id::Int)
    root = Link(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(
        ra=planet.Rp_e + 550e3 + 733.0 * id,
        rp=planet.Rp_e + 500e3 + 611.0 * id,
        i=51.6, ω=17.0, Ω=(360.0 * id) / 37, ν=(360.0 * id) / 41,
    )
    return SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, id)
end

function hki_configuration(planet, n_sats::Int, L::Int, M::Int)
    harmonics = GravitationalHarmonicsModel(L, M, HKI_GRAVITY_FILE, planet)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=90.0, orientation_sim=false, num_steps_to_save=12, data_rate=8.0
        ),
        environment_model=EnvironmentModel(
            planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false,
            ephemerides_model=SimpleEphemeridesModel(),
        ),
        dynamics_model=DynamicsModel([hki_spacecraft(planet, i) for i in 1:n_sats], (harmonics,)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=20.0,
        ),
        solver_config=SolverConfig(solver_mode=:tsit5),
    )
end

hki_bits(v) = reinterpret.(UInt64, Float64.(collect(v)))

function hki_run(args, mode::String)
    fields = [
        SaveField(:pos, (u, t, int) -> [SVector{3, Float64}(sc.pos) for sc in u.sc]; per_satellite=true),
        SaveField(:vel, (u, t, int) -> [SVector{3, Float64}(sc.vel) for sc in u.sc]; per_satellite=true),
    ]
    recorder = TrajectoryRecorder(args; save_fields=fields)
    result = withenv(
        "SPACEAGORA_RHS_EXECUTION_MODE" => mode,
        "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "0",
        "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "0",
    ) do
        run_simulation(
            args; return_solver_metadata=true, return_solution=true, visualization=false,
            extra_callbacks=(get_trajectory_recorder_callback(recorder),),
        )
    end
    @test result.retcode == "Success"
    return (times=collect(trajectory_times(recorder)), rows=trajectory_save_data(recorder))
end

@testset "a harmonics-only constellation's flat and per-satellite routes agree bit for bit" begin
    planet = SM_HKI.Earth()
    # 33 clears the default flat-route admission floor and straddles the 32
    # doubles/4-double-lane boundary the kernel-level tests above target.
    for (n_sats, L, M) in ((33, 30, 15), (8, 50, 6))
        args = hki_configuration(planet, n_sats, L, M)
        flat = hki_run(args, "flat")
        satellite = hki_run(args, "satellite")
        @test hki_bits(flat.times) == hki_bits(satellite.times)
        @test length(flat.rows) == length(satellite.rows)
        @test all(zip(satellite.rows, flat.rows)) do (a, b)
            all(key -> all(i -> hki_bits(a[key][i]) == hki_bits(b[key][i]), 1:n_sats), keys(a))
        end
        # The comparison is meaningful only if the constellation actually moved.
        @test any(row -> norm(row[:pos][1]) > 0.0, flat.rows)
    end
end

end # module HarmonicsKernelIdentityTests
