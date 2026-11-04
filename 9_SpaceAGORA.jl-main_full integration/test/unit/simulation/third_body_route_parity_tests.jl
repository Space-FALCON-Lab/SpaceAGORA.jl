module ThirdBodyRouteParityTests

# Guards for the three engine-side changes measured in
# docs/architecture/third_body_cost.md.
#
# 1. The RHS route at a thread budget of one. A constellation whose effectors
#    are all served by the flat route's serial pre-passes (the batchable
#    NBody/SRP/inverse-square kernels and the harmonics SIMD pre-pass) used to
#    fall to `:satellite_batch` as soon as it had a second effector, because
#    the `budget <= 1` guard came first. That gives up the batched coefficient
#    sweep and the once-per-evaluation shared body samples for no threading
#    reason at all: at allotment 1 the flat route spawns no tasks. The first
#    testset pins the admission and its boundary — one non-pre-pass effector in
#    the stack still routes to the per-satellite batch.
#
# 2. The per-spacecraft third-body sample. `sample_third_body_ephemerides`
#    built its position tuple with `ntuple` over a runtime length, which Julia
#    cannot infer, so every spacecraft of every derivative evaluation allocated
#    the tuple and the boxed closure behind it. The replacement maps over the
#    model's body-name tuple, whose length is part of the type. The second
#    testset pins that the two constructions produce the same positions bit for
#    bit, in the same order, and that the mapped one allocates less.
#
# 3. The solver family. `:auto_stiff` refused its explicit fast path for any
#    configuration carrying solar radiation pressure, so a constellation whose
#    dynamics are not stiff paid AutoTsit5's switch to Rodas5P and the whole
#    Rosenbrock W path. SRP is smooth away from the shadow boundary, so it now
#    belongs to the smooth-gravity list; the third testset pins the admission,
#    the boundary (a genuinely non-smooth effector still refuses) and the
#    escape hatch.

using Test
using StaticArrays
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SM_TBR = SpaceAGORA.SimulationModel
const SE_TBR = SpaceAGORA.SimulationEngine
const PE_TBR = SM_TBR.DynamicEffectors.PerturbationEffectors

# A degree-2 field written to a temporary file: the harmonics model only has a
# file-backed constructor, and nothing here depends on the coefficients'
# magnitude, only on the effector's type.
function tbr_write_degree2_coefficients(path::String)
    open(path, "w") do io
        println(io, "degree,order,C,S")
        println(io, "0,0,1.0,0.0")
        println(io, "1,0,0.0,0.0")
        println(io, "1,1,0.0,0.0")
        println(io, "2,0,-1.08262668e-3,0.0")
        println(io, "2,1,0.0,0.0")
        println(io, "2,2,1.5745e-6,-9.0387e-7")
    end
    return path
end

function tbr_spacecraft(planet, id::Int)
    root = Link(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(
        ra=planet.Rp_e + 550e3 + 1e3 * id,
        rp=planet.Rp_e + 500e3 + 1e3 * id,
        i=53.0, ω=24.0, Ω=10.0, ν=(360.0 * id) / 32,
    )
    return SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, id)
end

function tbr_configuration(planet, effectors::Tuple, n_sats::Int)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=10.0, orientation_sim=false, num_steps_to_save=10, data_rate=1.0
        ),
        environment_model=EnvironmentModel(
            planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false,
            # No SPICE: routing and sample construction are what this file is
            # about, and a fresh checkout has no kernels.
            ephemerides_model=SimpleEphemeridesModel(),
        ),
        dynamics_model=DynamicsModel([tbr_spacecraft(planet, i) for i in 1:n_sats], effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=20.0
        ),
    )
end

# The route decision as a solve at one thread would take it: no outer split
# advertised, no forced execution mode, and an inner budget of one.
const TBR_BUDGET1_ENV = [
    "SPACEAGORA_INNER_THREAD_BUDGET" => "1",
    "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
    "SPACEAGORA_RHS_EXECUTION_MODE" => nothing,
    "SPACEAGORA_PARALLEL_PROFILE" => nothing,
    "SPACEAGORA_RHS_PLAN_STEP_CACHE" => nothing,
    "SPACEAGORA_HARMONICS_BATCH_ENABLED" => nothing,
    "SPACEAGORA_EFFECTOR_FLAT_MIN_SATS" => nothing,
]

function tbr_plan(effectors::Tuple, n_sats::Int)
    planet = Earth()
    args = tbr_configuration(planet, effectors, n_sats)
    p = SM_TBR.ODEParams(n_sats=n_sats, args=args)
    return withenv(TBR_BUDGET1_ENV...) do
        SE_TBR._rhs_execution_plan(args, p, effectors, n_sats)
    end
end

# The interpolated third-body cache a GRAM/SPICE-backed run reads inside the
# RHS, built here from an analytic curve so the test needs no kernels.
function tbr_nbody_cache(body_query_names::Vector{String})
    ets = collect(0.0:600.0:6000.0)
    positions = Matrix{SVector{3, Float64}}(undef, length(ets), length(body_query_names))
    for (i, et) in pairs(ets), k in eachindex(body_query_names)
        angle = 1e-6 * et + 0.7 * k
        radius = k == 1 ? 1.49e11 : 3.84e8
        positions[i, k] = SVector{3, Float64}(
            radius * cos(angle), radius * sin(angle), 0.1 * radius * sin(0.5 * angle)
        )
    end
    return SM_TBR.NBodyEphemerisCache(
        PE_TBR._spice_query_name("Earth"),
        body_query_names,
        Dict(name => k for (k, name) in pairs(body_query_names)),
        ets,
        positions,
    )
end

# The pre-change construction, kept verbatim as the reference the replacement
# has to reproduce: `ntuple` over a runtime length, same body order, same
# cache-then-backend fallback.
function tbr_reference_positions(model, p, t::Float64)
    et = p.shared_buffers.et_start[] + t
    primary_body_name = SM_TBR.DynamicEffectors._spice_query_name(model.primary_body_name)
    memo_enabled = p.shared_buffers.spice_rhs_memo_enabled[]
    memo = p.shared_buffers.spice_rhs_memo
    cache_entry = p.shared_buffers.nbody_ephemeris_cache[]
    return ntuple(length(model.body_names)) do k
        body_name_spice = SM_TBR.DynamicEffectors._spice_query_name(model.body_names[k])
        pos = if cache_entry isa SM_TBR.NBodyEphemerisCache
            cached = SM_TBR.DynamicEffectors._nbody_body_position_from_cache_j2000_m(
                cache_entry, et, body_name_spice, primary_body_name,
            )
            cached === nothing ?
                PE_TBR._nbody_body_position_from_spice_j2000_m(
                    body_name_spice, et, primary_body_name, memo_enabled, memo,
                    p.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls,
                ) :
                cached
        else
            PE_TBR._nbody_body_position_from_spice_j2000_m(
                body_name_spice, et, primary_body_name, memo_enabled, memo,
                p.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls,
            )
        end
        return SVector{3, Float64}(pos)
    end
end

tbr_bits(v::SVector{3, Float64}) = ntuple(i -> reinterpret(UInt64, v[i]), 3)

@testset "third-body route and sample construction" begin

    @testset "a pre-pass-only stack takes the flat route at a thread budget of one" begin
        planet = Earth()
        mktempdir() do dir
            harmonics = GravitationalHarmonicsModel(
                2, 2, tbr_write_degree2_coefficients(joinpath(dir, "degree2.csv")), planet
            )
            nbody = NBodyGravityModel(
                body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet
            )
            srp = SolarRadiationPressureModel(1.2, 12.0)
            n_sats = 32   # above the default flat-route minimum of 24

            # Harmonics + n-body + SRP: every effector is served by a pre-pass,
            # so the flat route never builds its per-(satellite, effector)
            # queue and never spawns a task at allotment 1.
            plan = tbr_plan((harmonics, srp, nbody), n_sats)
            @test plan.mode === :flat_constellation_effector_queue
            @test plan.allotment == 1
            @test plan.dominant_axis === :flat_effector
            @test plan.effector_decision.use_threads == false
            @test plan.effector_decision.allotment == 1

            # The two-effector shapes the guard used to send to the
            # per-satellite batch as well.
            for stack in ((harmonics, nbody), (harmonics, srp))
                two = tbr_plan(stack, n_sats)
                @test two.mode === :flat_constellation_effector_queue
                @test two.allotment == 1
            end

            # One effector outside the pre-passes and the stack goes back to
            # the per-satellite batch, because the flat queue itself would have
            # to run. The fM aerodynamic model has a pre-pass of its own, so the
            # effector used here is its per-link-atmosphere variant, which
            # samples the atmosphere per link and stays on the queue.
            mixed = tbr_plan((harmonics, srp, nbody, AerodynamicCoefficientfM(per_link_atmosphere=true)), n_sats)
            @test mixed.mode === :satellite_batch
            @test mixed.allotment == 1

            # A constellation too small for the flat route's minimum is not
            # admitted either: the admission is about the batched sweep paying
            # for the route, not about the thread budget.
            small = tbr_plan((harmonics, srp, nbody), 4)
            @test small.mode === :satellite_batch
        end
    end

    @testset "the mapped third-body sample reproduces the ntuple construction bit for bit" begin
        planet = Earth()
        nbody = NBodyGravityModel(
            body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet
        )
        n_sats = 4
        args = tbr_configuration(planet, (nbody,), n_sats)
        p = SM_TBR.ODEParams(n_sats=n_sats, args=args)
        p.shared_buffers.nbody_ephemeris_cache[] = tbr_nbody_cache(
            [PE_TBR._spice_query_name("Sun"), PE_TBR._spice_query_name("Moon")]
        )
        p.shared_buffers.et_start[] = 0.0
        u = SE_TBR.build_initial_conditions(args)

        for t in (0.0, 1234.5, 2750.0, 6000.0)
            reference = tbr_reference_positions(nbody, p, t)
            for sat_idx in 1:n_sats
                sample = SE_TBR.sample_third_body_ephemerides(nbody, u.sc[sat_idx], p, sat_idx, t)
                @test sample.names === nbody.body_names
                @test length(sample.positions_ii) == length(reference)
                for k in eachindex(reference)
                    # Bit patterns, not ≈: the route parity this sample feeds
                    # is a byte-for-byte claim.
                    @test tbr_bits(sample.positions_ii[k]) == tbr_bits(reference[k])
                end
            end
        end

        # The point of mapping the tuple rather than building it from a
        # runtime length: the sample costs strictly less per spacecraft per
        # evaluation than the `ntuple` reference above, which is the whole
        # measured effect (-7.8% of a whole solve's allocation on the 256- and
        # 1024-spacecraft n-body rungs, docs/architecture/third_body_cost.md).
        # This is a comparison, not an absolute floor: both constructions carry
        # the same surrounding per-call cost, and only the difference between
        # them belongs to this change.
        sample_many(n) = begin
            last = SE_TBR.sample_third_body_ephemerides(nbody, u.sc[1], p, 1, 2750.0)
            for _ in 2:n
                last = SE_TBR.sample_third_body_ephemerides(nbody, u.sc[1], p, 1, 2750.0)
            end
            last
        end
        reference_many(n) = begin
            last = tbr_reference_positions(nbody, p, 2750.0)
            for _ in 2:n
                last = tbr_reference_positions(nbody, p, 2750.0)
            end
            last
        end
        sample_many(2); reference_many(2)   # warm up compilation
        @test @allocated(sample_many(256)) < @allocated(reference_many(256))
    end

    @testset "solar radiation pressure no longer disqualifies the auto-stiff fast path" begin
        planet = Earth()
        mktempdir() do dir
            harmonics = GravitationalHarmonicsModel(
                2, 2, tbr_write_degree2_coefficients(joinpath(dir, "degree2.csv")), planet
            )
            nbody = NBodyGravityModel(
                body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet
            )
            srp = SolarRadiationPressureModel(1.2, 12.0)
            cfg = SM_TBR.SolverConfig(solver_mode=:auto_stiff)

            # The P6 constellation's stack: harmonics + SRP + third-body
            # gravity. Every term is a smooth function of position (SRP's only
            # non-smooth feature is the shadow-boundary kink), so `:auto_stiff`
            # resolves to plain Tsit5 instead of AutoTsit5(Rodas5P).
            smooth = tbr_configuration(planet, (harmonics, srp, nbody), 8)
            @test SE_TBR._auto_stiff_smooth_gravity_reject_reason(cfg, smooth) === nothing
            @test SE_TBR._auto_stiff_smooth_gravity_eligible(cfg, smooth)

            # SRP alone, and SRP with the solar requirement as the only reason
            # it used to be refused.
            @test SE_TBR._auto_stiff_smooth_gravity_eligible(
                cfg, tbr_configuration(planet, (harmonics, srp), 8))

            # A genuinely non-smooth effector still refuses the fast path:
            # aerodynamics samples the atmosphere, whose density gradient is
            # what the implicit solver exists for here.
            aero = tbr_configuration(planet, (harmonics, srp, AerodynamicCoefficientfM()), 8)
            reason = SE_TBR._auto_stiff_smooth_gravity_reject_reason(cfg, aero)
            @test reason !== nothing
            @test occursin("AerodynamicCoefficientfM", reason)
            @test !SE_TBR._auto_stiff_smooth_gravity_eligible(cfg, aero)

            # The escape hatch a configuration that really wants the implicit
            # solver still has.
            off = SM_TBR.SolverConfig(solver_mode=:auto_stiff, auto_stiff_gravity_tsit5=false)
            @test !SE_TBR._auto_stiff_smooth_gravity_eligible(off, smooth)
        end
    end
end

end # module
