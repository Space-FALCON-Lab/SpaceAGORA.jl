module GramDensityPerturbationCallbackTests
# The opt-in GRAM perturbation callback without native GRAM: env parsing, the
# pass interpolant, the diagnostics log, and the accepted-step update in each
# mode. The walk instance is a stub type, so the per-pass reseed, the step
# sample and the pass sweep run against known draws. The callback reads only
# `p`, `u` and `t` off its integrator, so a NamedTuple stands in for one.
using Test, SpaceAGORA, SpaceAGORA.SimulationModel, StaticArrays

const SE = SpaceAGORA.SimulationEngine
const CB = SpaceAGORA.SimulationModel.SimulationCallbacks
const EM = SpaceAGORA.SimulationModel.EnvironmentModels

const MODE_ENV = "SPACEAGORA_GRAM_DENSITY_PERTURBATION"
const DT_ENV = "SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_DT_S"
const MAX_ENV = "SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_MAX_S"
const LOG_ENV = "SPACEAGORA_GRAM_DENSITY_PERTURBATION_LOG"
const RESEED_ENV = "SPACEAGORA_GRAM_DENSITY_PERTURBATION_PASS_RESEED"

# A walk whose k-th draw is r = 1 + 0.1 k, and which records its reseeds and
# which of its updates were marked as first after a clone or reseed.
mutable struct StubWalk
    seed::Int
    calls::Int
    first_updates::Vector{Bool}
end
StubWalk(seed::Int, calls::Int) = StubWalk(seed, calls, Bool[])
function EM._gram_walk_sample(w::StubWalk, h::Float64, lat::Float64, lon::Float64, t::Float64;
                              first_update::Bool=false)
    w.calls += 1
    push!(w.first_updates, first_update)
    return (1.0e-9 * (1.0 + 0.1 * w.calls), 1.0e-9, 0.05, 0.25)
end
EM._gram_walk_reseed!(w::StubWalk, seed::Int) = (w.seed = seed; nothing)

const N_SATS = 2

function perturbation_spacecraft(planet, id::Int)
    root = Link(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(ra=planet.Rp_e + 540e3, rp=planet.Rp_e + 500e3,
                          i=35.0, ω=40.0, Ω=10.0, ν=120.0 + 12.0 * (id - 1))
    return SpacecraftModel(Joint[], Link[root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, id)
end

# EI = 1000 km puts both ~500 km spacecraft inside the atmosphere; tests move
# them out by lowering the state's ei_m.
function perturbation_config(results_directory::String; resume::Bool=false)
    planet = Earth()
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false,
                                               normalize=false, save_csv=false,
                                               results_directory=results_directory,
                                               resume_from_checkpoint=resume),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, keplerian=true,
                                                   number_of_orbits=1, mission_time=600.0,
                                                   orientation_sim=false, num_steps_to_save=10,
                                                   data_rate=10.0),
        environment_model=EnvironmentModel(planet=planet, EI=1000.0,
                                           density_model=ExponentialAtmosphereModel(planet),
                                           ephemerides_model=SimpleEphemeridesModel(),
                                           thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
                                           topography=false, wind=false),
        dynamics_model=DynamicsModel([perturbation_spacecraft(planet, i) for i in 1:N_SATS],
                                     (InverseSquaredJ2GravityModel(), AerodynamicCoefficientfM())),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=20.0),
    )
end

function perturbation_params(args)
    p = SE.ODEParams(n_sats=N_SATS, args=args)
    SE._initialize_heat_rate_buffers!(p)
    SE._initialize_harmonics_workspace_buffers!(p)
    SE._initialize_save_cache_buffers!(p)
    p.shared_buffers.et_start[] =
        SimulationModel.ephemerides_time_seconds(args.initial_time, args.environment_model.ephemerides_model)
    return p
end

perturbation_callback(args, env::Pair...) =
    withenv(env...) do
        CB.get_gram_density_perturbation_callback(N_SATS, args)
    end

@testset "perturbation env parsing" begin
    for (raw, mode) in ((nothing, :off), ("", :off), (" OFF ", :off), ("0", :off), ("none", :off),
                        ("step", :step), ("A", :step), ("accepted_step", :step),
                        ("pass", :pass), ("b", :pass), ("lookahead", :pass),
                        ("naive_rhs", :naive_rhs), ("naive", :naive_rhs))
        withenv(MODE_ENV => raw) do
            @test CB._gram_density_perturbation_mode() === mode
        end
    end
    withenv(MODE_ENV => "bogus") do
        @test_throws ArgumentError CB._gram_density_perturbation_mode()
    end

    withenv(DT_ENV => nothing, MAX_ENV => nothing, RESEED_ENV => nothing) do
        @test CB._gram_density_perturbation_pass_dt_s() == 1.0
        @test CB._gram_density_perturbation_pass_max_s() == 7200.0
        @test CB._gram_density_perturbation_pass_reseed() == false
    end
    withenv(DT_ENV => "0.5", MAX_ENV => "60", RESEED_ENV => "1") do
        @test CB._gram_density_perturbation_pass_dt_s() == 0.5
        @test CB._gram_density_perturbation_pass_max_s() == 60.0
        @test CB._gram_density_perturbation_pass_reseed() == true
    end
    for bad in ("0", "-1", "Inf", "NaN", "abc")
        withenv(DT_ENV => bad) do
            @test_throws ArgumentError CB._gram_density_perturbation_pass_dt_s()
        end
        withenv(MAX_ENV => bad) do
            @test_throws ArgumentError CB._gram_density_perturbation_pass_max_s()
        end
    end
    withenv(RESEED_ENV => "maybe") do
        @test_throws ArgumentError CB._gram_density_perturbation_pass_reseed()
    end
end

@testset "per-pass seed and ratio" begin
    seeds = [CB._gram_pass_seed(base, pass) for base in (1, 1001, typemax(Int)) for pass in 1:50]
    @test all(s -> 1 <= s <= 16_777_215, seeds)
    @test CB._gram_pass_seed(1001, 3) == CB._gram_pass_seed(1001, 3)
    @test CB._gram_pass_seed(1001, 3) != CB._gram_pass_seed(1001, 4)
    @test CB._gram_pass_seed(1001, 3) != CB._gram_pass_seed(1002, 3)

    @test CB._gram_ratio(3.0, 2.0) == 1.5
    @test CB._gram_ratio(3.0, 0.0) == 1.0
    @test CB._gram_ratio(NaN, 2.0) == 1.0
    @test CB._gram_ratio(Inf, 2.0) == 1.0
end

@testset "pass interpolant and predicted altitude" begin
    st = CB._new_gram_density_perturbation_state(:pass, 1, 120e3, 2.0, 10.0, "")
    @test CB._gram_pass_factor(st, 1, 0.0) == 1.0              # not active
    st.pass_active[1] = true
    @test CB._gram_pass_factor(st, 1, 0.0) == 1.0              # no knots
    @test isnan(CB._gram_pass_predicted_alt(st, 1, 0.0))
    st.pass_t0[1] = 10.0
    st.pass_r[1] = [1.5]
    @test CB._gram_pass_factor(st, 1, 99.0) == 1.5             # one knot held
    st.pass_r[1] = [1.0, 2.0, 4.0]
    st.pass_alt[1] = [100e3, 90e3, 110e3]
    @test CB._gram_pass_factor(st, 1, 5.0) == 1.0              # before the first knot
    @test CB._gram_pass_factor(st, 1, 11.0) == 1.5             # halfway, knots 1-2
    @test CB._gram_pass_factor(st, 1, 13.0) == 3.0             # halfway, knots 2-3
    @test CB._gram_pass_factor(st, 1, 14.0) == 4.0             # last knot
    @test CB._gram_pass_factor(st, 1, 50.0) == 4.0             # held past the prediction
    @test CB._gram_pass_predicted_alt(st, 1, 5.0) == 100e3
    @test CB._gram_pass_predicted_alt(st, 1, 13.0) == 100e3
    @test isnan(CB._gram_pass_predicted_alt(st, 1, 14.0))
end

@testset "diagnostics log and summary" begin
    mktempdir() do dir
        st = CB._new_gram_density_perturbation_state(:step, 1, 120e3, 1.0, 10.0, "")
        @test CB._write_gram_perturbation_log(st) === nothing   # no path, no files
        @test isempty(readdir(dir))

        path = joinpath(dir, "nested", "log.csv")
        st = CB._new_gram_density_perturbation_state(:step, 2, 120e3, 1.0, 10.0, path, true)
        st.pass_count[2] = 3
        CB._gram_perturbation_log!(st, 1, 2, 4.5, 110e3, 1.25, 0.05, 2.0e-9, 0.5)
        CB._write_gram_perturbation_log(st)
        lines = readlines(path)
        @test lines[1] == "kind,sat,pass,t_s,alt_m,r,sigma_frac,mean_density_kgm3,aux"
        @test lines[2] == "1,2,3,4.5,110000.0,1.25,0.05,2.0e-9,0.5"
        @test length(lines) == 2
        summary = read(path * ".summary.toml", String)
        @test occursin("mode = \"step\"", summary)
        @test occursin("reseed = true", summary)
        @test occursin("pass_count = [0, 3]", summary)
        @test occursin("log_rows = 1", summary)
    end
end

@testset "callback construction and naive negative control" begin
    mktempdir() do dir
        args = perturbation_config(dir)
        u = SE.build_initial_conditions(args)
        @test perturbation_callback(args, MODE_ENV => "off") === nothing
        @test_throws ArgumentError perturbation_callback(args, MODE_ENV => "pass", DT_ENV => "-1")
        # A checkpoint does not carry the walk, so resuming with a mode on is refused.
        resumed = perturbation_config(dir; resume=true)
        for mode in ("step", "pass", "naive_rhs")
            @test_throws ArgumentError perturbation_callback(resumed, MODE_ENV => mode)
        end
        @test perturbation_callback(resumed, MODE_ENV => "off") === nothing

        cb = perturbation_callback(args, MODE_ENV => "naive_rhs", LOG_ENV => "results_directory")
        p = perturbation_params(args)
        integrator = (p=p, u=u, t=0.0)
        @test cb.affect!(integrator) === nothing                 # before initialize: no state
        cb.initialize(cb, u, 0.0, integrator)
        st = p.shared_buffers.gram_density_perturbation[]
        @test st.mode === :naive_rhs
        @test st.log_path == joinpath(dir, "gram_density_perturbation_log.csv")
        @test st.pass_count == [1, 1]                             # both start inside EI
        cb.affect!(integrator)
        @test st.pass_count == [1, 1]                             # no new entry
        # The mean model is not GRAM: the negative control leaves rho alone.
        wind = SVector(1.0, 2.0, 3.0)
        @test CB._apply_gram_density_perturbation(p, 1, 0.0, 100e3, 2.0, 200.0, wind) == (2.0, 200.0, wind)
        @test isempty(st.log_kind)
        cb.finalize(cb, u, 0.0, integrator)
        @test isfile(st.log_path)
        @test isfile(st.log_path * ".summary.toml")

        # The walk modes need a native GRAM model per spacecraft.
        for mode in ("step", "pass")
            cb_walk = perturbation_callback(args, MODE_ENV => mode)
            @test_throws ArgumentError cb_walk.initialize(cb_walk, u, 0.0, (p=perturbation_params(args), u=u, t=0.0))
        end
    end
end

# Installs a state with stub walks where initialize would have cloned native ones,
# attached to `p` as a run's first initialization at `t0` attaches it.
function install_stub_state!(cb, p, mode::Symbol, log_path::String; pass_dt=1.0, pass_max=5.0,
                             reseed::Bool=true, t0::Float64=0.0)
    st = CB._new_gram_density_perturbation_state(mode, N_SATS, p.args.environment_model.EI * 1e3,
                                                 pass_dt, pass_max, log_path, reseed)
    for i in 1:N_SATS
        st.walk_models[i] = StubWalk(0, 0)
        st.base_seeds[i] = 1000 + i
    end
    st.owner = p
    st.t0 = t0
    st.last_t = t0
    cb.affect!.state_ref[] = st
    p.shared_buffers.gram_density_perturbation[] = st
    return st
end

@testset "step mode: reseed at entry, hold the sampled factor, release above EI" begin
    mktempdir() do dir
        args = perturbation_config(dir)
        u = SE.build_initial_conditions(args)
        cb = perturbation_callback(args, MODE_ENV => "step")
        p = perturbation_params(args)
        integrator = (p=p, u=u, t=3.0)
        st = install_stub_state!(cb, p, :step, joinpath(dir, "step.csv"))

        p.shared_buffers.density_sample_t .= 0.0
        cb.affect!(integrator)
        @test st.pass_count == [1, 1]
        @test [w.seed for w in st.walk_models] == [CB._gram_pass_seed(1000 + i, 1) for i in 1:N_SATS]
        @test st.held_r ≈ [1.1, 1.1]
        @test st.walk_calls == [1, 1]
        @test all(isnan, p.shared_buffers.density_sample_t)      # staged samples invalidated
        @test CB._apply_gram_density_perturbation(p, 2, 3.0, 100e3, 2.0, 200.0, SVector(0.0, 0.0, 0.0))[1] ≈ 2.2

        p.shared_buffers.density_sample_t .= 0.0
        cb.affect!(integrator)
        @test st.held_r ≈ [1.2, 1.2]
        @test st.pass_count == [1, 1]
        @test all(isnan, p.shared_buffers.density_sample_t)

        # Leave the atmosphere: the factor drops to exactly 1 and no draw is made.
        st.ei_m = 0.0
        cb.affect!(integrator)
        @test st.held_r == [1.0, 1.0]
        @test st.walk_calls == [2, 2]
        @test st.in_atm_prev == [false, false]
        p.shared_buffers.density_sample_t .= 0.0
        cb.affect!(integrator)                                    # unchanged factor: no invalidation
        @test p.shared_buffers.density_sample_t == [0.0, 0.0]

        # Re-entry is pass 2, with its own seed.
        st.ei_m = 1000e3
        cb.affect!(integrator)
        @test st.pass_count == [2, 2]
        @test [w.seed for w in st.walk_models] == [CB._gram_pass_seed(1000 + i, 2) for i in 1:N_SATS]
        @test count(==(1), st.log_kind) == 6

        cb.finalize(cb, u, 3.0, integrator)
        @test length(readlines(joinpath(dir, "step.csv"))) == 7
    end
end

@testset "pass mode: predict and sweep at entry, interpolate, release above EI" begin
    mktempdir() do dir
        args = perturbation_config(dir)
        u = SE.build_initial_conditions(args)
        cb = perturbation_callback(args, MODE_ENV => "pass")
        p = perturbation_params(args)
        t0 = 7.0
        integrator = (p=p, u=u, t=t0)
        st = install_stub_state!(cb, p, :pass, joinpath(dir, "pass.csv"); pass_dt=1.0, pass_max=5.0)

        cb.affect!(integrator)
        # The drag-free prediction never leaves a 1000 km atmosphere, so the sweep
        # runs to the cap: ceil(5 / 1) + 1 knots, drawn in order.
        @test st.pass_active == [true, true]
        @test st.pass_t0 == [t0, t0]
        @test all(length.(st.pass_r) .== 6)
        @test st.pass_r[1] ≈ [1.0 + 0.1k for k in 1:6]
        @test st.walk_calls == [6, 6]
        @test all(a -> 0.0 < a < 1000e3, st.pass_alt[1])
        @test [w.seed for w in st.walk_models] == [CB._gram_pass_seed(1000 + i, 1) for i in 1:N_SATS]
        @test CB._gram_pass_factor(st, 1, t0 + 2.5) ≈ 1.35
        @test CB._apply_gram_density_perturbation(p, 1, t0 + 2.5, 100e3, 2.0, 200.0, SVector(0.0, 0.0, 0.0))[1] ≈ 2.7
        @test count(==(2), st.log_kind) == 12
        @test count(==(3), st.log_kind) == 2                     # applied factor at the entry step

        # A later accepted step inside the same pass reuses the prediction.
        cb.affect!((p=p, u=u, t=t0 + 1.0))
        @test st.walk_calls == [6, 6]
        @test count(==(3), st.log_kind) == 4

        # Climbing out deactivates the interpolant; the factor is 1 again.
        st.ei_m = 0.0
        cb.affect!(integrator)
        @test st.pass_active == [false, false]
        @test CB._gram_pass_factor(st, 1, t0 + 2.5) == 1.0

        # A prediction that climbs above EI stops at the first such knot.
        st.ei_m = 1000e3
        st.pass_max_s = 50.0
        cb.affect!(integrator)                                    # enter pass 2
        st.pass_active .= false
        st.ei_m = 0.0
        CB._build_gram_pass!(st, p, 1, SVector(7.0e6, 0.0, 0.0), SVector(0.0, 7.5e3, 0.0), t0)
        @test length(st.pass_r[1]) == 2
        # ... or reaches the surface.
        st.ei_m = 1000e3
        CB._build_gram_pass!(st, p, 1, SVector(1.0e3, 0.0, 0.0), SVector(0.0, 0.0, 0.0), t0)
        @test length(st.pass_r[1]) == 1

        cb.finalize(cb, u, t0, integrator)
        @test isfile(joinpath(dir, "pass.csv"))
    end
end

# Drives one run through accepted steps `(t, inside the atmosphere)`. After the
# steps listed in `boundaries` a checkpoint segment ends: its finalize, then the
# next segment's initialize at the same time, as the engine's checkpoint loop and
# the solver's callback reinitialization call them.
function drive_run!(cb, p, u, schedule; boundaries=Int[])
    for (k, (t, inside)) in enumerate(schedule)
        integrator = (p=p, u=u, t=t)
        p.shared_buffers.gram_density_perturbation[].ei_m = inside ? 1000e3 : 0.0
        cb.affect!(integrator)
        if k in boundaries
            cb.finalize(cb, u, t, integrator)
            cb.initialize(cb, u, t, integrator)
        end
    end
    t_end = last(schedule)[1]
    cb.finalize(cb, u, t_end, (p=p, u=u, t=t_end))
    return p.shared_buffers.gram_density_perturbation[]
end

const CONTINUED_FIELDS = (:pass_count, :walk_calls, :held_r, :pass_active, :pass_t0, :pass_r, :pass_alt,
                          :in_atm_prev, :base_seeds, :walk_fresh, :last_t, :log_kind, :log_sat, :log_pass,
                          :log_t, :log_alt, :log_r, :log_sigma, :log_mean, :log_aux)

@testset "checkpoint segments continue the run ($mode, reseed=$reseed)" for mode in (:step, :pass),
                                                                             reseed in (false, true)
    mktempdir() do dir
        args = perturbation_config(dir)
        u = SE.build_initial_conditions(args)
        # Three passes; segment boundaries inside pass 2 (t = 5) and between passes (t = 7).
        schedule = [(1.0, true), (2.0, true), (3.0, false), (4.0, true), (5.0, true), (6.0, true),
                    (7.0, false), (8.0, true), (9.0, true)]
        runs = map((:continuous, :segmented)) do kind
            cb = perturbation_callback(args, MODE_ENV => string(mode))
            p = perturbation_params(args)
            log = joinpath(dir, "$(kind).csv")
            st = install_stub_state!(cb, p, mode, log; reseed=reseed)
            final = drive_run!(cb, p, u, schedule; boundaries=kind === :segmented ? [5, 7] : Int[])
            @test final === st                                   # one state for the whole run
            (st=st, p=p, log=log)
        end
        a, b = runs
        for f in CONTINUED_FIELDS
            @test isequal(getfield(a.st, f), getfield(b.st, f))
        end
        @test [(w.seed, w.calls, w.first_updates) for w in a.st.walk_models] ==
              [(w.seed, w.calls, w.first_updates) for w in b.st.walk_models]
        @test b.st.pass_count == [3, 3]
        @test [w.seed for w in b.st.walk_models] ==
              (reseed ? [CB._gram_pass_seed(1000 + i, 3) for i in 1:N_SATS] : zeros(Int, N_SATS))
        wind = SVector(0.0, 0.0, 0.0)
        for t in (8.0, 8.5, 9.0)
            @test CB._apply_gram_density_perturbation(b.p, 1, t, 100e3, 1.0, 200.0, wind) ==
                  CB._apply_gram_density_perturbation(a.p, 1, t, 100e3, 1.0, 200.0, wind)
        end
        # The log written at the end holds every segment's rows.
        @test read(b.log, String) == read(a.log, String)
        @test read(b.log * ".summary.toml", String) == read(a.log * ".summary.toml", String)
        # Each walk's first update after its clone is marked, and with reseeding
        # the first update of every pass; no other update is.
        draws = mode === :step ? [2, 3, 2] : [6, 6, 6]
        expected = reduce(vcat, [vcat(reseed || k == 1, fill(false, n - 1)) for (k, n) in enumerate(draws)])
        @test all(w -> w.first_updates == expected, b.st.walk_models)
    end
end

@testset "checkpoint continuation: start time, run identity and staged densities" begin
    mktempdir() do dir
        args = perturbation_config(dir)
        u = SE.build_initial_conditions(args)
        cb = perturbation_callback(args, MODE_ENV => "naive_rhs")
        p = perturbation_params(args)
        cb.initialize(cb, u, 0.0, (p=p, u=u, t=0.0))
        st = p.shared_buffers.gram_density_perturbation[]
        @test st.owner === p && st.t0 == 0.0 && st.last_t == 0.0
        st.ei_m = 0.0
        cb.affect!((p=p, u=u, t=2.0))                             # leave the atmosphere
        st.ei_m = 1000e3
        cb.affect!((p=p, u=u, t=4.0))                             # pass 2
        @test st.pass_count == [2, 2] && st.last_t == 4.0

        # The next segment starts where this one ended: the same state, nothing
        # sampled, staged densities invalidated.
        p.shared_buffers.density_sample_t .= 0.0
        cb.initialize(cb, u, 4.0, (p=p, u=u, t=4.0))
        @test p.shared_buffers.gram_density_perturbation[] === st
        @test st.pass_count == [2, 2]
        @test all(isnan, p.shared_buffers.density_sample_t)

        # Neither the previous segment's end nor the run's start: refused, unchanged.
        @test_throws ArgumentError cb.initialize(cb, u, 3.0, (p=p, u=u, t=3.0))
        @test p.shared_buffers.gram_density_perturbation[] === st

        # The run's start again (a solve over again from the beginning): fresh.
        cb.initialize(cb, u, 0.0, (p=p, u=u, t=0.0))
        fresh = p.shared_buffers.gram_density_perturbation[]
        @test fresh !== st
        @test fresh.pass_count == [1, 1]

        # Another run's parameters never continue this run's state.
        p2 = perturbation_params(args)
        cb.initialize(cb, u, 0.0, (p=p2, u=u, t=0.0))
        other = p2.shared_buffers.gram_density_perturbation[]
        @test other !== fresh && other.owner === p2
        @test p.shared_buffers.gram_density_perturbation[] === fresh
    end
end

@testset "checkpoint segments through the solver continue an active pass" begin
    mktempdir() do dir
        args = perturbation_config(dir)
        u0 = SE.build_initial_conditions(args)
        stationary!(du, u, p, t) = fill!(du, 0.0)
        # The engine's checkpoint loop: one problem per segment, the same params
        # and callbacks; the next segment starts at the last one's final time.
        function run_segments(spans, cached::Bool, tag::String)
            p = perturbation_params(args)
            cb = perturbation_callback(args, MODE_ENV => "pass")
            st = install_stub_state!(cb, p, :pass, joinpath(dir, tag * ".csv"))
            cache = cached ? SE.SolverIntegratorCache() : nothing
            cfg = SE.SolverConfig(solver_mode=:tsit5)
            u = deepcopy(u0)
            for span in spans
                prob = SE.ODEProblem(stationary!, u, span, p; callback=CB.CallbackSet(cb))
                sol, _ = SE._solve_with_solver_policy(prob, cfg, args, 1e-8, 1e-8;
                                                      solver_cache=cache, needs_full_solution=false)
                @test SE.SciMLBase.successful_retcode(sol.retcode)
                @test sol.t[end] == span[2]
                u = deepcopy(sol.u[end])
            end
            return p, st
        end
        p_c, st_c = run_segments([(0.0, 2.0)], true, "continuous")
        @test st_c.pass_count == [1, 1] && st_c.pass_active == [true, true]
        for cached in (true, false)
            p_s, st_s = run_segments([(0.0, 1.0), (1.0, 2.0)], cached, "segmented_$(cached)")
            @test p_s.shared_buffers.gram_density_perturbation[] === st_s
            @test st_s.last_t == 2.0
            @test st_s.pass_count == st_c.pass_count
            @test st_s.pass_t0 == st_c.pass_t0
            @test st_s.pass_r == st_c.pass_r
            @test st_s.walk_calls == st_c.walk_calls
            @test [w.seed for w in st_s.walk_models] == [w.seed for w in st_c.walk_models]
            @test CB._gram_pass_factor(st_s, 1, 2.0) == CB._gram_pass_factor(st_c, 1, 2.0)
        end
    end
end

@testset "multirate rejects captured perturbation modes before its first subsolve" begin
    mktempdir() do dir
        args = perturbation_config(dir)
        cfg = SE.SolverConfig(solver_mode=:multirate, multirate_slow_dt_s=1.0,
            multirate_fast_substeps=2, multirate_slow_solver=:tsit5, multirate_fast_solver=:tsit5)
        for mode in ("step", "pass", "naive_rhs", "off"), in_set in (false, true)
            p = perturbation_params(args)
            u = SE.build_initial_conditions(args)
            rhs_calls = Ref(0)
            stationary!(du, u, p, t) = (rhs_calls[] += 1; fill!(du, 0.0))
            # Build under one mode, then solve with ENV changed to off. The
            # solver must validate the captured callback, not the current ENV.
            cb = perturbation_callback(args, MODE_ENV => mode)
            ordinary_cb = CB.DiscreteCallback((u, t, i) -> false, i -> nothing)
            callbacks = in_set ? (cb === nothing ? CB.CallbackSet(ordinary_cb) :
                CB.CallbackSet(ordinary_cb, cb)) : cb
            prob = SE.SplitODEProblem(stationary!, stationary!, u, (0.0, 2.0), p;
                                      callback=callbacks)
            withenv(MODE_ENV => "off") do
                if mode == "off"
                    sol, _ = SE._solve_with_solver_policy(prob, cfg, args, 1e-8, 1e-8)
                    @test SE.SciMLBase.successful_retcode(sol.retcode)
                    @test sol.t[end] == 2.0
                    @test rhs_calls[] > 0
                else
                    err = try
                        SE._solve_with_solver_policy(prob, cfg, args, 1e-8, 1e-8)
                        nothing
                    catch caught
                        caught
                    end
                    @test err isa ArgumentError
                    @test occursin("multirate does not support", sprint(showerror, err))
                    @test occursin("DENSITY_PERTURBATION=$(mode)", sprint(showerror, err))
                    @test rhs_calls[] == 0
                    @test cb.affect!.state_ref[] === nothing
                end
                @test p.shared_buffers.gram_density_perturbation[] === nothing
            end
        end
    end
end
end
