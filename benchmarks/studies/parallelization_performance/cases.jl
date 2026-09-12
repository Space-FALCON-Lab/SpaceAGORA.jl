# Load through the real, Pkg-loaded SpaceAGORA package (rather than raw
# `include`-ing the src/ files directly into Main, as this study used to)
# so package extensions actually attach -- in particular
# SpaceAGORAGRAMSuiteExt, which is what supplies the `GRAMAtmosphereModel`/
# `GRAMAtmosphereModelSurrogate` keyword constructors used by the GRAM-axis
# router-evaluation cases below. A raw-included copy of simulation_model.jl
# is a second, independent module instantiation with its own type identities
# that Julia's extension mechanism (keyed to the real SpaceAGORA Base.PkgId)
# cannot see, so those constructors never resolved under the old include-only
# setup even after `import GRAMSuite`. `using SpaceAGORA` still exposes every
# name this file (and modes.jl/execution.jl/trajectory_parity.jl) previously
# got from the raw includes, since SpaceAGORA.jl itself includes/re-exports
# the same source files.
using SpaceAGORA
using SpaceAGORA.SimulationModel
const SimulationEngine = SpaceAGORA.SimulationEngine
using SpaceAGORA.SimulationModel: GRAMAtmosphereModel, GRAMAtmosphereModelSurrogate

# Mirrors examples/common.jl's `ensure_gramsuite_loaded!` (that file isn't
# included here to avoid pulling in its Pkg.activate bootstrapping, which
# this study's cli.jl/worker subprocess launch already handles via
# `--project=$(PPC_REPO_ROOT)`). Idempotent: safe to call more than once.
function ppc_ensure_gramsuite_loaded!()
    isdefined(SpaceAGORA, :GRAMSuite) && return nothing
    vendored_gramsuite = joinpath(PPC_REPO_ROOT, "data", "GRAMSuite.jl")
    if Base.find_package("GRAMSuite") === nothing && isdir(vendored_gramsuite)
        pushfirst!(LOAD_PATH, vendored_gramsuite)
    end
    @eval import GRAMSuite
    return nothing
end

# Case names that build a live GRAMAtmosphereModel. Must be loaded here, at
# file-include time, rather than lazily inside those cases' branches in
# ppc_single_config below: `@eval import GRAMSuite` inside a function only
# advances the *global* world age, which is invisible to any frame already on
# the call stack above that `eval` without Base.invokelatest -- and GRAM's
# extension-provided methods are called from many places (the constructor,
# the density-callback's `getDensity`, ...) scattered too deep/far apart to
# invokelatest individually. Each worker subprocess handles exactly one
# --case=X for its whole lifetime (see ppc_worker_cmd), so ARGS already says
# up front whether this process will ever need GRAMSuite loaded.
const PPC_GRAM_LIVE_CASES = (
    "multi_16_gram_live", "multi_4_gram_live", "montecarlo_mars_gram_live",
    # Matches the whole heavy_<N>sat_gram_nbody_l50 ladder via `occursin` below,
    # so new satellite counts do not need adding here.
    "gram_nbody_l50",
    # Matches the B10 atmosphere-fidelity ladder's GRAM rungs
    # (atmo256_gram_surrogate_10min / atmo256_gram_live_10min /
    # atmo256_gram_live_nbody_10min) the same way. The surrogate rung needs
    # GRAMSuite loaded too: GRAMAtmosphereModelSurrogate's keyword constructor
    # builds a live GRAMAtmosphereModel first and then resolves the surrogate
    # file through GRAMSuite.gram_default_surrogate_file.
    "atmo256_gram",
    # Registering a case that builds a live GRAMAtmosphereModel means adding it
    # here as well as to the catalog. Miss this and the extension never attaches,
    # so the keyword constructor resolves to the stub in density_models.jl and
    # the worker dies with "no method matching GRAMAtmosphereModel(; planet_name)"
    # -- which reads like a broken case rather than a missing import.
    "campaign_gnc_6dof_gram",
)
any(a -> any(c -> occursin(c, a), PPC_GRAM_LIVE_CASES), ARGS) && ppc_ensure_gramsuite_loaded!()

# GRAMAtmosphereModel is cached per planet (not rebuilt per call/sample): the
# vendored GRAMSuite.jl itself dynamically (re)defines `set_library!` inside
# its own constructor path (data/GRAMSuite.jl/GRAM Suite 2.0/Julia/generic.jl),
# so a *second* construction for the same planet in one process hits the same
# world-age barrier as above ("MethodError: no method matching set_library!" /
# Julia's own "access to binding in a world prior to its definition" warning)
# -- independent of anything in this file or SpaceAGORAGRAMSuiteExt. This bites
# montecarlo_mars_gram_live specifically, since ppc_single_config is called
# fresh per MC sample within one worker process (ppc_run_sample_batch), so
# mc_samples > 1 would otherwise try to build a second Mars GRAMAtmosphereModel.
# Reusing one instance is also the officially-supported concurrent-access
# pattern (GRAMAtmosphereModel's docstring: `instance_lock` "serializes native
# GRAM calls against this wrapper instance"), so this cache is simultaneously
# the crash fix and the intended way to exercise GRAM contention under
# threaded outer-parallelism -- vs. outer_process, where each sample gets its
# own process and thus its own uncontended instance. That contrast (shared
# lock-serialized instance under threads vs. per-process isolation) is exactly
# the "native-library contention" axis the manuscript's execution-architecture
# contribution claims to route on.
const _PPC_GRAM_MODEL_CACHE = Dict{String, Any}()
const _PPC_GRAM_MODEL_CACHE_LOCK = ReentrantLock()
# `get!`'s do-block form is not thread-safe against concurrent *first*
# population: under outer_threads/full_smart with mc_samples > 1, every MC
# sample's ppc_single_config call lands on a different thread at roughly the
# same time (Threads.@threads in ppc_run_sample_batch), so with a plain Dict
# multiple threads can all see the entry missing and race to construct
# GRAMAtmosphereModel concurrently -- each hitting the same
# set_library!-in-the-wrong-world-age crash the cache above exists to avoid,
# just via a race instead of a second sequential call. Locking around the
# whole check-and-construct makes population itself serialize (matching the
# type's own instance_lock intent for concurrent *use*).
function ppc_gram_atmosphere_model(planet_name::String)
    lock(_PPC_GRAM_MODEL_CACHE_LOCK) do
        get!(_PPC_GRAM_MODEL_CACHE, planet_name) do
            GRAMAtmosphereModel(planet_name=planet_name)
        end
    end
end

const PPC_SPICE_PATH = joinpath(PPC_REPO_ROOT, "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE")
const PPC_EARTH_HARMONICS_FILE = joinpath(PPC_REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")

# Internal (non-exported) home of `_lvlh_cascade_torque`, reached the same way
# the CYGNSS closed-loop-twin driver reaches it
# (spaceagora-private-telemetry/CYGNSS/campaign_drivers/productized/cygnss_constellation_scaling_worker.jl),
# for the actuator/control-effector router-evaluation case below.
const PE = SimulationModel.DynamicEffectors.PerturbationEffectors

Base.@kwdef struct PPCCaseSpec
    name::String
    family::String
    description::String
    builder::Function
    montecarlo::Bool = false
    orientation::Bool = false
    default_samples::Int = 1
end

function ppc_spacecraft(
    planet;
    id::Int=1,
    ra_alt_m::Float64=550e3,
    rp_alt_m::Float64=500e3,
    i_deg::Float64=35.0,
    omega_deg::Float64=40.0,
    raan_deg::Float64=10.0,
    nu_deg::Float64=170.0,
    with_panel::Bool=false,
    panel_count::Int=1,
    root_mass::Float64=500.0,
    root_area::Float64=12.0,
    prop_mass::Float64=0.0,
    orientation_state=nothing
)
    root = Link(root=true, m=root_mass, ref_area=root_area)
    links = Link[root]
    if with_panel
        for panel_idx in 1:panel_count
            theta = 2π * (panel_idx - 1) / max(1, panel_count)
            panel = Link(
                root=false,
                m=8.0,
                ref_area=3.0,
                r=MVector{3, Float64}(1.8 * cos(theta), 1.8 * sin(theta), 0.0)
            )
            push!(links, panel)
        end
    end
    ra = planet.Rp_e + ra_alt_m
    rp = planet.Rp_e + rp_alt_m
    ic = if orientation_state === nothing
        InitialCondition(ra=ra, rp=rp, i=i_deg, ω=omega_deg, Ω=raan_deg, ν=nu_deg)
    else
        q0, w0 = orientation_state
        a = (ra + rp) / 2.0
        e = (ra - rp) / (ra + rp)
        InitialCondition(a, e, i_deg, omega_deg, raan_deg, nu_deg, q0, w0)
    end
    dry_mass = sum(link.m for link in links)
    return SpacecraftModel(Joint[], links, root, true, dry_mass, prop_mass, root.inertia, 0, 0, ic, id)
end

function ppc_constellation(
    planet,
    n::Int;
    with_panel::Bool=false,
    panel_count::Int=1,
    orientation_state=nothing,
)
    sats = SpacecraftModel[]
    for i in 1:n
        push!(sats, ppc_spacecraft(
            planet;
            id=i,
            ra_alt_m=540e3 + 2e3 * (i - 1),
            rp_alt_m=500e3 + 1e3 * (i - 1),
            nu_deg=120.0 + 240.0 * (i - 1) / max(1, n),
            with_panel=with_panel,
            panel_count=panel_count,
            orientation_state=orientation_state
        ))
    end
    return sats
end

function ppc_harmonics_model(planet, degree::Int)
    if isfile(PPC_EARTH_HARMONICS_FILE)
        return GravitationalHarmonicsModel(degree, degree, PPC_EARTH_HARMONICS_FILE, planet)
    end
    return InverseSquaredJ2GravityModel()
end

function ppc_build_config(;
    planet,
    spacecraft,
    mission_time_s::Float64,
    orientation_sim::Bool,
    dynamic_effectors::Tuple,
    density_model=NoAtmosphereModel(),
    guidance_effectors::Tuple=(),
    guidance_rates::Vector{Float64}=Float64[],
    control_effectors::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
    dt_max_orbit::Float64=10.0,
    reltol_orbit::Float64=1e-9,
    abstol_orbit::Float64=1e-9,
    num_steps_to_save::Int=300,
    results::Bool=false,
    data_rate::Float64=10.0
)
    # `results` and `data_rate` are the output-cadence axis (B14).
    #
    # `num_steps_to_save` above is NOT that axis and never was: nothing outside
    # analysis/verification/telemetry_verification/ reads it, and the solver runs
    # on save_everystep with no saveat, so varying it changes literally nothing
    # about the work performed. The real mechanism is
    # get_data_saving_callback (src/simulation/callbacks/event_callbacks.jl),
    # which builds `SavingCallback(...; saveat=data_rate)` -- but it is only
    # installed when `simulation_settings.results` is true
    # (src/simulation/callbacks/density_callbacks/assembly.jl), and this harness
    # hardcoded results=false, so every benchmark case to date has done no
    # trajectory output at all. Setting results=true both installs that callback
    # and (via `needs_full_solution` in simulation/engine/execution.jl) turns on
    # per-step solution storage, so the cadence rungs must all share the same
    # `results` value for `data_rate` to be the only thing that varies between
    # them.
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=results,
            verbose=false,
            generate_plots=false,
            normalize=false,
            save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=mission_time_s,
            orientation_sim=orientation_sim,
            num_steps_to_save=num_steps_to_save,
            data_rate=data_rate
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=density_model,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel(spacecraft, dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=guidance_effectors, guidance_rates=guidance_rates),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=control_effectors, control_rates=control_rates),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=reltol_orbit,
            abstol_orbit=abstol_orbit,
            dt_max_orbit=dt_max_orbit
        )
    )
end

@inline function ppc_mission_time(
    profile::String;
    test::Float64=10.0,
    smoke::Float64=120.0,
    full::Float64=1800.0
)::Float64
    profile == "test" && return test
    profile == "smoke" && return smoke
    return full
end

function ppc_single_config(case_name::String, cfg::PPCConfig; seed::Int=cfg.seed, mc_index::Int=1)
    planet = Earth("", PPC_SPICE_PATH)
    mars = Mars("", PPC_SPICE_PATH)
    q0 = normalize(SVector{4, Float64}(0.15, -0.05, 0.2, 0.96))
    w0 = SVector{3, Float64}(0.01, -0.02, 0.015)
    rng = MersenneTwister(seed + 1009 * mc_index)

    if case_name == "single_inverse_square_vacuum"
        return ppc_build_config(
            planet=planet,
            spacecraft=[ppc_spacecraft(planet)],
            mission_time_s=ppc_mission_time(cfg.profile),
            orientation_sim=false,
            dynamic_effectors=(InverseSquaredGravityModel(),)
        )
    elseif occursin(r"^gravity_[0-9]+sat_inverse_square_vacuum$", case_name)
        n = parse(Int, match(r"^gravity_([0-9]+)sat", case_name).captures[1])
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, n),
            mission_time_s=ppc_mission_time(cfg.profile),
            orientation_sim=false,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=20.0
        )
    elseif case_name == "single_harmonics_l20_vacuum"
        return ppc_build_config(
            planet=planet,
            spacecraft=[ppc_spacecraft(planet)],
            mission_time_s=ppc_mission_time(cfg.profile),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 20),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=10.0
        )
    elseif occursin(r"^gravity_[0-9]+sat_l20_vacuum$", case_name)
        n = parse(Int, match(r"^gravity_([0-9]+)sat", case_name).captures[1])
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, n),
            mission_time_s=ppc_mission_time(cfg.profile),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 20),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=20.0
        )
    elseif occursin(r"^gravity_[0-9]+sat_l50_vacuum_1hr$", case_name)
        # Ad hoc spacecraft-count x thread-count sweep case: L50 harmonics (heavier
        # per-satellite RHS cost than L20/inverse-square, to give outer-thread
        # parallelism the most room to show a real effect) over a fixed 1-hour
        # simulated mission, at N parameterized by the case name.
        n = parse(Int, match(r"^gravity_([0-9]+)sat", case_name).captures[1])
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, n),
            mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=300.0, full=3600.0),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 50),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=20.0
        )
    elseif case_name == "single_harmonics_l50_vacuum" || case_name == "montecarlo_high_accuracy"
        ra_jitter = case_name == "montecarlo_high_accuracy" ? randn(rng) * 8e3 : 0.0
        rp_jitter = case_name == "montecarlo_high_accuracy" ? randn(rng) * 8e3 : 0.0
        return ppc_build_config(
            planet=planet,
            spacecraft=[ppc_spacecraft(planet; ra_alt_m=550e3 + ra_jitter, rp_alt_m=500e3 + rp_jitter)],
            mission_time_s=ppc_mission_time(cfg.profile),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 50),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=5.0,
            reltol_orbit=1e-10,
            abstol_orbit=1e-10
        )
    elseif case_name == "srp_heavy_high_area"
        return ppc_build_config(
            planet=planet,
            spacecraft=[ppc_spacecraft(planet; with_panel=true, panel_count=8, root_area=40.0, orientation_state=(q0, w0))],
            mission_time_s=ppc_mission_time(cfg.profile),
            orientation_sim=true,
            dynamic_effectors=(InverseSquaredGravityModel(), SolarRadiationPressureModel(1.2, 120.0), SolarRadiationPressureModel(1.8, 220.0)),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=2.0
        )
    elseif case_name == "articulated_1sat_fullstack"
        return ppc_build_config(
            planet=planet,
            spacecraft=[ppc_spacecraft(planet; with_panel=true, panel_count=28, orientation_state=(q0, w0), root_area=10.0)],
            mission_time_s=ppc_mission_time(cfg.profile; smoke=90.0, full=1200.0),
            orientation_sim=true,
            dynamic_effectors=(ppc_harmonics_model(planet, 20), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=1.0
        )
    elseif case_name == "multi_16_aero_surrogate_cached"
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 16),
            mission_time_s=ppc_mission_time(cfg.profile; smoke=90.0, full=1200.0),
            orientation_sim=false,
            dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=5.0
        )
    elseif case_name == "multi_64_high_fidelity" || case_name == "multi_128_high_fidelity"
        n = case_name == "multi_128_high_fidelity" ? 128 : 64
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, n),
            mission_time_s=ppc_mission_time(cfg.profile; smoke=90.0, full=1200.0),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 20), SolarRadiationPressureModel(1.2, 12.0), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=10.0
        )
    elseif case_name == "callback_128_aero_thermal"
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 128),
            mission_time_s=ppc_mission_time(cfg.profile; smoke=90.0, full=1200.0),
            orientation_sim=false,
            dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=3.0
        )
    elseif case_name == "montecarlo_multi_sat"
        spacecraft = SpacecraftModel[]
        for i in 1:4
            push!(spacecraft, ppc_spacecraft(
                planet;
                id=i,
                ra_alt_m=540e3 + randn(rng) * 6e3,
                rp_alt_m=500e3 + randn(rng) * 6e3,
                nu_deg=130.0 + 45.0 * i + randn(rng) * 2.0
            ))
        end
        return ppc_build_config(
            planet=planet,
            spacecraft=spacecraft,
            mission_time_s=ppc_mission_time(cfg.profile),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 20),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=10.0
        )
    elseif case_name == "montecarlo_mars_aerobraking"
        return ppc_build_config(
            planet=mars,
            spacecraft=[ppc_spacecraft(
                mars;
                ra_alt_m=4500e3 + randn(rng) * 100e3,
                rp_alt_m=max(110e3, 135e3 + randn(rng) * 10e3),
                i_deg=93.0,
                omega_deg=80.0,
                raan_deg=30.0,
                nu_deg=180.0 + randn(rng) * 4.0
            )],
            mission_time_s=ppc_mission_time(cfg.profile; smoke=120.0, full=1800.0),
            orientation_sim=false,
            dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(mars),
            dt_max_orbit=1.0
        )

    # ── Router-evaluation axis coverage (point 8: spacecraft count, atmosphere/GRAM
    # usage, force/actuator model count, interacting vs. independent propagation,
    # thread/process budgets, output cadence and mission duration) ─────────────────

    elseif case_name == "multi_4_aero_surrogate_cached"
        # Rounds out the many_sat_high_fidelity spacecraft-count ladder at the low
        # end (existing: 16/64/128) for the router phase's spacecraft-count axis.
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 4),
            mission_time_s=ppc_mission_time(cfg.profile; smoke=90.0, full=1200.0),
            orientation_sim=false,
            dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=5.0
        )
    elseif case_name == "multi_256_high_fidelity"
        # High end of the same ladder (existing: 64/128).
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 256),
            mission_time_s=ppc_mission_time(cfg.profile; smoke=90.0, full=1200.0),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 20), SolarRadiationPressureModel(1.2, 12.0), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=10.0
        )
    elseif case_name == "multi_16_gram_live"
        # Interacting constellation against the *native* GRAM atmosphere model
        # (SpaceAGORAGRAMSuiteExt's GRAMAtmosphereModel), not the analytic
        # ExponentialAtmosphereModel every other atmosphere case in this catalog
        # uses. Exercises real GRAM-call contention (native-library mutex/cache
        # traffic) under outer-loop constellation parallelism -- the "atmosphere
        # and GRAM usage" axis, on the interacting side. GRAMSuite is already
        # loaded (see PPC_GRAM_LIVE_CASES above the case catalog) -- world-age
        # safe because that happens before this function is even compiled.
        #
        # KNOWN ISSUE (2026-08-18): at the real full mission duration (1200s),
        # this case leaks memory unboundedly -- observed growing ~150-300 MB/min
        # with no plateau, eventually OOM-killing the host. Reproduced with
        # solver_mode=tsit5 too (rules out the auto_stiff/Rodas5P autoswitch
        # path), and the MERRA2 data file is read exactly once (rules out a
        # repeated-file-read cause), so this looks like a per-call resource leak
        # inside the vendored GRAMSuite.jl native binding itself, scaling with
        # GRAM call volume (16 satellites x many RHS evaluations). The
        # single-satellite montecarlo_mars_gram_live case below does NOT show
        # this (35s, clean) at its own real full mission duration (1800s), so
        # it's specific to this case, not live GRAM in general. Deliberately
        # excluded from B6 (paper_parallelization_benchmarks/cli.jl) until this
        # is root-caused. Do not re-add it to a phase without re-validating at
        # the real (not "test"-profile) mission duration first -- a short smoke
        # run will not reproduce this.
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 16),
            mission_time_s=ppc_mission_time(cfg.profile; smoke=90.0, full=1200.0),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 20), AerodynamicCoefficientfM()),
            density_model=ppc_gram_atmosphere_model("earth"),
            dt_max_orbit=5.0
        )
    elseif case_name == "multi_4_gram_live"
        # Added to check whether multi_16_gram_live's leak (see its own comment
        # above) is proportional to satellite count, as a possible small-N
        # stopgap for the interacting side of the atmosphere-GRAM axis. Result
        # (2026-08-18, bounded OS-timeout + memory-watched reproduction at the
        # real full mission duration, 1200s): memory does NOT leak here -- RSS
        # held in the 1.6-1.9 GB range with no runaway growth, unlike N=16's
        # unbounded climb to 30+ GB. But it's still far too SLOW to be a
        # practical stopgap: it did not finish a single solve within a 10-minute
        # hard timeout (vs. ~35s for the single-satellite montecarlo_mars_gram_live
        # at its own real 1800s mission). So this is a *different* failure mode
        # from N=16's leak, not a smaller instance of the same one -- rules out
        # "proportional to N" as the leak's cause, but doesn't give us a usable
        # interacting-GRAM case either.
        #
        # Root-caused (2026-08-18) via solver statistics (sol.stats) at short
        # mission durations, plain solver_mode=tsit5 (rules out auto_stiff/
        # Rodas5P as the cause): a single Earth satellite against this same
        # live GRAMAtmosphereModel needs only 50 accepted steps / 352 f_evals
        # for a 10s mission -- completely healthy. This 4-satellite case needs
        # 14,540 accepted steps / 101,836 f_evals for the same 10s window
        # (average step ~0.7ms, ~7000x smaller than dt_max_orbit=5.0's cap).
        # That ~290x blowup from 4x more satellites rules out "Earth GRAM /
        # MERRA2 data is inherently noisy" (N=1 is clean) and points instead
        # at the multi-satellite *sharing* of one GRAMAtmosphereModel instance:
        # density queries for different satellites interleave within a single
        # RHS evaluation, and if the native library carries call-order-
        # dependent internal state (consistent with the set_library!
        # single-construction-only bug documented on ppc_gram_atmosphere_model
        # above -- the same underlying "not designed for concurrent users"
        # class of problem), satellite 2's query could be subtly influenced by
        # satellite 1's preceding one. A step-adaptive integrator reads that as
        # extreme non-smoothness and shrinks its step chasing noise it can
        # never resolve. Plausible fix -- one GRAMAtmosphereModel instance per
        # satellite -- is exactly what the world-age bug currently forbids, so
        # this and multi_16_gram_live's leak look like two symptoms of the same
        # upstream GRAMSuite limitation rather than two independent bugs. NOT
        # wired into any phase.
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 4),
            mission_time_s=ppc_mission_time(cfg.profile; smoke=90.0, full=1200.0),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 20), AerodynamicCoefficientfM()),
            density_model=ppc_gram_atmosphere_model("earth"),
            dt_max_orbit=5.0
        )
    elseif case_name == "campaign_gnc_6dof_gram"
        # Everything at once: embedded GNC, coupled 6-DOF attitude, a live native
        # atmosphere, and campaign-scale propagation. No other case in this
        # catalog carries all four -- montecarlo_mars_gram_live has the campaign
        # and the native atmosphere but is 3-DOF and uncontrolled,
        # heavy_256sat_coupled6dof_2hr has the attitude but an analytic
        # atmosphere and one sample, and stack<N>_e6_actuated has the GNC but
        # neither GRAM nor a campaign. This case exists to time the combination.
        #
        # Earth rather than Mars because the actuator is a magnetorquer, which
        # needs a global magnetic field to push against; Mars has none, so the
        # Mars aerobraking geometry the other GRAM campaign uses would make the
        # control loop physically meaningless.
        #
        # Four spacecraft per sample, not more: the control-effector rung is
        # quadratic in N (one MagneticMomentumManagerModel per satellite, each
        # evaluating a commanded torque), and every density query on top of that
        # is a native GRAM call under a global lock. The campaign axis carries
        # the scale here, which is the point -- GRAM is not thread-safe, so the
        # samples are what the process route can actually distribute.
        gnc_n = 4
        gnc_harmonics = ppc_harmonics_model(planet, 20)
        gnc_srp = SolarRadiationPressureModel(1.2, 12.0)
        gnc_spacecraft = ppc_constellation(
            planet, gnc_n;
            with_panel=true, panel_count=4, orientation_state=(q0, w0),
        )
        gnc_controller = PE.LVLHCascadeAttitudeControlModel(
            q_cmd_lb=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
            k_out=SVector{3, Float64}(0.05, 0.05, 0.05),
            w_max=0.01,
            k_rate=SVector{3, Float64}(50.0, 50.0, 50.0),
            tau_max=0.01,
        )
        gnc_torque = (t, r, v, q, w) -> PE._lvlh_cascade_torque(gnc_controller, r, v, q, w)
        gnc_l_pi = MMatrix{3, 3, Float64}(1.0I)
        gnc_b_field = (t, r_ii) -> PE.get_magnetic_field_dipole(r_ii, gnc_l_pi)
        gnc_mgr = MagneticMomentumManagerModel[
            MagneticMomentumManagerModel(
                sat_idx=i, mu_gain=7.0e-4, m_max_am2=1.8,
                h_wheels_0=SVector{3, Float64}(0.0, 0.0, 0.0),
                commanded_torque=gnc_torque, b_field_ii=gnc_b_field
            )
            for i in 1:gnc_n
        ]
        return ppc_build_config(
            planet=planet,
            spacecraft=gnc_spacecraft,
            mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=120.0, full=600.0),
            orientation_sim=true,
            dynamic_effectors=(gnc_harmonics, gnc_srp, AerodynamicCoefficientfM(), gnc_controller),
            density_model=ppc_gram_atmosphere_model("earth"),
            control_effectors=Tuple(gnc_mgr),
            control_rates=fill(1.0, gnc_n),
            dt_max_orbit=10.0
        )
    elseif case_name == "montecarlo_mars_gram_live"
        # Independent-trial counterpart to multi_16_gram_live: same Mars
        # aerobraking geometry as montecarlo_mars_aerobraking, but native GRAM
        # density instead of ExponentialAtmosphereModel. This is also the
        # workload the review's point 7 flagged the router choosing poorly on
        # (R5 ~25% slower than the best tested route for the GRAM case) --
        # giving it its own case here (rather than only the surrogate/analytic
        # version) lets the expanded thread/process ladder actually probe that.
        return ppc_build_config(
            planet=mars,
            spacecraft=[ppc_spacecraft(
                mars;
                ra_alt_m=4500e3 + randn(rng) * 100e3,
                rp_alt_m=max(110e3, 135e3 + randn(rng) * 10e3),
                i_deg=93.0,
                omega_deg=80.0,
                raan_deg=30.0,
                nu_deg=180.0 + randn(rng) * 4.0
            )],
            mission_time_s=ppc_mission_time(cfg.profile; smoke=120.0, full=1800.0),
            orientation_sim=false,
            dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
            density_model=ppc_gram_atmosphere_model("mars"),
            dt_max_orbit=1.0
        )
    elseif case_name == "multi_8sat_magnetorquer_attitude"
        # Every other case in this catalog has control_effectors=() -- no case
        # exercises the discrete-rate control-effector path (stateful per-satellite
        # accumulator on the ControlModel callback, distinct from the continuous
        # dynamic_effectors RHS terms) at constellation scale. This case adds a
        # real actuator: an LVLH cascade attitude controller (dynamic_effectors,
        # shared across satellites -- see _lvlh_cascade_torque, which is
        # satellite-state-driven rather than sat_idx-keyed) commanding per-satellite
        # MagneticMomentumManagerModel magnetorquer unloading (control_effectors,
        # one stateful instance per sat_idx). The controller gains below are generic
        # small-sat-scale values for exercising this code path under routing, NOT
        # flight-fit numbers -- unlike the CYGNSS closed-loop-twin driver in
        # spaceagora-private-telemetry, this case makes no flight-fidelity claim.
        # b_field_ii reuses the library's own tilted-dipole model
        # (get_magnetic_field_dipole) with a fixed (non-rotating) ECEF frame --
        # a deliberate simplification appropriate for a routing benchmark over
        # missions this short (tens of minutes), not a flight-accurate B-field.
        n_sats = 8
        spacecraft = SpacecraftModel[
            ppc_spacecraft(
                planet;
                id=i,
                ra_alt_m=540e3 + 2e3 * (i - 1),
                rp_alt_m=500e3 + 1e3 * (i - 1),
                nu_deg=120.0 + 240.0 * (i - 1) / n_sats,
                orientation_state=(q0, w0)
            )
            for i in 1:n_sats
        ]
        controller = PE.LVLHCascadeAttitudeControlModel(
            q_cmd_lb=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
            k_out=SVector{3, Float64}(0.05, 0.05, 0.05),
            w_max=0.01,
            k_rate=SVector{3, Float64}(50.0, 50.0, 50.0),
            tau_max=0.01,
        )
        commanded_torque = (t, r, v, q, w) -> PE._lvlh_cascade_torque(controller, r, v, q, w)
        l_pi_fixed = MMatrix{3, 3, Float64}(1.0I)
        b_field_ii = (t, r_ii) -> PE.get_magnetic_field_dipole(r_ii, l_pi_fixed)
        mgr_models = MagneticMomentumManagerModel[
            MagneticMomentumManagerModel(
                sat_idx=i, mu_gain=7.0e-4, m_max_am2=1.8,
                h_wheels_0=SVector{3, Float64}(0.0, 0.0, 0.0),
                commanded_torque=commanded_torque, b_field_ii=b_field_ii
            )
            for i in 1:n_sats
        ]
        return ppc_build_config(
            planet=planet,
            spacecraft=spacecraft,
            mission_time_s=ppc_mission_time(cfg.profile; smoke=90.0, full=1200.0),
            orientation_sim=true,
            dynamic_effectors=(InverseSquaredGravityModel(), controller),
            density_model=NoAtmosphereModel(),
            guidance_effectors=(),
            control_effectors=Tuple(mgr_models),
            control_rates=fill(1.0, n_sats),
            dt_max_orbit=2.0
        )
    elseif case_name == "gravity_16sat_l20_vacuum_longmission"
        # Mission-duration axis: same physics as gravity_16sat_l20_vacuum, ~4x the
        # standard "full" mission length.
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 16),
            mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=480.0, full=7200.0),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 20),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=20.0
        )
    # ── B10: atmosphere / GRAM fidelity ladder ────────────────────────────────
    #
    # Five rungs at a fixed N=256 constellation, fixed 600 s mission and fixed L50
    # harmonics. The ONLY thing that varies across the ladder is how atmospheric
    # density is obtained, so any routing difference between rungs is attributable
    # to the density path rather than to problem size. This is the axis review
    # point 7 singles out (R5 ~25% slower than the best tested route, and slower
    # than serial, on the GRAM workload), so it carries the most weight of the six.
    #
    # 600 s rather than the 1 h the vacuum ladder uses: that is the duration the
    # heavy_<N>sat_gram_nbody_l50 family is already validated at for N up to 256,
    # and live GRAM at N=256 is by far the slowest thing in this catalog. Keeping
    # every rung at the same duration is what makes them comparable.
    #
    # The rungs, in increasing order of native-library serialisation:
    #   vacuum       NoAtmosphereModel, no aero effector -- the floor.
    #   exponential  analytic ExponentialAtmosphereModel + aero. Adds the density
    #                callback and drag, with no native library involved.
    #   gram_surrogate  GRAMAtmosphereModelSurrogate: precomputed GRAM table read
    #                from <gram_root>/<Planet>/<planet>_surrogate.jls, so the
    #                atmosphere is GRAM's but the queries do not enter native GRAM.
    #   gram_live    GRAMAtmosphereModel: real native GRAM calls, one shared
    #                instance serialised by its own instance_lock.
    #   gram_live_nbody  adds Sun/Moon NBodyGravityModel on top, i.e. a *second*
    #                serialised native library (CSPICE under SPICE_LOCK). §8.4 of
    #                docs/architecture/parallelization_measurement_and_allocation_findings.md
    #                measured this configuration as the one that does produce real
    #                thread scaling where the smaller GRAM cases were flat.
    elseif startswith(case_name, "atmo256_")
        atmo_n = 256
        atmo_time = ppc_mission_time(cfg.profile; test=10.0, smoke=90.0, full=600.0)
        atmo_harmonics = ppc_harmonics_model(planet, 50)
        atmo_nbody = NBodyGravityModel(
            body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet
        )
        density, effectors = if case_name == "atmo256_vacuum_10min"
            (NoAtmosphereModel(), (atmo_harmonics,))
        elseif case_name == "atmo256_exponential_10min"
            (ExponentialAtmosphereModel(planet), (atmo_harmonics, AerodynamicCoefficientfM()))
        elseif case_name == "atmo256_gram_surrogate_10min"
            (GRAMAtmosphereModelSurrogate(planet_name="earth"),
             (atmo_harmonics, AerodynamicCoefficientfM()))
        elseif case_name == "atmo256_gram_live_10min"
            (ppc_gram_atmosphere_model("earth"), (atmo_harmonics, AerodynamicCoefficientfM()))
        elseif case_name == "atmo256_gram_live_nbody_10min"
            (ppc_gram_atmosphere_model("earth"),
             (atmo_harmonics, atmo_nbody, AerodynamicCoefficientfM()))
        else
            throw(ArgumentError("Unknown atmosphere-ladder case '$case_name'."))
        end
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, atmo_n),
            mission_time_s=atmo_time,
            orientation_sim=false,
            dynamic_effectors=effectors,
            density_model=density,
            dt_max_orbit=5.0
        )

    # ── B11: force- and actuator-model count ladder ───────────────────────────
    #
    # Six rungs at a fixed constellation size and fixed 600 s mission, each adding
    # exactly one model to the rung before it. Nothing else in this catalog sweeps
    # model count as an isolated variable -- the existing multi-effector cases each
    # pick one arbitrary stack -- yet it is the variable §7.1 of the findings doc
    # identifies as driving the heterogeneous-effector-tuple cost the flat work
    # queue has to schedule around.
    #
    #   e1_harm      L20 harmonics only                         (1 dynamic effector)
    #   e2_srp       + SolarRadiationPressureModel               (2)
    #   e3_aero      + AerodynamicCoefficientfM, exponential density (3, + density callback)
    #   e4_nbody     + Sun/Moon NBodyGravityModel                (4, + SPICE under lock)
    #   e5_6dof      same four effectors, orientation_sim=true   (attitude states on)
    #   e6_actuated  + LVLH cascade controller and per-satellite magnetorquer
    #                control effectors                           (5 dynamic + N control)
    #
    # e5 and e6 are split deliberately: e5-vs-e4 isolates turning attitude
    # propagation on, and e6-vs-e5 isolates adding the actuator/control-callback
    # path, so neither is confounded with the other. The controller construction
    # mirrors multi_8sat_magnetorquer_attitude below -- see that case's comment for
    # why these gains are routing-benchmark values and make no flight-fidelity
    # claim.
    #
    # SATELLITE COUNT IS IN THE CASE NAME, and the ladder is deliberately run at a
    # smaller N than the other axes, because the actuator rung is quadratic in N.
    # get_control_callbacks (src/simulation/callbacks/control_callbacks.jl) builds
    # one PeriodicCallback per control effector and each one loops over *all*
    # num_sats, so N per-satellite MagneticMomentumManagerModels cost N^2
    # calcControlEffect! calls per tick -- of which N are the intended ones and
    # N(N-1) are calls with a sat_idx the model does not own. That is invisible at
    # the 8 satellites multi_8sat_magnetorquer_attitude uses (64 calls) and
    # dominant at 256 (65,536): measured post-warm-up, e6 at N=256 takes 33.4 s for
    # a *10 s* mission against 0.13 s for the otherwise-identical e5, which would
    # put a full-duration run in the hours.
    #
    # This is a real defect in the control path rather than an artifact of the
    # benchmark -- the shared-model-with-per-spacecraft-vectors pattern that
    # BaseThrusterModel enforces is the shape that avoids it -- but fixing it is
    # out of scope here. The ladder is sized so every rung is measurable instead,
    # and the quadratic itself is reported.
    elseif occursin(r"^stack[0-9]+_e[0-9]", case_name)
        stack_n = parse(Int, match(r"^stack([0-9]+)_", case_name).captures[1])
        # Duration is keyed off N because the actuator rung's N^2 term (above)
        # means one mission length cannot put both sub-ladders in a measurable
        # band. Measured post-warm-up per 10 s of simulated mission, serial:
        # e1 0.0009 s / e5 0.0157 s / e6 0.583 s at N=32, and e5 0.13 s at N=256.
        # At N>=128 a 1 h mission puts the model-count rungs at roughly 3-47 s;
        # at N<128 a 600 s mission puts the actuator rung at ~35 s. Each
        # sub-ladder is internally single-variable, which is what the axis needs.
        stack_time = stack_n >= 128 ?
            ppc_mission_time(cfg.profile; test=10.0, smoke=300.0, full=3600.0) :
            ppc_mission_time(cfg.profile; test=10.0, smoke=120.0, full=600.0)
        stack_harmonics = ppc_harmonics_model(planet, 20)
        stack_srp = SolarRadiationPressureModel(1.2, 12.0)
        stack_nbody = NBodyGravityModel(
            body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet
        )
        stack_rung = match(r"^stack[0-9]+_(e[0-9]+_[a-z0-9]+)$", case_name)
        stack_rung === nothing &&
            throw(ArgumentError("Unknown effector-ladder case '$case_name'."))
        stack_rung = stack_rung.captures[1]
        stack_attitude = stack_rung in ("e5_6dof", "e6_actuated")
        stack_effectors, stack_density = if stack_rung == "e1_harm"
            ((stack_harmonics,), NoAtmosphereModel())
        elseif stack_rung == "e2_srp"
            ((stack_harmonics, stack_srp), NoAtmosphereModel())
        elseif stack_rung == "e3_aero"
            ((stack_harmonics, stack_srp, AerodynamicCoefficientfM()),
             ExponentialAtmosphereModel(planet))
        elseif stack_rung in ("e4_nbody", "e5_6dof", "e6_actuated")
            ((stack_harmonics, stack_srp, stack_nbody, AerodynamicCoefficientfM()),
             ExponentialAtmosphereModel(planet))
        else
            throw(ArgumentError("Unknown effector-ladder case '$case_name'."))
        end
        stack_spacecraft = ppc_constellation(
            planet, stack_n;
            with_panel=stack_attitude,
            panel_count=stack_attitude ? 4 : 1,
            orientation_state=stack_attitude ? (q0, w0) : nothing,
        )
        if stack_rung == "e6_actuated"
            stack_controller = PE.LVLHCascadeAttitudeControlModel(
                q_cmd_lb=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
                k_out=SVector{3, Float64}(0.05, 0.05, 0.05),
                w_max=0.01,
                k_rate=SVector{3, Float64}(50.0, 50.0, 50.0),
                tau_max=0.01,
            )
            stack_torque = (t, r, v, q, w) -> PE._lvlh_cascade_torque(stack_controller, r, v, q, w)
            stack_l_pi = MMatrix{3, 3, Float64}(1.0I)
            stack_b_field = (t, r_ii) -> PE.get_magnetic_field_dipole(r_ii, stack_l_pi)
            stack_mgr = MagneticMomentumManagerModel[
                MagneticMomentumManagerModel(
                    sat_idx=i, mu_gain=7.0e-4, m_max_am2=1.8,
                    h_wheels_0=SVector{3, Float64}(0.0, 0.0, 0.0),
                    commanded_torque=stack_torque, b_field_ii=stack_b_field
                )
                for i in 1:stack_n
            ]
            return ppc_build_config(
                planet=planet,
                spacecraft=stack_spacecraft,
                mission_time_s=stack_time,
                orientation_sim=true,
                dynamic_effectors=(stack_effectors..., stack_controller),
                density_model=stack_density,
                control_effectors=Tuple(stack_mgr),
                control_rates=fill(1.0, stack_n),
                dt_max_orbit=10.0
            )
        end
        return ppc_build_config(
            planet=planet,
            spacecraft=stack_spacecraft,
            mission_time_s=stack_time,
            orientation_sim=stack_attitude,
            dynamic_effectors=stack_effectors,
            density_model=stack_density,
            dt_max_orbit=10.0
        )

    # ── B12: interacting vs. independent propagation, matched work ────────────
    #
    # Review point 8 asks for "interacting versus independent propagation" as an
    # axis. B6 gestured at it by putting constellation cases and Monte Carlo cases
    # in the same phase, but those differ in physics, duration and satellite count
    # as well as in coupling, so nothing there isolates coupling itself.
    #
    # These two builders are the matched pair. Identical planet, orbit geometry,
    # effector stack, density model, tolerances and mission duration; the only
    # difference is how the N satellite-hours are arranged:
    #
    #   interact_<N>sat_1hr   one simulation propagating N satellites -- shared
    #                         config, shared environment queries, one solver, one
    #                         coupled state vector. Outer parallelism is across
    #                         satellites inside one RHS call.
    #   independent_1sat_1hr  one satellite, run at mc_samples = N -- N fully
    #                         independent solves. Outer parallelism is across
    #                         samples (threads or processes), and the RHS is
    #                         single-satellite.
    #
    # Run at N = 16/64/256 on both sides, the pair is a direct measurement of the
    # manuscript's third contribution: whether the router correctly distinguishes
    # shared-state multi-spacecraft propagation from process-isolable independent
    # campaigns at matched total work.
    elseif occursin(r"^interact_[0-9]+sat_1hr$", case_name)
        coupling_n = parse(Int, match(r"^interact_([0-9]+)sat", case_name).captures[1])
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, coupling_n),
            mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=300.0, full=3600.0),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 50), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=20.0
        )
    elseif occursin(r"^mcgrid_[0-9]+sat_[0-9]+mc$", case_name)
        # Joint inner/outer routing. Every other Monte Carlo case in the catalog
        # carries exactly one spacecraft per sample, and every constellation case
        # carries exactly one sample, so in all of them one of the two axes is
        # empty and the route follows from the workload shape without a trade
        # being made. These rungs put real work on both axes at once: S samples of
        # N spacecraft each, where the same core budget can be spent fanning out
        # samples across processes, spreading spacecraft across threads inside one
        # sample, or splitting between the two.
        #
        # The rungs are sized to hold N*S constant within a grid, so the total
        # work is identical across a grid and only the aspect ratio changes. A router that
        # is genuinely choosing should answer differently at each rung -- outer for
        # many small samples, inner for few large ones -- and 32x32 is deliberately
        # the ambiguous middle. Identical answers across all three would say the
        # split is derived from the outer decision rather than selected, which is
        # what _campaign_route_plan's inner_thread_budget = fld(nthreads, workers)
        # suggests and what this case exists to test.
        #
        # Physics matches interact_<N>sat_1hr exactly (degree 50 harmonics,
        # exponential atmosphere, aerodynamics, 1 h) so the per-spacecraft cost is
        # the already-measured 0.072 s/spacecraft-hour. The 1024 spacecraft-hour
        # rungs land near 74 s serial; the 128 spacecraft-hour rungs near 9 s,
        # which is still 3x PPB_NOISE_FLOOR_SERIAL_S.
        mcgrid_n = parse(Int, match(r"^mcgrid_([0-9]+)sat", case_name).captures[1])
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, mcgrid_n),
            mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=300.0, full=3600.0),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 50), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=20.0
        )
    elseif case_name == "independent_1sat_1hr"
        # Independent side of the pair above. Deliberately NOT jittered the way
        # the montecarlo_* cases are: an initial-condition spread would give each
        # sample a different step-count and make per-sample cost vary, which would
        # show up as load imbalance and confound the coupling comparison. Every
        # sample here is the same problem, so total work is exactly N x the
        # single-satellite cost and the interacting/independent contrast is clean.
        return ppc_build_config(
            planet=planet,
            spacecraft=[ppc_spacecraft(planet; id=1)],
            mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=300.0, full=3600.0),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 50), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=20.0
        )

    # ── B14: mission duration and output cadence ──────────────────────────────
    #
    # Duration rungs: 15 min / 1 h / 6 h at N=1024, L50, vacuum. The 1 h and 6 h
    # points already exist as gravity_1024sat_l50_vacuum_1hr and
    # heavy_1024sat_l50_6hr; this adds the short end.
    elseif case_name == "gravity_1024sat_l50_vacuum_15min"
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 1024),
            mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=120.0, full=900.0),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 50),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=20.0
        )

    # Cadence rungs: identical physics and duration to
    # gravity_1024sat_l50_vacuum_1hr, varying only how much trajectory output the
    # run produces. See ppc_build_config's comment for why this goes through
    # `results`/`data_rate` and not `num_steps_to_save` (which is inert).
    #
    #   cadence_none   results=false -- no saving callback, no per-step solution
    #                  storage. This is what every other case in the catalog does,
    #                  so it is the ladder's zero point rather than a new rung.
    #   cadence_60s / _10s / _1s   results=true with saveat = 60 / 10 / 1 s.
    #
    # At 1024 satellites the saving callback's per-satellite _save_snapshot is a
    # genuinely serial spine, so this ladder measures an Amdahl term the router
    # has no visibility into at all -- which is the finding either way it comes
    # out.
    elseif startswith(case_name, "cadence_1024sat_")
        cadence_results, cadence_rate = if case_name == "cadence_1024sat_none"
            (false, 10.0)
        elseif case_name == "cadence_1024sat_10s"
            (true, 10.0)
        elseif case_name == "cadence_1024sat_1s"
            (true, 1.0)
        else
            throw(ArgumentError("Unknown output-cadence case '$case_name'."))
        end
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 1024),
            mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=300.0, full=3600.0),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 50),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=20.0,
            results=cadence_results,
            data_rate=cadence_rate
        )

    # ── Heavy scaling cases ────────────────────────────────────────────────────
    #
    # Every pre-existing constellation case in this catalog resolves in well
    # under a second of actual integration once JIT compilation is excluded
    # (measured post-warm-up on a 12-core reference box: 0.017 s for
    # gravity_16sat_l20_vacuum, 0.18 s for multi_64_high_fidelity, 2.9 s for
    # multi_256_high_fidelity). At that size the fixed per-solve costs -- problem
    # assembly, solver setup, the serial spine of step control and solution
    # saving -- dominate, so no thread count can move the wall time and every
    # routing profile measures the same thing. The cases below are sized so the
    # serial baseline is ~8-10 s and the parallelisable RHS is the majority of
    # it, which is the regime where outer/inner routing decisions are actually
    # observable. Reference measurements (12 physical cores, post-warm-up):
    #
    #   heavy_1024sat_l50_6hr        7.9 s serial -> 2.0 s @ 4 threads (3.9x)
    #   heavy_4096sat_l50_1hr        9.7 s serial -> 1.9 s @ 12 threads (5.1x)
    #
    # The 1024- vs 4096-satellite pair is the important contrast: 1024
    # saturates at ~4 workers, 4096 keeps scaling to 12+, so the two together
    # show *where* the per-RHS fork/join overhead stops being amortised rather
    # than just asserting a single speedup number.

    elseif case_name == "heavy_1024sat_l50_6hr"
        # RHS-bound vacuum constellation: one heavy effector (L50 harmonics, the
        # batched-SIMD path) and no callbacks, so essentially all of the solve
        # is the satellite-parallel region. This is the cleanest read on
        # outer-loop scaling in the catalog.
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 1024),
            mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=600.0, full=21600.0),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 50),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=20.0
        )
    elseif case_name == "heavy_4096sat_l50_1hr"
        # Same physics as above at 4x the constellation size and 1/6th the
        # mission length, so the serial baseline is comparable but each RHS call
        # carries 4x the work. Scaling continues past the point where the
        # 1024-satellite case flattens.
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 4096),
            mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=300.0, full=3600.0),
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 50),),
            density_model=NoAtmosphereModel(),
            dt_max_orbit=20.0
        )
    elseif case_name == "heavy_1024sat_fullstack_1hr"
        # Multi-effector counterpart: harmonics + SRP + aero over an analytic
        # atmosphere, so the flat effector queue has a heterogeneous effector mix
        # to schedule and the density/thermal callbacks carry real work. This is
        # the case where inner (callback) parallelism has something to do, as
        # opposed to the vacuum cases where it does not.
        #
        # One hour, not the six the vacuum cases use: the atmosphere path is
        # roughly 20x the cost of L50 harmonics at the same constellation size
        # (measured 163.6 s serial at 6 h vs 7.9 s for heavy_1024sat_l50_6hr), so
        # a matching duration would make this single case dominate the phase's
        # runtime. At 1 h the serial baseline is ~27 s -- still comfortably in
        # the regime where the RHS dominates fixed setup cost.
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 1024),
            mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=300.0, full=3600.0),
            orientation_sim=false,
            dynamic_effectors=(
                ppc_harmonics_model(planet, 20),
                SolarRadiationPressureModel(1.2, 12.0),
                AerodynamicCoefficientfM(),
            ),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=20.0
        )
    elseif case_name == "heavy_256sat_coupled6dof_2hr"
        # Coupled 6-DOF at constellation scale: attitude propagation on for every
        # satellite on top of harmonics/SRP/aero, which is the configuration the
        # manuscript's architecture claim is actually about and which no other
        # many-satellite case in this catalog exercises (they all run
        # orientation_sim=false). Measured 88.7 s serial -- the heaviest case in
        # this family, and deliberately kept at 2 h (~1.3 orbits) rather than
        # trimmed further, since a sub-orbit arc would not exercise the
        # attitude/aero coupling over a full periapsis passage.
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, 256; with_panel=true, panel_count=4, orientation_state=(q0, w0)),
            mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=300.0, full=7200.0),
            orientation_sim=true,
            dynamic_effectors=(
                ppc_harmonics_model(planet, 20),
                SolarRadiationPressureModel(1.2, 12.0),
                AerodynamicCoefficientfM(),
            ),
            density_model=ExponentialAtmosphereModel(planet),
            dt_max_orbit=10.0
        )
    elseif occursin(r"^heavy_[0-9]+sat_gram_nbody_l50$", case_name)
        # Native-library contention stressor: the only case in this catalog that
        # puts a *live* GRAM atmosphere, L50 spherical harmonics, and third-body
        # (Sun + Moon) gravity on a constellation at once.
        #
        # Purpose is to make the two serialised native libraries the bottleneck
        # rather than a rounding error. The existing GRAM cases do not: at 4 and
        # 16 satellites (multi_4_gram_live / multi_16_gram_live) the whole solve
        # is 1.2-5.4 s and neither thread count nor SPACEAGORA_GRAM_ISOLATED_POOL
        # moves it at all, because there simply are not enough density queries
        # in flight for GRAM's per-instance lock to be contended. This case
        # raises the query count (satellite count x steps) and adds a second
        # native serialisation point on top: NBodyGravityModel drives SPICE
        # ephemeris lookups under SPICE_LOCK, and CSPICE is not thread-safe.
        #
        # Both L50 harmonics and the two third bodies are there to make each RHS
        # call expensive in *thread-safe* work as well, so the measurement can
        # separate "threads do not help because the native calls serialise" from
        # "threads do not help because there is nothing to do".
        #
        # Satellite count is parameterised by the case name
        # (heavy_64sat_gram_nbody_l50, heavy_128sat_gram_nbody_l50, ...) so the
        # contention can be swept without editing this file. Start small: GRAM
        # constellations have a history of pathological step-size collapse and
        # native memory growth at higher N (see multi_16_gram_live's comment),
        # and while both appear fixed on this branch, that is one clean run's
        # worth of evidence, not a guarantee.
        n = parse(Int, match(r"^heavy_([0-9]+)sat", case_name).captures[1])
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, n),
            mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=90.0, full=600.0),
            orientation_sim=false,
            dynamic_effectors=(
                ppc_harmonics_model(planet, 50),
                NBodyGravityModel(body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet),
                AerodynamicCoefficientfM(),
            ),
            density_model=ppc_gram_atmosphere_model("earth"),
            dt_max_orbit=5.0
        )
    elseif case_name == "montecarlo_heavy_aerobraking"
        # Process-route stressor: same Mars aerobraking physics as
        # montecarlo_mars_aerobraking but a 12x longer arc -- ~1.05 s of
        # integration per sample (measured), against which pmap dispatch is a
        # ~1% term. With the original arc length the throughput-vs-workers curve
        # was dominated by per-worker process startup and JIT rather than by the
        # trials themselves.
        return ppc_build_config(
            planet=mars,
            spacecraft=[ppc_spacecraft(
                mars;
                ra_alt_m=4500e3 + randn(rng) * 100e3,
                rp_alt_m=max(110e3, 135e3 + randn(rng) * 10e3),
                i_deg=93.0,
                omega_deg=80.0,
                raan_deg=30.0,
                nu_deg=180.0 + randn(rng) * 4.0
            )],
            mission_time_s=ppc_mission_time(cfg.profile; test=10.0, smoke=600.0, full=21600.0),
            orientation_sim=false,
            # Inverse-square, not ppc_harmonics_model: that helper hardcodes the
            # Earth GGM05C coefficient file (PPC_EARTH_HARMONICS_FILE) and would
            # silently apply Earth's harmonics to Mars. Same choice as
            # montecarlo_mars_aerobraking.
            dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
            density_model=ExponentialAtmosphereModel(mars),
            dt_max_orbit=1.0
        )
    end

    throw(ArgumentError("Unknown parallelization performance case '$case_name'."))
end

function ppc_case_catalog()::Dict{String, PPCCaseSpec}
    cases = Dict{String, PPCCaseSpec}()
    add!(name, family, description; montecarlo=false, orientation=false, default_samples=1) = begin
        cases[name] = PPCCaseSpec(
            name=name,
            family=family,
            description=description,
            builder=ppc_single_config,
            montecarlo=montecarlo,
            orientation=orientation,
            default_samples=default_samples
        )
    end
    add!("single_inverse_square_vacuum", "gravity_only", "1 spacecraft, inverse-square gravity, no atmosphere")
    for n in (4, 16, 64, 256, 1024, 2048)
        add!("gravity_$(n)sat_inverse_square_vacuum", "gravity_only", "$(n) spacecraft, inverse-square gravity, no atmosphere")
        add!("gravity_$(n)sat_l20_vacuum", "gravity_only", "$(n) spacecraft, L20 harmonics, no atmosphere")
    end
    # Ad hoc spacecraft-count x thread-count sweep, L50 harmonics, 1hr mission.
    for n in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096)
        add!("gravity_$(n)sat_l50_vacuum_1hr", "gravity_only", "$(n) spacecraft, L50 harmonics, no atmosphere, 1hr mission")
    end
    add!("single_harmonics_l20_vacuum", "gravity_only", "1 spacecraft, L20 harmonics, no atmosphere")
    add!("single_harmonics_l50_vacuum", "few_sat_high_fidelity", "1 spacecraft, L50 harmonics, no atmosphere")
    add!("srp_heavy_high_area", "few_sat_high_fidelity", "1 high-area articulated spacecraft with stacked SRP", orientation=true)
    add!("articulated_1sat_fullstack", "few_sat_high_fidelity", "1 articulated spacecraft with harmonics, aero, thermal, and attitude", orientation=true)
    add!("multi_16_aero_surrogate_cached", "many_sat_high_fidelity", "16 spacecraft with aero and analytic density")
    add!("multi_64_high_fidelity", "many_sat_high_fidelity", "64 spacecraft with harmonics, SRP, aero, and analytic density")
    add!("multi_128_high_fidelity", "many_sat_high_fidelity", "128 spacecraft with harmonics, SRP, aero, and analytic density")
    add!("callback_128_aero_thermal", "many_sat_high_fidelity", "128 spacecraft callback-heavy aero and thermal stress")
    add!("montecarlo_high_accuracy", "monte_carlo", "Monte Carlo high-accuracy gravity seeds", montecarlo=true)
    add!("montecarlo_multi_sat", "monte_carlo", "Monte Carlo 4-spacecraft gravity seeds", montecarlo=true)
    add!("montecarlo_mars_aerobraking", "monte_carlo", "Monte Carlo Mars aerobraking seeds", montecarlo=true)

    # Router-evaluation axis coverage (B6 / point 8).
    add!("multi_4_aero_surrogate_cached", "many_sat_high_fidelity", "4 spacecraft with aero and analytic density")
    add!("multi_256_high_fidelity", "many_sat_high_fidelity", "256 spacecraft with harmonics, SRP, aero, and analytic density")
    add!("multi_16_gram_live", "gram_live", "16 spacecraft, harmonics and aero, native GRAM atmosphere (interacting)")
    add!("multi_4_gram_live", "gram_live", "4 spacecraft, harmonics and aero, native GRAM atmosphere (interacting, small-N)")
    add!("campaign_gnc_6dof_gram", "full_stack_campaign",
        "4 spacecraft, LVLH attitude control + magnetorquers, coupled 6-DOF, live Earth GRAM, run as N independent samples",
        montecarlo=true, orientation=true)
    add!("montecarlo_mars_gram_live", "gram_live", "Monte Carlo Mars aerobraking, native GRAM atmosphere (independent)", montecarlo=true)
    add!("multi_8sat_magnetorquer_attitude", "actuator", "8 spacecraft, LVLH attitude control and magnetorquer unloading actuator", orientation=true)
    add!("gravity_16sat_l20_vacuum_longmission", "duration_cadence", "16 spacecraft, L20 harmonics, ~4x mission duration")

    # Heavy scaling cases (see the ppc_single_config branches for sizing rationale).
    add!("heavy_1024sat_l50_6hr", "heavy_scaling", "1024 spacecraft, L50 harmonics, 6hr mission, RHS-bound vacuum")
    add!("heavy_4096sat_l50_1hr", "heavy_scaling", "4096 spacecraft, L50 harmonics, 1hr mission, RHS-bound vacuum")
    add!("heavy_1024sat_fullstack_1hr", "heavy_scaling", "1024 spacecraft, harmonics + SRP + aero + analytic density, 1hr mission")
    add!("heavy_256sat_coupled6dof_2hr", "heavy_scaling", "256 spacecraft, coupled 6-DOF with harmonics + SRP + aero, 2hr mission", orientation=true)
    add!("montecarlo_heavy_aerobraking", "heavy_scaling", "Monte Carlo Mars aerobraking, 12x longer arc per sample", montecarlo=true)
    for n in (16, 32, 64, 128, 256)
        add!("heavy_$(n)sat_gram_nbody_l50", "heavy_scaling",
             "$(n) spacecraft, live GRAM atmosphere + L50 harmonics + Sun/Moon third-body gravity (native-library contention)")
    end

    # ── Expanded router evaluation (B9-B14, review point 8) ───────────────────
    # Single-variable ladders, one per axis point 8 names. See the corresponding
    # ppc_single_config branches for the sizing and isolation rationale.

    # Atmosphere / GRAM fidelity (B10): N and duration fixed, density path varies.
    add!("atmo256_vacuum_10min", "atmosphere_ladder", "256 spacecraft, L50 harmonics, no atmosphere (fidelity-ladder floor)")
    add!("atmo256_exponential_10min", "atmosphere_ladder", "256 spacecraft, L50 harmonics + aero, analytic exponential density")
    add!("atmo256_gram_surrogate_10min", "atmosphere_ladder", "256 spacecraft, L50 harmonics + aero, precomputed GRAM surrogate table")
    add!("atmo256_gram_live_10min", "atmosphere_ladder", "256 spacecraft, L50 harmonics + aero, live native GRAM atmosphere")
    add!("atmo256_gram_live_nbody_10min", "atmosphere_ladder", "256 spacecraft, live GRAM + Sun/Moon third-body (two serialised native libraries)")

    # Force / actuator model count (B11): N and duration fixed, model count varies.
    # Registered at several N because the actuator rung is quadratic in satellite
    # count (see the ppc_single_config branch), so the ladder has to be run at
    # whichever N keeps every rung both measurable and affordable; B11 picks one.
    for n in (32, 64, 256)
        add!("stack$(n)_e1_harm", "effector_ladder", "$(n) spacecraft, L20 harmonics only (1 dynamic effector)")
        add!("stack$(n)_e2_srp", "effector_ladder", "$(n) spacecraft, harmonics + SRP (2 dynamic effectors)")
        add!("stack$(n)_e3_aero", "effector_ladder", "$(n) spacecraft, harmonics + SRP + aero (3 effectors, density callback)")
        add!("stack$(n)_e4_nbody", "effector_ladder", "$(n) spacecraft, harmonics + SRP + third-body + aero (4 effectors)")
        add!("stack$(n)_e5_6dof", "effector_ladder", "$(n) spacecraft, same 4 effectors with attitude propagation on", orientation=true)
        add!("stack$(n)_e6_actuated", "effector_ladder", "$(n) spacecraft, 6-DOF plus LVLH controller and per-satellite magnetorquer actuators", orientation=true)
    end

    # Interacting vs. independent propagation at matched work (B12).
    for n in (16, 64, 256)
        add!("interact_$(n)sat_1hr", "coupling_mode", "$(n) spacecraft propagated as one coupled simulation (interacting side)")
    end
    add!("independent_1sat_1hr", "coupling_mode", "1 spacecraft, identical physics, run as N independent samples (independent side)", montecarlo=true)

    # Joint inner/outer routing grids. N x S is held constant within each grid so
    # its rungs are the same total work at different aspect ratios. The sample
    # count is carried in the name and in default_samples, and execution honours
    # it for this family: a joint_routing case runs at its own S unless the caller
    # names one with --mc-samples. Every other family still takes the sample
    # ladder from PPCConfig.mc_samples, so nothing outside this grid is affected.
    # See the mcgrid_ builder for why this axis is missing from every other case.
    #
    # Two totals. The 1024 spacecraft-hour grid is the reference and the most
    # expensive family in the catalog at ~74 s serial per rung, which is why it is
    # still unrun. The 128 spacecraft-hour grid repeats the same aspect-ratio
    # sweep at ~9 s serial, an eighth of the cost and still well clear of the 3 s
    # measurability floor, so the joint-routing axis can be exercised in a normal
    # working session rather than an overnight slot. Every rung of both grids
    # carries multiple spacecraft per sample, which is the property that makes the
    # inner axis non-empty.
    for (mcgrid_n, mcgrid_s) in ((64, 16), (32, 32), (16, 64),   # 1024 sc-h
                                 (32,  4), (16,  8), ( 8, 16))   #  128 sc-h
        add!(
            "mcgrid_$(mcgrid_n)sat_$(mcgrid_s)mc", "joint_routing",
            "$(mcgrid_n) spacecraft x $(mcgrid_s) samples, " *
            "$(mcgrid_n * mcgrid_s) spacecraft-hours, both routing axes live",
            montecarlo=true, default_samples=mcgrid_s
        )
    end

    # Mission duration and output cadence (B14).
    add!("gravity_1024sat_l50_vacuum_15min", "duration_cadence", "1024 spacecraft, L50 harmonics, 15 min mission (short end of duration ladder)")
    # MEASURED, and the reason this ladder is three rungs and not six.
    #
    # At N=1024 over a 10 s mission, post-warm-up and serial: 0.051 s with output
    # off, and 3.13 / 3.09 / 3.17 s at saveat = 60 / 10 / 1 s. Turning trajectory
    # output on costs 62x; the requested cadence then changes that by under 1%,
    # i.e. it is flat. Pinning SPACEAGORA_SOLVER_SAVE_EVERYSTEP=0 and
    # SPACEAGORA_SOLVER_SAVE_ON=0 does not move it either (3.25 / 3.25 / 3.24 s),
    # so the cost is neither per-step solution storage nor the solver's own save
    # path -- it is a fixed cost of having the saving callback in the CallbackSet
    # at all, and it does not scale with how often the callback actually saves.
    #
    # So the honest form of point 8's "output cadence" axis is on-vs-off, not a
    # cadence sweep. The 1 s rung is retained as the evidence that cadence really
    # is flat once the full mode and thread ladder runs over it; adding 60 s and a
    # per-step variant on top would spend three times the budget re-measuring the
    # same number, which is the mistake B6 made.
    add!("cadence_1024sat_none", "duration_cadence", "1024 spacecraft, L50 harmonics, no trajectory output (cadence-ladder zero point)")
    add!("cadence_1024sat_10s", "duration_cadence", "1024 spacecraft, L50 harmonics, trajectory saved every 10 s")
    add!("cadence_1024sat_1s", "duration_cadence", "1024 spacecraft, L50 harmonics, trajectory saved every 1 s (cadence-flatness check)")

    return cases
end

function ppc_resolve_cases(requested::Vector{String})::Vector{String}
    catalog = ppc_case_catalog()
    names = isempty(requested) ? sort(collect(keys(catalog))) : requested
    unknown = [name for name in names if !haskey(catalog, name)]
    isempty(unknown) || throw(ArgumentError("Unknown case(s): $(join(unknown, ", "))"))
    return names
end
