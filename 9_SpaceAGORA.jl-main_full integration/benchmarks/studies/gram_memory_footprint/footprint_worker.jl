# One measured subprocess for the GRAM memory-footprint study. One invocation =
# one (density path, constellation size) point: it loads SpaceAGORA (and, for a
# native path, GRAMSuite), reads its resident set BEFORE the workload exists,
# builds an Earth LEO constellation, solves one short mission, and prints a
# single machine-parseable line:
#
#   GMF_RESULT ok=<bool> rss_base_mb=... rss_built_mb=... rss_solved_mb=... \
#              rss_retained_mb=... maxrss_mb=... vmrss_base_mb=... \
#              vmrss_retained_mb=... solve_s=...
#
# `rss_base_mb` is the package + SPICE + native-GRAM image with no spacecraft in
# it; `rss_retained_mb - rss_base_mb` is what N spacecraft on this path actually
# cost, which is exactly the term `native_gram_worker_extra_bytes` charges.
# Both ends of that subtraction are read after a full collection, so the figure
# is live state rather than whatever garbage the last solve happened to leave;
# `maxrss_mb - rss_base_mb` is the same quantity taken at the high-water mark,
# and the routing constant is fitted against that conservative form.
#
# Configuration arrives entirely through the environment (set by run_footprint.jl):
#
#   GMF_PATH       none | gram_point | gram_freeze_step | gram_lookahead | gram_surrogate
#   GMF_N_SATS     constellation size
#   GMF_MISSION_S  mission length [s] (default 300)
#   GMF_ALT_KM     circular altitude [km] (default 150 for GRAM paths, 550 for none)
#   GMF_GRAVITY    invsq | l20 | l50 (default l20)
#   GMF_SOLVE_REPEATS  how many solves of the SAME size to run before reading
#                  the peak (default 1). The paper_scenarios workers run four
#                  (one warmup + three timed), and a process-pool worker runs
#                  many samples in its lifetime, so a repeat count above one
#                  measures the heap the allocator settles at rather than the
#                  first solve's peak.
#   GMF_ISOLATE    1 = run_simulation's own `isolate_state` deep copy (its
#                  production default, and what a campaign member pays);
#                  0 = no copy, which is what the paper_scenarios S1/S2 rows
#                  were measured with. Default 1, the conservative arm.
#
# The SPACEAGORA_* switches that select the native path (freeze-per-step, the
# vacuum-predicted look-ahead cache) are set by the controller as real process
# environment, so they are already in place when the callbacks read them.
#
# Run with --threads=1: a process-pool worker -- the thing the routing estimate
# prices -- is started `--threads=1`, so that is the configuration whose
# footprint the constant must describe.

const GMF_REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const GMF_PATH = get(ENV, "GMF_PATH", "gram_point")
const GMF_N_SATS = parse(Int, get(ENV, "GMF_N_SATS", "16"))
const GMF_MISSION_S = parse(Float64, get(ENV, "GMF_MISSION_S", "300.0"))
const GMF_GRAVITY = get(ENV, "GMF_GRAVITY", "l20")
const GMF_ALT_KM = parse(Float64, get(ENV, "GMF_ALT_KM", GMF_PATH == "none" ? "550.0" : "150.0"))
const GMF_ISOLATE = get(ENV, "GMF_ISOLATE", "1") == "1"
const GMF_SOLVE_REPEATS = max(1, parse(Int, get(ENV, "GMF_SOLVE_REPEATS", "1")))

const GMF_NATIVE_PATHS = ("gram_point", "gram_freeze_step", "gram_lookahead")
const GMF_NEEDS_GRAM = GMF_PATH in GMF_NATIVE_PATHS || GMF_PATH == "gram_surrogate"

include(joinpath(GMF_REPO_ROOT, "examples", "common.jl"))

"""Resident set of this process in bytes, read from `/proc/self/statm` (Linux)
or `Sys.maxrss()` elsewhere. This is the live resident set, not the peak."""
function gmf_rss_bytes()::Int
    if Sys.islinux() && isfile("/proc/self/statm")
        fields = split(read("/proc/self/statm", String))
        pages = length(fields) >= 2 ? tryparse(Int, fields[2]) : nothing
        pages === nothing || return pages * 4096
    end
    return Int(Sys.maxrss())
end

"""`VmRSS` from `/proc/self/status` in bytes, or `-1` where it is unreadable.
Recorded alongside `gmf_rss_bytes` as an independent reading of the same
quantity, so a disagreement between the two is visible in the CSV."""
function gmf_vmrss_bytes()::Int
    Sys.islinux() && isfile("/proc/self/status") || return -1
    for line in eachline("/proc/self/status")
        startswith(line, "VmRSS:") || continue
        kb = tryparse(Int, first(split(strip(split(line, ":"; limit=2)[2]))))
        kb === nothing || return kb * 1024
    end
    return -1
end

gmf_mb(bytes::Integer)::Float64 = bytes < 0 ? -1.0 : bytes / (1 << 20)

const GMF_SPICE_AVAILABLE = isfile(joinpath(SPICE_PATH, "pck", "pck00011.tpc"))
if GMF_NEEDS_GRAM && !GMF_SPICE_AVAILABLE
    error("GMF_PATH=$(GMF_PATH) needs the SPICE kernels under $(SPICE_PATH); none found.")
end
GMF_NEEDS_GRAM && ensure_gramsuite_loaded!()

gmf_earth() = GMF_SPICE_AVAILABLE ? Earth("", SPICE_PATH) : Earth()

# The run's epoch, shared by the configuration and by any density model built
# for it. A GRAMAtmosphereModelSurrogate carries a fixed table and refuses to be
# re-epoched at solve time, so its native fallback has to be constructed at this
# epoch up front; the other models accept the engine's realignment.
gmf_initial_time() = InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)

function gmf_density_model()
    GMF_PATH == "none" && return NoAtmosphereModel()
    if GMF_PATH == "gram_surrogate"
        return SimulationModel.GRAMAtmosphereModelSurrogate(
            planet_name="earth", initial_time=gmf_initial_time()
        )
    end
    # gram_point, gram_freeze_step and gram_lookahead all use the same native
    # model object; the calling method is selected purely through SPACEAGORA_*.
    return SimulationModel.GRAMAtmosphereModel(planet_name="earth")
end

function gmf_gravity_effector(planet)
    GMF_GRAVITY == "invsq" && return InverseSquaredGravityModel()
    lm = GMF_GRAVITY == "l50" ? 50 : 20
    return GravitationalHarmonicsModel(
        lm, lm, joinpath(GMF_REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv"), planet
    )
end

# Same Earth LEO constellation geometry as benchmarks/studies/paper_scenarios's
# scenario_worker.jl, so a footprint here is comparable with an S1/S2 row there.
function gmf_build_config(n_sats::Int)
    planet = gmf_earth()
    gravity = gmf_gravity_effector(planet)
    effectors = GMF_PATH == "none" ? (gravity,) : (gravity, AerodynamicCoefficientfM())
    alt_m = GMF_ALT_KM * 1e3

    spacecraft = SpacecraftModel[]
    for i in 1:n_sats
        root = Link(root=true, m=500.0, ref_area=12.0)
        phase = 50.0 * (i - 1) / max(n_sats, 1)
        ic = InitialCondition(
            ra=planet.Rp_e + alt_m + phase,
            rp=planet.Rp_e + alt_m + phase,
            i=53.0,
            ω=0.0,
            Ω=10.0,
            ν=360.0 * (i - 1) / max(n_sats, 1)
        )
        push!(spacecraft, SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, i))
    end

    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=GMF_MISSION_S,
            orientation_sim=false,
            num_steps_to_save=20
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=300.0,
            density_model=gmf_density_model(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false,
            ephemerides_model=GMF_SPICE_AVAILABLE ? SpiceEphemeridesModel() : SimpleEphemeridesModel()
        ),
        dynamics_model=DynamicsModel(spacecraft, effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=gmf_initial_time(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=2.0)
    )
end

function main()
    # A one-spacecraft solve first: it pays the solver's compilation and the
    # native GRAM image's first-call allocations, neither of which is a
    # per-spacecraft cost, so the baseline below already contains them and the
    # measured slope is the workload alone.
    warm = gmf_build_config(1)
    SpaceAGORA.run_simulation(warm; isolate_state=GMF_ISOLATE, return_solver_metadata=true)
    warm = nothing
    GC.gc(); GC.gc()

    rss_base = gmf_rss_bytes()
    vmrss_base = gmf_vmrss_bytes()

    args = gmf_build_config(GMF_N_SATS)
    GC.gc()
    rss_built = gmf_rss_bytes()

    local result
    solve_s = 0.0
    ok = true
    for _ in 1:GMF_SOLVE_REPEATS
        solve_s = @elapsed result = SpaceAGORA.run_simulation(
            args; isolate_state=GMF_ISOLATE, return_solver_metadata=true
        )
        ok &= String(result.retcode) == "Success"
    end
    retcode = String(result.retcode)

    rss_solved = gmf_rss_bytes()
    # Retained: the same reading after a full collection, so it is live state
    # and not the garbage of the last step. The baseline above is taken the same
    # way, so the two are subtractable.
    result = nothing
    GC.gc(); GC.gc()
    rss_retained = gmf_rss_bytes()
    vmrss_retained = gmf_vmrss_bytes()
    maxrss = Int(Sys.maxrss())

    println("GMF_RESULT ok=$(ok) " *
            "rss_base_mb=$(round(gmf_mb(rss_base); digits=1)) " *
            "rss_built_mb=$(round(gmf_mb(rss_built); digits=1)) " *
            "rss_solved_mb=$(round(gmf_mb(rss_solved); digits=1)) " *
            "rss_retained_mb=$(round(gmf_mb(rss_retained); digits=1)) " *
            "maxrss_mb=$(round(gmf_mb(maxrss); digits=1)) " *
            "vmrss_base_mb=$(round(gmf_mb(vmrss_base); digits=1)) " *
            "vmrss_retained_mb=$(round(gmf_mb(vmrss_retained); digits=1)) " *
            "solve_s=$(round(solve_s; digits=3)) " *
            "isolate=$(GMF_ISOLATE ? 1 : 0) " *
            "solve_repeats=$(GMF_SOLVE_REPEATS) " *
            "retcode=$(retcode)")
    flush(stdout)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
