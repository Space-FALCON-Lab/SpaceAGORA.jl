# Does every member of the P6 native-GRAM traces actually call native GRAM, and is
# every one of them inside the atmosphere so the look-ahead cache builds?
#
# Audits the two harness cases exactly as the harness builds them
# (`ppc_single_config` in parallelization_performance/cases.jl):
#
#   aero_<N>sat_l50_gram_lookahead_100s   P6 trace 5, one N-member constellation
#   aero_<N>sat_l50_gram_process_100s     P6 trace 6, N one-spacecraft samples
#                                         (sample k built with mc_index = k)
#
# and, for contrast, the definition both used before (`ppc_constellation`,
# entry interface 120 km), rebuilt here from the same helpers.
#
# Part A, per member, at t = 0:
#   alt_km          norm(position) - planet.Rp_e, the expression the engine's own
#                   `_initialize_in_atmosphere_flags!` uses
#   below_2000km    whether that altitude is under the cut-off above which
#                   `getDensity(::GRAMAtmosphereModel, ...)` returns vacuum
#   in_atmosphere   the engine's own flag, set by calling
#                   `_initialize_in_atmosphere_flags!` on the harness's config;
#                   it is what gates the vacuum-predicted look-ahead cache
#   native_calls    native GRAM calls made by ONE production density evaluation
#                   for that member, counted as `gram_density` acquisitions of
#                   the native-lock counters around `getDensity`. The GRAM gate
#                   in that function depends on altitude alone, so it is
#                   evaluated at the member's altitude with latitude, longitude
#                   and elapsed time fixed at 0; 1 means the member calls GRAM
#                   on every density callback, 0 that it never does.
#
# Part B, end to end, trace 5 at the largest N: a short real solve through the
# harness config under each density path, reporting total `gram_density`
# acquisitions. On the freeze-per-step path every callback evaluates every
# member once, so acquisitions / N must come out an integer equal to the number
# of density callbacks if and only if every member called GRAM every time.
#
# Usage (one process, memory-capped):
#   systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q -- \
#     julia --project=. --threads=1 benchmarks/studies/gram_thread_scaling/p6_case_audit.jl \
#       [--sizes=256,4096] [--e2e-mission=10] [--out=results/p6_case_audit.csv]

using Printf

const GTS_DIR = @__DIR__
const GTS_REPO_ROOT = normpath(joinpath(GTS_DIR, "..", "..", ".."))
const GTS_PPC_DIR = joinpath(GTS_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")
include(joinpath(GTS_PPC_DIR, "cli.jl"))
include(joinpath(GTS_PPC_DIR, "modes.jl"))
include(joinpath(GTS_PPC_DIR, "cases.jl"))
ppc_ensure_gramsuite_loaded!()

const SE = SpaceAGORA.SimulationEngine
const EM = SpaceAGORA.SimulationModel.EnvironmentModels
const RS = SpaceAGORA.RuntimeServices

audit_arg(key, default) = begin
    for a in ARGS
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end
const AUDIT_SIZES = [parse(Int, x) for x in split(audit_arg("sizes", "256,4096"), ",")]
const AUDIT_E2E_MISSION = parse(Float64, audit_arg("e2e-mission", "10"))
const AUDIT_OUT = audit_arg("out", joinpath(GTS_DIR, "results", "p6_case_audit.csv"))

# profile="full" so the case builders use the mission length in the case name.
const AUDIT_CFG = PPCConfig(profile="full")

function previous_definition(n::Int; sample::Union{Nothing, Int}=nothing)
    planet = Earth("", PPC_SPICE_PATH)
    sc = sample === nothing ? ppc_constellation(planet, n) : [ppc_spacecraft(planet; id=1)]
    return ppc_build_config(
        planet=planet, spacecraft=sc, mission_time_s=100.0, orientation_sim=false,
        dynamic_effectors=(ppc_harmonics_model(planet, 50), AerodynamicCoefficientfM()),
        density_model=ppc_gram_atmosphere_model("earth"), dt_max_orbit=5.0
    )
end

"""Per-member (alt_km, in_atmosphere) for every spacecraft in one config."""
function member_states(args)
    ic = SE.build_initial_conditions(args)
    p = ODEParams(n_sats=length(args.dynamics_model.spacecraft), args=args)
    SE._initialize_in_atmosphere_flags!(p, ic)
    Rp = args.environment_model.planet.Rp_e
    out = Tuple{Float64, Bool}[]
    for i in 1:length(args.dynamics_model.spacecraft)
        alt = norm(SE._state_position_ii(ic, i)) - Rp
        push!(out, (alt, p.shared_buffers.in_atmosphere[i]))
    end
    return out, p
end

"""Native GRAM calls made by one production density evaluation at altitude h."""
function native_calls_at(h::Float64, args, p)
    RS.reset_native_lock_stats!()
    EM.getDensity(args.environment_model.density_model, h, 0.0, 0.0, 0.0, true, p)
    return RS.native_lock_stats_snapshot().sites.gram_density.acquisitions
end

const ROWS = String[]
const SUMMARY = String[]

function audit!(definition::String, case::String, n::Int, members)
    # members: vector of (member index, args, p, alt, in_atm)
    calls = 0; below = 0; inatm = 0
    for (m, args, p, alt, in_atm) in members
        nc = native_calls_at(alt, args, p)
        b = alt < 2000.0e3
        calls += nc > 0; below += b; inatm += in_atm
        push!(ROWS, join((definition, case, n, m, @sprintf("%.3f", alt * 1e-3),
                          b, in_atm, nc), ","))
    end
    line = @sprintf("%-9s %-38s N=%-5d below_2000km=%5d  in_atmosphere=%5d  calling_GRAM=%5d",
                    definition, case, n, below, inatm, calls)
    println(line); flush(stdout)
    push!(SUMMARY, line)
end

function main()
    # GRAM's injected ephemeris needs the SPICE kernels the engine loads; a bare
    # getDensity before any solve fails on the leapseconds kernel. One tiny solve
    # loads them for the process.
    SE.run_simulation(ppc_single_config("aero_16sat_l50_gram_lookahead_100s",
                                        PPCConfig(profile="test")); isolate_state=false)

    for n in AUDIT_SIZES
        c5 = "aero_$(n)sat_l50_gram_lookahead_100s"
        c6 = "aero_$(n)sat_l50_gram_process_100s"

        # Trace 5, current and previous.
        for (def, args) in (("current", ppc_single_config(c5, AUDIT_CFG)),
                            ("previous", previous_definition(n)))
            st, p = member_states(args)
            audit!(def, c5, n, [(i, args, p, st[i][1], st[i][2]) for i in eachindex(st)])
        end

        # Trace 6: one config per sample, exactly as the harness builds it.
        for def in ("current", "previous")
            members = Any[]
            for k in 1:n
                args = def == "current" ? ppc_single_config(c6, AUDIT_CFG; mc_index=k) :
                                          previous_definition(n; sample=k)
                st, p = member_states(args)
                push!(members, (k, args, p, st[1][1], st[1][2]))
            end
            audit!(def, c6, n, members)
        end
    end

    # Part B: end to end on the largest trace-5 case.
    n = maximum(AUDIT_SIZES)
    c5 = "aero_$(n)sat_l50_gram_lookahead_100s"
    println("\n-- end to end: $c5 at a $(AUDIT_E2E_MISSION) s mission, 1 thread --")
    base = ppc_single_config(c5, AUDIT_CFG)
    for (path, env) in (
        ("freeze", ["SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "1", "SPACEAGORA_VACUUM_GRAM_CACHE" => "0"]),
        ("lookahead", ["SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "0", "SPACEAGORA_VACUUM_GRAM_CACHE" => "1",
                       "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "20",
                       "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => string(AUDIT_E2E_MISSION + 500.0),
                       "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "1e8"]),
    )
        args = SimulationConfiguration(
            simulation_settings=base.simulation_settings,
            mission_configuration=MissionConfiguration(
                mission_type=base.mission_configuration.mission_type,
                keplerian=base.mission_configuration.keplerian,
                number_of_orbits=base.mission_configuration.number_of_orbits,
                mission_time=AUDIT_E2E_MISSION,
                orientation_sim=base.mission_configuration.orientation_sim,
                num_steps_to_save=base.mission_configuration.num_steps_to_save,
                data_rate=base.mission_configuration.data_rate),
            environment_model=base.environment_model, dynamics_model=base.dynamics_model,
            guidance_model=base.guidance_model, navigation_model=base.navigation_model,
            control_model=base.control_model, initial_time=base.initial_time,
            integration_tolerances=base.integration_tolerances)
        RS.reset_native_lock_stats!()
        withenv(env..., "SPACEAGORA_GRAM_ISOLATED_POOL" => "off") do
            SE.run_simulation(args; isolate_state=false)
        end
        acq = RS.native_lock_stats_snapshot().sites.gram_density.acquisitions
        line = @sprintf("e2e %-9s N=%d gram_density acquisitions=%d  per member=%.3f",
                        path, n, acq, acq / n)
        println(line); push!(SUMMARY, line)
    end

    mkpath(dirname(AUDIT_OUT))
    open(AUDIT_OUT, "w") do io
        println(io, "definition,case,n,member,alt_km,below_2000km,in_atmosphere,native_calls")
        foreach(r -> println(io, r), ROWS)
    end
    open(replace(AUDIT_OUT, ".csv" => "_summary.txt"), "w") do io
        foreach(l -> println(io, l), SUMMARY)
    end
    println("\nwrote $(length(ROWS)) member rows to $(AUDIT_OUT)")
end

main()
