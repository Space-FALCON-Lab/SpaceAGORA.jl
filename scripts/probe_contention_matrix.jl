#!/usr/bin/env julia
#
# P1 gate data: what the contention probe reports across the workload families.
#
#   julia --project=. --threads=1 scripts/probe_contention_matrix.jl
#   PCM_N=1024 julia --project=. scripts/probe_contention_matrix.jl
#
# The allocator's scaling law needs two workload-side inputs it cannot count:
# allocation pressure and native-lock occupancy. This measures both on one warm
# RHS pass per family and prints them side by side, so the question the gate
# asks -- do these separate the families that scale from the families that
# invert -- is answered by a table rather than by argument.
#
# Two separations are being tested, and they are different claims:
#
#   ALLOCATION separates plateau from inversion. `gravity` iterates a
#   homogeneous effector 1-tuple and should allocate nothing per pass;
#   `fullstack` iterates a heterogeneous tuple whose element type is abstract
#   and should allocate on every satellite-effector pair. Those are the two
#   ends of the measured scaling behaviour (~4x plateau against inversion below
#   serial at twelve threads), and if the probe cannot tell them apart then the
#   contention coefficient is not predictable from it.
#
#   OCCUPANCY separates workloads the density-model family mislabels.
#   `spice_nbody` runs in vacuum -- it signs as `dens=none` -- and holds the
#   shared native lock on every third-body sample. `gram` names GRAM and, if it
#   is sampling a cached or out-of-atmosphere path, may hold it almost never.
#   A signature keyed on the density model gets both backwards; measured rho
#   does not.
#
# Single-threaded by default and deliberately: the probe reads process-wide
# allocation counters, so a concurrent worker would be attributed to the case
# under measurement.

using SpaceAGORA
using StaticArrays
using Printf

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const PC = SM.ParallelCost
const RS = SpaceAGORA.RuntimeServices

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const HARM_FILE = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
const N_SATS = parse(Int, get(ENV, "PCM_N", "256"))
const HARM_DEGREE = parse(Int, get(ENV, "PCM_L", "20"))
# Epoch step between probe passes. Non-zero so that epoch-memoised samples
# (third-body ephemeris above all) miss on every pass, as they do in a solve.
const MISSION_S = parse(Float64, get(ENV, "PCM_MISSION_S", "3600.0"))
# Spread the probe passes across the mission, not past its end: the ephemeris
# look-ahead caches are built over the mission span, and epochs beyond it miss
# every cache and report SPICE traffic a solve never pays.
const T_STEP = parse(Float64, get(ENV, "PCM_T_STEP", string(MISSION_S / 8)))

function ensure_gramsuite_loaded!()
    isdefined(SpaceAGORA, :GRAMSuite) && return true
    vendored = joinpath(REPO_ROOT, "data", "GRAMSuite.jl")
    try
        if Base.find_package("GRAMSuite") === nothing && isdir(vendored)
            pushfirst!(LOAD_PATH, vendored)
        end
        @eval import GRAMSuite
        return true
    catch err
        @info "GRAMSuite not loadable; the GRAM family will be skipped." exception = err
        return false
    end
end

include(joinpath(@__DIR__, "probe_contention_matrix_support.jl"))

# Case construction lives in the support file so the matrix and the thread
# ladder measure the same workload. A private copy here drifted once already:
# it omitted the ephemeris look-ahead caches that execution.jl populates, and
# the resulting SPICE traffic was reported as a property of the workload.
const build_case = build_probe_case

# Each family is built lazily and independently: an asset that is missing on
# this machine must cost its own row, not the table.
const GRAM_AVAILABLE = Ref(false)

function families()
    out = Pair{String, Function}[]
    earth = SM.Earth()
    harm() = SM.GravitationalHarmonicsModel(HARM_DEGREE, HARM_DEGREE, HARM_FILE, earth)

    push!(out, "gravity" => () -> build_case(
        earth, (harm(),), SM.NoAtmosphereModel(), SM.SimpleEphemeridesModel();
        ra_alt_m = 400e3, rp_alt_m = 380e3))

    push!(out, "aero" => () -> build_case(
        earth, (harm(), SM.AerodynamicCoefficientfM()),
        SM.ExponentialAtmosphereModel(earth), SM.SimpleEphemeridesModel();
        ra_alt_m = 200e3, rp_alt_m = 160e3))

    # SRP samples the Sun through SPICE, so this family needs a kernel-backed
    # planet and ephemeris. Built on SimpleEphemeridesModel it fails outright
    # ("At least one SPK file needs to be loaded by SPKLEF"), which is a
    # configuration error rather than a measurement.
    push!(out, "fullstack" => () -> begin
        spice_earth = SM.Earth(HARM_FILE)
        build_case(spice_earth,
                   (SM.GravitationalHarmonicsModel(HARM_DEGREE, HARM_DEGREE, HARM_FILE, spice_earth),
                    SM.SolarRadiationPressureModel(1.2, 12.0), SM.AerodynamicCoefficientfM()),
                   SM.ExponentialAtmosphereModel(spice_earth), SM.SpiceEphemeridesModel();
                   ra_alt_m = 200e3, rp_alt_m = 160e3)
    end)

    # The case the density-family signature gets wrong: vacuum, so it signs as
    # `dens=none`, yet every third-body sample enters the shared native lock.
    push!(out, "spice_nbody" => () -> begin
        spice_earth = SM.Earth(HARM_FILE)
        nbody = SM.NBodyGravityModel(body_names = ("Sun", "Moon"),
                                     primary_body_name = "Earth", planet = spice_earth)
        build_case(spice_earth,
                   (SM.GravitationalHarmonicsModel(HARM_DEGREE, HARM_DEGREE, HARM_FILE, spice_earth), nbody),
                   SM.NoAtmosphereModel(), SM.SpiceEphemeridesModel();
                   ra_alt_m = 400e3, rp_alt_m = 380e3)
    end)

    if GRAM_AVAILABLE[]
        push!(out, "gram" => () -> build_case(
            earth, (harm(), SM.AerodynamicCoefficientfM()),
            SM.GRAMAtmosphereModel(planet_name = "earth"), SM.SimpleEphemeridesModel();
            ra_alt_m = 200e3, rp_alt_m = 160e3))
    end
    return out
end

function run_family(name::String, builder::Function)
    case = try
        builder()
    catch err
        @warn "Could not build family; skipping." family = name exception = err
        return nothing
    end

    counts = try
        PC.constellation_work_counts(case.args, N_SATS)
    catch
        nothing
    end
    workspace = counts === nothing ? 0.0 : counts.simd_workspace_bytes_per_sat

    before = RS.native_lock_stats_snapshot()
    ci = try
        SE.probe_contention_inputs!(case.p, case.u0, case.args;
                                    k = 7, t_step = T_STEP,
                                    workspace_bytes_per_sat = workspace)
    catch err
        @warn "Probe failed; skipping." family = name exception = err
        return nothing
    end
    if !ci.valid
        # A bare "no evidence" hides whether the workload is unmeasurable or
        # simply broken, and those need different responses.
        try
            SE.spacecraft_dynamics!(zero(case.u0), case.u0, case.p, 0.0)
            @warn "Probe reported no evidence though the RHS evaluates." family = name
        catch err
            @warn "Probe reported no evidence; the RHS itself fails." family = name exception = err
        end
        return nothing
    end

    after = RS.native_lock_stats_snapshot()
    top_site = :none
    top_hold = 0
    for site in RS.NATIVE_LOCK_SITES
        # Differenced across the probe window: building the case furnishes
        # kernels and constructs planets, and charging that to the workload
        # would report every SPICE family as spice_body-dominated.
        held = getproperty(after.sites, site).hold_ns - getproperty(before.sites, site).hold_ns
        if held > top_hold
            top_hold = held
            top_site = site
        end
    end

    return (
        name = name,
        n_effectors = length(case.args.dynamics_model.dynamic_effectors),
        pass_ms = ci.pass_ns / 1e6,
        alloc_kib = ci.alloc_bytes_per_pass / 1024,
        alloc_mibs = PC.alloc_bytes_per_second(ci) / 1024^2,
        rho = PC.lock_duty_cycle(ci),
        wait_hold = PC.lock_wait_hold_ratio(ci),
        acq = ci.lock_acquisitions,
        top_site = top_site,
        in_domain = counts === nothing ? false : counts.in_domain,
    )
end

function main()
    println("contention probe matrix   N=", N_SATS, "  L=", HARM_DEGREE,
            "  threads=", Threads.nthreads())
    println()
    @printf("%-12s %5s %10s %12s %12s %9s %10s %8s %-16s %s\n",
            "family", "effs", "pass_ms", "alloc_KiB", "alloc_MiB/s", "rho",
            "wait/hold", "lock_acq", "top_lock_site", "in_domain")
    println("-"^116)
    rows = Any[]
    for (name, builder) in families()
        RS.reset_native_lock_stats!()
        row = run_family(name, builder)
        row === nothing && continue
        push!(rows, row)
        @printf("%-12s %5d %10.3f %12.1f %12.2f %9.4f %10.3f %8d %-16s %s\n",
                row.name, row.n_effectors, row.pass_ms, row.alloc_kib, row.alloc_mibs,
                row.rho, row.wait_hold, row.acq, String(row.top_site), row.in_domain)
    end
    println()

    # The gate, stated as the comparison it is rather than left to the reader.
    byname = Dict(r.name => r for r in rows)
    if haskey(byname, "gravity") && haskey(byname, "fullstack")
        g, f = byname["gravity"], byname["fullstack"]
        println("allocation separation  gravity=", round(g.alloc_kib, digits = 1),
                " KiB/pass  vs  fullstack=", round(f.alloc_kib, digits = 1),
                " KiB/pass  ->  ", f.alloc_kib > 4 * max(g.alloc_kib, 1e-9) ? "SEPARATED" : "NOT SEPARATED")
    end
    if haskey(byname, "spice_nbody")
        # The claim under test is NOT that spice_nbody holds the lock more than
        # GRAM does -- it does not, and it should not. It is that the
        # density-model family cannot tell spice_nbody from gravity or aero,
        # because all three sign as a non-GRAM density, while measured
        # occupancy separates them cleanly.
        s = byname["spice_nbody"]
        peers = [r for r in rows if r.name in ("gravity", "aero")]
        worst = isempty(peers) ? 0.0 : maximum(r.rho for r in peers)
        println("occupancy separation   spice_nbody rho=", round(s.rho, digits = 4),
                " (signs as dens=none)  vs  same-signature peers rho<=", round(worst, digits = 4),
                "  ->  ", s.rho > max(2 * worst, 1e-4) ? "SEPARATED" : "NOT SEPARATED")
        if haskey(byname, "gram")
            println("                       gram rho=", round(byname["gram"].rho, digits = 4),
                    " -- the family a density-keyed signature does get right")
        end
    end
end

GRAM_AVAILABLE[] = ensure_gramsuite_loaded!()
main()
