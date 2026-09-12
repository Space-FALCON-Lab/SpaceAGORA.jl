# GRAM multi-satellite micro-benchmark.
#
# Compares, for an N-satellite constellation stepping through time:
#   A. baseline  : one native GRAM call per (step, satellite) -- GRAM computes
#                  its own ephemeris via CSPICE on every call.
#   B. injected  : one native ephemeris computation per *step* (shared across
#                  the constellation), then N CSPICE-free GRAM calls.
#   C. injected+threads : B, with one cloned GRAM atmosphere per thread and a
#                  per-instance lock instead of the process-global lock.
#
# Run: julia --project=<repo> -t <N> gram_micro2.jl

using Printf, Base.Threads

const REPO = "/home/space-falcon-1/Documents/SpaceAGORA.jl"

import Pkg
Pkg.activate(REPO; io=devnull)
if Base.find_package("GRAMSuite") === nothing
    pushfirst!(LOAD_PATH, joinpath(REPO, "data", "GRAMSuite.jl"))
end
using GRAMSuite

gramsuite_gram() = getfield(GRAMSuite, :_GRAM_WRAPPER)[]
gcall(G, name::Symbol, args...; kw...) =
    Base.invokelatest(Base.invokelatest(getfield, G, name), args...; kw...)

function rss_mb()
    for line in eachline("/proc/self/status")
        startswith(line, "VmRSS:") && return parse(Int, split(line)[2]) / 1024
    end
    return NaN
end

build_model(planet) = GRAMSuite.GRAMAtmosphereModel(
    planet_name = planet,
    gram_root_directory = joinpath(REPO, "data", "GRAMSuite.jl", "GRAM Suite 2.0"),
    gram_library_path = get(ENV, "GRAM_LIB_OVERRIDE", ""),
    initial_time = GRAMSuite.InitialTime(year=2020, month=6, day=1, hour=12, minute=0, second=0.0),
)

# ---------------------------------------------------------------------------
# Constellation ground track: nsat satellites, nsteps epochs, 5 s cadence.
# ---------------------------------------------------------------------------
@inline function track_point(sat::Int, k::Int, nsat::Int)
    t   = 5.0 * k
    lat = deg2rad(50.0 * sin(2pi * (k / 1100.0) + 2pi * sat / nsat))
    lon = deg2rad(mod(360.0 * (k / 1100.0) + 360.0 * sat / nsat, 360.0))
    h   = 300.0e3 + 20.0e3 * sin(2pi * (k / 550.0) + sat)
    return h, lat, lon, t
end

# ---------------------------------------------------------------------------
# Ephemeris injection: body-level fields shared, per-satellite fields derived.
# ---------------------------------------------------------------------------
@inline function inject_state(EphT, base, planet::String, lat_rad::Float64, lon_rad::Float64)
    lon_deg = rad2deg(lon_rad)

    # GRAM's own local-solar-time formula (common/source/Ephemeris.cpp:438-448)
    st = 12.0 + (lon_deg - base.subsolarLongitude) / 15.0
    if planet == "venus" || planet == "uranus"
        st = 24.0 - st
    end
    st = st > 24.0 ? st - 24.0 : (st < 0.0 ? st + 24.0 : st)

    sslat = deg2rad(base.subsolarLatitude)
    sslon = deg2rad(base.subsolarLongitude)
    csza  = sin(lat_rad)*sin(sslat) + cos(lat_rad)*cos(sslat)*cos(lon_rad - sslon)
    sza   = rad2deg(acos(clamp(csza, -1.0, 1.0)))

    return EphT(st, base.longitudeSun, base.subsolarLatitude, base.subsolarLongitude,
                base.orbitalRadius, base.oneWayLightTime, sza, base.secondsPerSol)
end

# Function handles bound once from the runtime-included wrapper module, so the
# hot loops make ordinary (not invokelatest / not dynamic-getproperty) calls.
struct GramAPI{A,B,C,D,E,T}
    set_position!::A
    update!::B
    get_dynamics_state::C
    get_ephemeris_state::D
    set_ephemeris_state!::E
    EphT::Type{T}
end

function bind_api(G)
    gf(n) = Base.invokelatest(getfield, G, n)
    return GramAPI(gf(:set_position!), gf(:update!), gf(:get_dynamics_state),
                   gf(:get_ephemeris_state), gf(:set_ephemeris_state!),
                   gf(:EphemerisStateC))
end

# One native GRAM evaluation (GRAM computes ephemeris internally via CSPICE).
@inline function native_eval!(api, atm, h, lat, lon, t)
    api.set_position!(atm; height=h*1e-3, latitude=rad2deg(lat),
                      longitude=rad2deg(lon), elapsed_time=t)
    err = api.update!(atm)
    err == 0 || error("GRAM update failed: $err")
    return api.get_dynamics_state(atm).density
end

# One injected GRAM evaluation (no CSPICE inside GRAM).
@inline function injected_eval!(api, atm, base, planet, h, lat, lon, t)
    api.set_ephemeris_state!(atm, inject_state(api.EphT, base, planet, lat, lon))
    api.set_position!(atm; height=h*1e-3, latitude=rad2deg(lat),
                      longitude=rad2deg(lon), elapsed_time=t)
    err = api.update!(atm)
    err == 0 || error("GRAM update failed: $err")
    return api.get_dynamics_state(atm).density
end

# ---------------------------------------------------------------------------
# A. baseline
# ---------------------------------------------------------------------------
function run_baseline(api, model, nsteps, nsat; collect_rho=false)
    atm = model.gram_atmosphere
    out = collect_rho ? Vector{Float64}(undef, nsteps*nsat) : Float64[]
    GC.gc(true); rss0 = rss_mb(); acc = 0.0
    t0 = time_ns()
    for k in 1:nsteps, sat in 1:nsat
        h, lat, lon, t = track_point(sat, k, nsat)
        rho = native_eval!(api, atm, h, lat, lon, t)
        acc += rho
        collect_rho && (out[(k-1)*nsat + sat] = rho)
    end
    dt = (time_ns()-t0)*1e-9
    GC.gc(true)
    return (; dt, rss = rss_mb()-rss0, acc, out)
end

# ---------------------------------------------------------------------------
# B. injected, single-threaded
# ---------------------------------------------------------------------------
function run_injected(api, model, nsteps, nsat, planet; collect_rho=false)
    atm = model.gram_atmosphere
    out = collect_rho ? Vector{Float64}(undef, nsteps*nsat) : Float64[]
    GC.gc(true); rss0 = rss_mb(); acc = 0.0
    t0 = time_ns()
    for k in 1:nsteps
        # ---- one native ephemeris computation for the whole constellation ----
        h1, lat1, lon1, t1 = track_point(1, k, nsat)
        native_eval!(api, atm, h1, lat1, lon1, t1)
        base = api.get_ephemeris_state(atm)
        # ---- N CSPICE-free evaluations ----
        for sat in 1:nsat
            h, lat, lon, t = track_point(sat, k, nsat)
            rho = injected_eval!(api, atm, base, planet, h, lat, lon, t)
            acc += rho
            collect_rho && (out[(k-1)*nsat + sat] = rho)
        end
    end
    dt = (time_ns()-t0)*1e-9
    GC.gc(true)
    return (; dt, rss = rss_mb()-rss0, acc, out)
end

# ---------------------------------------------------------------------------
# C. injected + per-thread cloned atmospheres
# ---------------------------------------------------------------------------
function run_injected_threaded(api, G, model, nsteps, nsat, planet, nworkers)
    atm0 = model.gram_atmosphere

    # Clone one native atmosphere per worker (copyAtmosphere_C -- no re-init).
    copy_atmosphere = Base.invokelatest(getfield, G, :copy_atmosphere)
    tclone = time_ns()
    atms = [copy_atmosphere(atm0) for _ in 1:nworkers]
    clone_s = (time_ns()-tclone)*1e-9

    # Partition satellites into exactly `nworkers` contiguous chunks, one per
    # cloned atmosphere. Indexing by threadid() is WRONG here: two threads
    # sharing an atmosphere interleave set_ephemeris_state!/update!, so one
    # thread consumes the other's userSuppliedEphemeris flag and the loser
    # falls through to native CSPICE -- which corrupts SPICE's trcpkg call
    # stack and hard-crashes the process. One instance per worker, always.
    chunks = [collect(w:nworkers:nsat) for w in 1:nworkers]

    # Diagnostic: MICRO_REFRESH=frozen computes the body ephemeris once up front
    # instead of once per step, so the timed region contains no native (CSPICE)
    # GRAM call at all. Isolates "workers racing each other" from "workers
    # racing the main thread's per-step refresh".
    frozen = get(ENV, "MICRO_REFRESH", "perstep") == "frozen"
    base_frozen = if frozen
        h1, lat1, lon1, t1 = track_point(1, 1, nsat)
        native_eval!(api, atm0, h1, lat1, lon1, t1)
        api.get_ephemeris_state(atm0)
    else
        nothing
    end

    accs = zeros(Float64, nworkers)
    GC.gc(true); rss0 = rss_mb()
    t0 = time_ns()
    for k in 1:nsteps
        base = if frozen
            base_frozen
        else
            h1, lat1, lon1, t1 = track_point(1, k, nsat)
            native_eval!(api, atm0, h1, lat1, lon1, t1)
            api.get_ephemeris_state(atm0)
        end

        @sync for w in 1:nworkers
            Threads.@spawn begin
                atm_w = atms[w]
                s = 0.0
                for sat in chunks[w]
                    h, lat, lon, t = track_point(sat, k, nsat)
                    s += injected_eval!(api, atm_w, base, planet, h, lat, lon, t)
                end
                accs[w] += s
            end
        end
    end
    dt = (time_ns()-t0)*1e-9
    GC.gc(true)
    return (; dt, rss = rss_mb()-rss0, acc = sum(accs), clone_s)
end

# ---------------------------------------------------------------------------
function run_all(model, planet, nsat, nsteps)
    nthr = nthreads()
    G = gramsuite_gram()
    api = bind_api(G)
    @printf("planet=%s  nsat=%d  nsteps=%d  threads=%d\n", planet, nsat, nsteps, nthr)
    @printf("model built, RSS %.0f MB\n\n", rss_mb())

    # warmup / MERRA2 file read
    run_baseline(api, model, 5, nsat); run_injected(api, model, 5, nsat, planet)

    println("--- A. baseline: native GRAM per (step, satellite) ---")
    a = run_baseline(api, model, nsteps, nsat)
    @printf("   %7.3f s   %8.0f eval/s   RSS %+.1f MB\n\n",
            a.dt, nsteps*nsat/a.dt, a.rss)

    println("--- B. ephemeris-injected, 1 native ephemeris per step ---")
    b = run_injected(api, model, nsteps, nsat, planet)
    @printf("   %7.3f s   %8.0f eval/s   RSS %+.1f MB   speedup %.2fx\n\n",
            b.dt, nsteps*nsat/b.dt, b.rss, a.dt/b.dt)

    println("--- C. injected + per-thread cloned atmospheres ---")
    ladder = haskey(ENV,"MICRO_WORKERS") ? [parse(Int,x) for x in split(ENV["MICRO_WORKERS"],",")] : unique([1, 2, 4, 8, min(nthr,16)])
    for w in ladder
        w > nthr && continue
        c = run_injected_threaded(api, G, model, nsteps, nsat, planet, w)
        @printf("   workers=%2d  %7.3f s  %8.0f eval/s  RSS %+.1f MB  clone %.3f s  vs-baseline %.2fx\n",
                w, c.dt, nsteps*nsat/c.dt, c.rss, c.clone_s, a.dt/c.dt)
    end
    println()

    println("--- D. accuracy parity (per-step injection vs native) ---")
    ap = run_baseline(api, model, 200, nsat; collect_rho=true)
    bp = run_injected(api, model, 200, nsat, planet; collect_rho=true)
    rel = abs.(bp.out .- ap.out) ./ max.(abs.(ap.out), 1e-300)
    @printf("   n=%d  max rel diff %.3e   mean rel diff %.3e\n",
            length(rel), maximum(rel), sum(rel)/length(rel))
end

# Build the model at top level: this Base.include()s the vendored GRAM wrapper
# and so advances the world age. Everything after runs in the newer world.
let planet = get(ENV, "MICRO_PLANET", "earth"),
    nsat   = parse(Int, get(ENV, "MICRO_NSAT", "16")),
    nsteps = parse(Int, get(ENV, "MICRO_NSTEPS", "600"))
    model = build_model(planet)
    Base.invokelatest(run_all, model, planet, nsat, nsteps)
end
