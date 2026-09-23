# WS11b step 2: does the isolated GRAM pool return the same density, to the bit,
# as the single locked GRAM instance?
#
# The pool replaces one shared native GRAM model, serialized on the process-wide
# lock, with `workers` independent `deepcopy`ed models, each behind its own lock.
# Whether that is admissible is a pure value question: for the same (altitude,
# latitude, longitude, elapsed time) and the same epoch, an instance that has
# seen a different sequence of earlier calls must still return the same density,
# temperature and wind.
#
# Three comparisons, all over the same grid of states, all exact (`===` on the
# bits, via `reinterpret`), never a tolerance:
#
#   instance   every pool instance evaluated serially, one point at a time, in
#              the same order, against the template instance. Isolates "a second
#              GRAM instance disagrees with the first" from anything threading
#              does.
#   batch      the shipped locked batch call (`getDensityBatch!`, which for a
#              GRAMAtmosphereModel loops `getDensity` under the global lock)
#              against the shipped pool batch call
#              (`_gram_isolated_pool_batch_eval!`). This is exactly the A/B the
#              density callback performs at runtime.
#   replay     the locked batch call run twice. A difference here means native
#              GRAM is itself path-dependent and no pool arrangement can be
#              bit-identical; it is the control that tells the other two rows
#              apart from a property of GRAM.
#
# Usage:
#   julia --project=. --threads=4 \
#       benchmarks/studies/gram_thread_scaling/density_grid_parity.jl \
#       [--npoints=512] [--workers=4] [--planet=earth]

using Printf
using StaticArrays

const GTS_DIR = @__DIR__
const GTS_REPO_ROOT = normpath(joinpath(GTS_DIR, "..", "..", ".."))
const GTS_PPC_DIR = joinpath(GTS_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(GTS_PPC_DIR, "cli.jl"))
include(joinpath(GTS_PPC_DIR, "modes.jl"))
include(joinpath(GTS_PPC_DIR, "cases.jl"))

# Top level, not inside a function: `ppc_ensure_gramsuite_loaded!` runs
# `@eval import GRAMSuite`, and everything compiled in the same world age as
# that call cannot see the resulting bindings. cases.jl calls it the same way.
ppc_ensure_gramsuite_loaded!()

const CB = SpaceAGORA.SimulationModel.SimulationCallbacks
const EM = SpaceAGORA.SimulationModel.EnvironmentModels

gts_arg(args, key, default) = begin
    for a in args
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end

const GTS_NPOINTS = parse(Int, gts_arg(ARGS, "npoints", "512"))
const GTS_WORKERS = parse(Int, gts_arg(ARGS, "workers", string(max(2, Threads.nthreads()))))
const GTS_PLANET = gts_arg(ARGS, "planet", "earth")

# Bit comparison. `==` on Float64 calls -0.0 equal to 0.0 and every NaN unequal
# to itself, so neither direction of it answers "is this the same Float64".
@inline bits(x::Float64) = reinterpret(UInt64, x)
@inline same_bits(a::Float64, b::Float64) = bits(a) === bits(b)

function first_bit_difference(a::Vector{Float64}, b::Vector{Float64})
    for i in eachindex(a)
        same_bits(a[i], b[i]) || return (i, a[i], b[i])
    end
    return nothing
end

function first_bit_difference(a::Vector{SVector{3, Float64}}, b::Vector{SVector{3, Float64}})
    for i in eachindex(a)
        for k in 1:3
            same_bits(a[i][k], b[i][k]) || return (i, a[i][k], b[i][k])
        end
    end
    return nothing
end

"""
    gts_grid(n)

A deterministic sweep of states in the band where the aero constellations this
study measures actually sample the atmosphere. Altitude spans the drag regime
(the 64-spacecraft reference constellation spans roughly 500 to 670 km, and
the sweep reaches down to 150 km), latitude the full pole-to-pole range,
longitude a full revolution, and elapsed time the 100 s mission those traces
use. The three angular sweeps are given mutually irrational strides so that no
two grid points repeat a (latitude, longitude) pair -- a repeat would let a
cached native lookup hide a real disagreement.
"""
function gts_grid(n::Int)
    hs = Vector{Float64}(undef, n)
    lats = Vector{Float64}(undef, n)
    lons = Vector{Float64}(undef, n)
    ts = Vector{Float64}(undef, n)
    for i in 1:n
        x = (i - 1) / max(1, n - 1)
        hs[i] = 150.0e3 + x * 550.0e3
        lats[i] = -0.5pi + pi * mod(x * sqrt(2.0), 1.0)
        lons[i] = -pi + 2pi * mod(x * sqrt(3.0), 1.0)
        ts[i] = x * 100.0
    end
    return hs, lats, lons, ts
end

function gts_params(planet_name::String)
    planet = Earth("", PPC_SPICE_PATH)
    args = ppc_build_config(
        planet=planet,
        spacecraft=ppc_constellation(planet, 1),
        mission_time_s=100.0,
        orientation_sim=false,
        dynamic_effectors=(ppc_harmonics_model(planet, 50), AerodynamicCoefficientfM()),
        density_model=ppc_gram_atmosphere_model(planet_name),
        dt_max_orbit=5.0
    )
    return args, ODEParams(n_sats=GTS_NPOINTS, args=args)
end

function main()
    planet_name = GTS_PLANET
    planet_name == "earth" || error("Only the earth case is wired here; got '$planet_name'.")
    args, p = gts_params(planet_name)
    model = args.environment_model.density_model
    hs, lats, lons, ts = gts_grid(GTS_NPOINTS)
    n = GTS_NPOINTS

    @printf("grid      n=%d  h=[%.1f, %.1f] km  t=[%.1f, %.1f] s  threads=%d  workers=%d\n",
            n, hs[1] * 1e-3, hs[end] * 1e-3, ts[1], ts[end], Threads.nthreads(), GTS_WORKERS)

    alloc() = (zeros(Float64, n), zeros(Float64, n),
               [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n])

    # -- locked reference (the shipped path) ---------------------------------
    rho_ref, T_ref, w_ref = alloc()
    getDensityBatch!(rho_ref, T_ref, w_ref, model, hs, lats, lons, ts, true, p)
    nonzero = count(!iszero, rho_ref)
    @printf("locked    rho[1]=%.17g  nonzero=%d/%d\n", rho_ref[1], nonzero, n)
    nonzero > 0 || error("Every reference density is zero: the grid never reaches GRAM.")

    failures = String[]

    # -- replay control ------------------------------------------------------
    rho_2, T_2, w_2 = alloc()
    getDensityBatch!(rho_2, T_2, w_2, model, hs, lats, lons, ts, true, p)
    for (label, d) in (("rho", first_bit_difference(rho_ref, rho_2)),
                       ("T", first_bit_difference(T_ref, T_2)),
                       ("wind", first_bit_difference(w_ref, w_2)))
        if d === nothing
            @printf("replay    %-4s identical\n", label)
        else
            @printf("replay    %-4s DIFFERS at i=%d: %.17g vs %.17g\n", label, d[1], d[2], d[3])
            push!(failures, "replay/$label")
        end
    end

    # -- per-instance control ------------------------------------------------
    # `_ensure_gram_isolated_pool!` is the production builder: deepcopy plus the
    # single-threaded warm-up call that keeps fresh clones out of CSPICE.
    models, locks = CB._ensure_gram_isolated_pool!(p, model, GTS_WORKERS)
    @printf("pool      built %d instances\n", length(models))
    for k in eachindex(models)
        rho_k, T_k, w_k = alloc()
        for i in 1:n
            rho_k[i], T_k[i], w_k[i] = CB._gram_isolated_pool_density_state(
                models[k], hs[i], lats[i], lons[i], ts[i], true, p, locks[k]
            )
        end
        d = first_bit_difference(rho_ref, rho_k)
        if d === nothing
            @printf("instance  %d rho identical to the locked reference\n", k)
        else
            @printf("instance  %d rho DIFFERS at i=%d: locked %.17g vs pool %.17g (rel %.3e)\n",
                    k, d[1], d[2], d[3], abs(d[3] - d[2]) / max(abs(d[2]), eps()))
            push!(failures, "instance$k/rho")
        end
        dT = first_bit_difference(T_ref, T_k)
        dT === nothing || (push!(failures, "instance$k/T");
                           @printf("instance  %d T DIFFERS at i=%d: %.17g vs %.17g\n", k, dT[1], dT[2], dT[3]))
        dw = first_bit_difference(w_ref, w_k)
        dw === nothing || (push!(failures, "instance$k/wind");
                           @printf("instance  %d wind DIFFERS at i=%d: %.17g vs %.17g\n", k, dw[1], dw[2], dw[3]))
    end

    # -- the shipped pool batch call, threaded -------------------------------
    rho_p, T_p, w_p = alloc()
    pooled = withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "on",
                     "SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS" => string(GTS_WORKERS)) do
        CB._gram_isolated_pool_batch_eval!(
            rho_p, T_p, w_p, model, hs, lats, lons, ts, true, p;
            allotment_hint=GTS_WORKERS
        )
    end
    if !pooled
        println("batch     pool path declined the call (threads=$(Threads.nthreads()), workers=$(GTS_WORKERS)); nothing measured")
        push!(failures, "batch/declined")
    else
        for (label, d) in (("rho", first_bit_difference(rho_ref, rho_p)),
                           ("T", first_bit_difference(T_ref, T_p)),
                           ("wind", first_bit_difference(w_ref, w_p)))
            if d === nothing
                @printf("batch     %-4s identical\n", label)
            else
                @printf("batch     %-4s DIFFERS at i=%d: locked %.17g vs pool %.17g (rel %.3e)\n",
                        label, d[1], d[2], d[3], abs(d[3] - d[2]) / max(abs(d[2]), eps()))
                push!(failures, "batch/$label")
            end
        end
    end

    println()
    if isempty(failures)
        println("RESULT: bit-identical on every comparison")
    else
        println("RESULT: NOT bit-identical -- $(join(failures, ", "))")
    end
    return isempty(failures)
end

main() || exit(1)
