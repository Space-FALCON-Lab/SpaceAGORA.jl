# Standalone single-point CLI for leo_constellation_size_scaling_point.jl's
# run_scaling_point. Thin wrapper kept for ad-hoc single-point runs/debugging;
# the full sweep (leo_constellation_size_scaling.jl) calls run_scaling_point
# directly in-process instead of spawning this as a subprocess per point --
# see that script's header and leo_constellation_size_scaling_point.jl's
# header for why a subprocess per point is no longer necessary.
#
#   julia --project=. --threads=<N> leo_constellation_size_scaling_worker.jl <n_sats> <mode> [route]
#
# mode is one of:
#   standard  -- real/native GRAM, no vacuum-predicted lookahead cache
#                (SPACEAGORA_DENSITY_FREEZE_PER_STEP=1, density calls
#                serialized -- see leo_2048_constellation_gram_scaling.jl's
#                own header for why concurrent direct native GRAM calls are
#                unsafe). Measured faster than the lookahead cache for this
#                600s mission length (see project findings): the lookahead
#                cache's upfront drag-free-trajectory + spline-fit cost isn't
#                amortized by a mission this short.
#   surrogate -- GRAM offline surrogate (GRAMAtmosphereModelSurrogate):
#                precomputed interpolation grid, no native calls, no lock.
#   no_gram   -- NoAtmosphereModel(), and the aerodynamic effector dropped
#                entirely -- baseline gravity-only dynamics cost, isolating
#                constellation-size scaling with zero atmosphere-model
#                overhead.
#
# route (optional, defaults to "monolithic") is one of:
#   monolithic -- one coupled state vector for all N_SATS satellites, RHS execution
#                 mode "flat" (batched harmonics SIMD kernel, flat effector queue).
#   ensemble   -- SpaceAGORA.run_constellation_ensemble: each satellite becomes an
#                 independent single-satellite run_simulation call dispatched to a
#                 worker task, matching leo_2048_constellation_gram_scaling.jl's own
#                 "ensemble" mode -- RHS execution mode "serial" per member (no
#                 multi-satellite batching to gain within a 1-satellite member), outer
#                 parallelism across members instead of inner parallelism within one.

mode = length(ARGS) >= 2 ? ARGS[2] : ""
mode in ("standard", "surrogate", "no_gram") || error(
    "Usage: julia --project=. --threads=<N> $(basename(@__FILE__)) <n_sats> <standard|surrogate|no_gram> [monolithic|ensemble]"
)
route = length(ARGS) >= 3 ? ARGS[3] : "monolithic"
route in ("monolithic", "ensemble") || error("route must be \"monolithic\" or \"ensemble\", got \"$(route)\"")
n_sats = parse(Int, ARGS[1])

include(joinpath(@__DIR__, "leo_constellation_size_scaling_point.jl"))

run_scaling_point(n_sats, mode, route)
