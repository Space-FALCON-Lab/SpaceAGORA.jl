# Analytic work counts for the built-in dynamic effectors.
#
# Discipline for this file: an effector gets an `effector_cost_terms` method ONLY
# when its kernel has been read and the count verified against the actual loop
# structure. Everything else deliberately falls through to the `::Any` default
# and is measured by the run-start probe instead. A wrong analytic count is worse
# than no analytic count -- the probe is honest about what it does not know,
# whereas a fabricated formula is confidently wrong and silently mis-routes.

@inline function _alf_recurrence_iterations(L::Int, M::Int)::Int
    # Associated-Legendre recurrence rows in _harmonics_flat_batch_kernel! phase 4.
    #
    # The kernel fills rows lazily behind a `max_recur_row` high-water mark, so
    # every row from 3 to L+2 is computed exactly once even though the loop body
    # offers two chances to compute it (as `row` at degree l = r-1 with
    # jmax = min(M'+1, l-1), or as `next_row` at degree l = r-2 with
    # jmax_next = min(M'+1, l)). Both spellings evaluate to min(M'+1, r-2) for
    # row r, which is what makes the closed form below exact rather than an
    # upper bound.
    Mp = max(M, 1) + 1
    total = 0
    for r in 3:(L + 2)
        total += min(Mp, r - 2)
    end
    return total
end

"""
    effector_cost_terms(model::GravitationalHarmonicsModel) -> WorkCounts

Counts for `_harmonics_flat_batch_kernel!`, phase by phase.

The count is driven by `active_orders_by_degree`, not by `(L+1)(L+2)/2`. That
matters more than it looks: a zonal-only field (no tesseral or sectoral terms,
which is what a J2-through-J50 model is) has no active orders at any degree, so
its cost is **linear** in L, while a full field is quadratic. Using the
closed-form triangular count for both would overstate a zonal L=50 model by
better than an order of magnitude and route it as though it were expensive.
"""
function effector_cost_terms(model::GravitationalHarmonicsModel)::WorkCounts
    L = Int(model.L)
    M = Int(model.M)
    L >= 1 || return WorkCounts(scalar_items = 2.0)

    active_total = 0
    @inbounds for row in eachindex(model.active_orders_by_degree)
        row <= L + 1 || continue
        row >= 2 || continue
        active_total += length(model.active_orders_by_degree[row])
    end

    alf = _alf_recurrence_iterations(L, M)

    # Vectorized work counted in FLOPS PER LANE, not in loop passes.
    #
    # Counting passes was the first version and it understated the real kernel
    # by 21x, because the passes are nowhere near equal: the sub-diagonal
    # recurrence is 2 flops per lane while the active-order accumulation is 24,
    # and the heavy one is the one that dominates a full field. A model that
    # weights them equally is not slightly off, it is measuring a different
    # quantity -- and it failed cold validation at 25% decision accuracy with the
    # error concentrated exactly where the mix of pass types changes, which is
    # across degree and order.
    #
    # Weights below are counted off `_harmonics_flat_batch_kernel!` body by
    # body. Multiplies and adds count 1 each; a muladd counts 2, matching the
    # calibration kernel's own accounting so the rate and the count are in the
    # same units.
    #
    #   phase 2 seed          A[b,2,1] = u*sqrt3                          1
    #   phase 2 row           A = u * s2n3 * A                            2
    #   phase 3 step          R,I complex rotate: 4 mul + 2 add           6
    #   phase 4 ALF           A = u*N1*A - N2*A                           4
    #   phase 4 rho reset     rho *= rho; rr = rho/RE; 4 zero stores      2
    #   phase 4 zonal         2 x (VR * A * D0) accumulate                6
    #   phase 4 active order  D,E,F (12) + mA (1) + 4 accumulates (11)   24
    #   phase 4 tail          4 x (rr * sum) accumulate                    8
    simd_terms = Float64(
        1 + 2 * (L + 1) +
        6 * (M + 1) +
        4 * alf +
        (2 + 6 + 8) * L +
        24 * active_total
    )

    # Phase 1 gather and phase 5 back-transform, counted in flops per satellite
    # like the vectorized terms -- not as "two passes", which is what they were
    # and which understated them by roughly forty-fold.
    #
    # These are the only parts of the kernel that are NOT vectorized (plain
    # `@inbounds for b = 1:B`, no `@turbo`), and they are dense: a 3x3
    # matrix-vector product, a norm with its square root, a reciprocal, and a
    # second matrix-vector product on the way out. For a zonal field they
    # dominate the whole kernel, which is why treating them as negligible made
    # the model predict that a 1024-satellite L20 zonal solve was too cheap to
    # be worth parallelising -- when measurement says the opposite.
    #
    #   phase 1: lpi * pos (15), norm incl. sqrt (16), reciprocal (10),
    #            s/t/u scaling (3), rho and rho_np1 (3)               ~47
    #   phase 5: g_pp assembly (6), central term incl. inv_r^3 (9),
    #            lpi' * g_pp (15), mass scaling (3), 3 accumulates (3) ~36
    scalar_items = 83.0

    # N1/N2 per ALF step; C/VR01/VR11 for each zonal term; C/S/VR01/VR11 per
    # active order. Counted as reads, not bytes: the coefficient matrices are
    # column-major and the kernel walks a row at fixed `row` over varying `j`,
    # so each read strides by (L+2) elements and the traffic per read is set by
    # cache-line granularity and table size rather than by sizeof(Float64).
    # ns_per_coeff_touch absorbs that, which is why it is calibrated rather than
    # derived.
    coeff_touches = Float64(2 * alf + 3 * L + 4 * active_total)

    # Footprint the row walk strides over, which decides the cache level and
    # therefore the per-touch rate. Six arrays are interleaved at the same
    # (row, j): C, S, VR01, VR11 in the order accumulation and N1, N2 in the ALF
    # recurrence, each an (L+2) x (M+2) Float64 matrix. Counting all six is the
    # working set the walk actually pulls through cache, not the size of any one
    # of them. The cold-prediction check against the real kernel is what would
    # expose this being the wrong aggregate.
    table_bytes = 6.0 * (L + 2) * (M + 2) * 8.0

    # The A workspace the batch kernel carries per satellite: (L+3) x (M+2)
    # Float64. Multiplied by batch width this is what the SIMD lane rate is
    # indexed on, because at wide batch it leaves cache entirely.
    workspace_per_sat = Float64((L + 3) * (M + 2) * 8)

    return WorkCounts(
        simd_terms = simd_terms,
        scalar_items = scalar_items,
        coeff_touches = coeff_touches,
        coeff_table_bytes = table_bytes,
        simd_workspace_bytes_per_sat = workspace_per_sat,
    )
end

# Closed-form central-body terms: a handful of flops per satellite, no table.
@inline effector_cost_terms(::InverseSquaredGravityModel)::WorkCounts =
    WorkCounts(scalar_items = 1.0)
@inline effector_cost_terms(::InverseSquaredJ2GravityModel)::WorkCounts =
    WorkCounts(scalar_items = 2.0)
@inline effector_cost_terms(::ConstantGravityModel)::WorkCounts =
    WorkCounts(scalar_items = 1.0)

"""
    effector_cost_terms(model::NBodyGravityModel) -> WorkCounts

One third-body acceleration per body per satellite. The ephemeris lookups
themselves are hoisted out of the satellite loop and cached per pass (see
`pos_primary_k_all` in `perturbations.jl`), so they are not counted per
satellite; only the per-body arithmetic is.
"""
@inline function effector_cost_terms(model::NBodyGravityModel)::WorkCounts
    return WorkCounts(scalar_items = Float64(max(1, length(model.body_names))))
end

# Flat-queue node predicates, mirrored from SimulationEngine.
#
# SOURCE OF TRUTH: `_batchable_effector` / `_harmonics_prepass_effector` /
# `_count_flat_queue_only_effectors` in src/simulation/engine/dynamics_rhs.jl.
# They are duplicated here rather than called because ParallelCost is included
# from core/simulation_model.jl, which the engine is built on top of -- calling
# upward would invert the dependency.
#
# Duplication is made safe by test/unit/parallel/cost_work_counts_tests.jl,
# which asserts these agree with the engine's versions for every built-in
# effector, including the `gravity_gradient` variants. Divergence is a test
# failure, not a silent mis-count.
@inline _cost_batchable_effector(::Any)::Bool = false
@inline _cost_batchable_effector(::NBodyGravityModel)::Bool = true
@inline _cost_batchable_effector(::SolarRadiationPressureModel)::Bool = true
@inline _cost_batchable_effector(e::InverseSquaredGravityModel)::Bool = !e.gravity_gradient
@inline _cost_batchable_effector(e::InverseSquaredJ2GravityModel)::Bool = !e.gravity_gradient

@inline _cost_harmonics_prepass_effector(::Any)::Bool = false
@inline _cost_harmonics_prepass_effector(::GravitationalHarmonicsModel)::Bool = true

"""
    flat_queue_node_effector(effector; partition_active = false) -> Bool

Whether an effector produces flat-queue work items.

Effectors resolved by a pre-pass -- the batchable ones that write straight into
`totals` from position buffers, and harmonics with its own SIMD batch -- have
already written their contribution by the time the queue runs, and the queue
skips them. Counting them as nodes overstates `queue_nodes` by a factor of
(total effectors)/(queue-only effectors), which for a vacuum harmonics-only
constellation means counting N nodes where the queue does zero work.

`partition_active` reverses the exclusion: under a solver partition (split
IMEX) the pre-pass results are not applicable and every selected effector goes
through the queue, which is what the `partition === nothing` guard in
`_flat_selection_mask` encodes.
"""
@inline function flat_queue_node_effector(effector; partition_active::Bool = false)::Bool
    partition_active && return true
    return !(_cost_batchable_effector(effector) || _cost_harmonics_prepass_effector(effector))
end

"""
    model_in_cost_domain(density_model) -> Bool

Whether the linear cost model can represent this density model's contribution.

Returns false for native point GRAM, whose queries serialize behind a
process-wide lock. Lock contention is superlinear in concurrency and has no
representation as a sum of per-unit rates, so no calibration of `w_k` can make
`sum w_k n_k` fit it -- the measured signature is a flat curve from 2 to 10
threads and a 2.6x jump at 12, where the RHS plan flips, not a smooth rise. A
model that answered here would be confidently wrong exactly where the cost of
being wrong is highest (+131% regret on `atmo256_gram_live_10min`).

Out-of-domain workloads are not mispredicted, they are declined: the predictor
abstains and the caller keeps its existing heuristic route.

The lock-free GRAM surrogate is in domain -- it has no shared lock, which is
why it is a separate policy source from native GRAM in the first place.
"""
@inline function model_in_cost_domain(density_model)::Bool
    return !(density_model isa EnvironmentModels.GRAMAtmosphereModel)
end

"""
    constellation_work_counts(args, n_active_sats; probe = nothing) -> WorkCounts

Sum the per-pass work of every dynamic effector in `args` over a constellation
of `n_active_sats` active satellites.

Effectors that declare `effector_cost_terms` contribute analytic counts.
Effectors that do not are looked up in `probe` -- a `Dict` of measured
per-satellite nanoseconds keyed by effector index, produced at run start by
`probe_effector_costs!`.

An effector that is in neither clears `in_domain`. It is tempting to let it
contribute zero and carry on, but zero is a *number*, and downstream the
difference between "this effector is free" and "nobody measured this effector"
would be invisible: the predictor would confidently rank a workload whose
dominant cost it never saw. That is the same silent-underestimate failure the
declare-or-probe design exists to avoid, so an unmeasured unknown effector
disqualifies the whole prediction rather than quietly biasing it. In a wired-up
run this arises only when the probe itself failed -- a throwing effector, or an
RHS warm-up that could not complete -- which is exactly when the model should
decline to answer.
"""
function constellation_work_counts(
    args,
    n_active_sats::Int;
    probe::Union{Nothing, Dict{Int, Float64}} = nothing,
    partition_active::Bool = false,
)::WorkCounts
    n_sats = max(0, n_active_sats)
    effectors = args.dynamics_model.dynamic_effectors
    total = WorkCounts()
    queue_effectors = 0

    @inbounds for idx in eachindex(effectors)
        effector = effectors[idx]
        flat_queue_node_effector(effector; partition_active=partition_active) &&
            (queue_effectors += 1)

        declared = effector_cost_terms(effector)
        if declared !== nothing
            total = total + declared
            continue
        end
        measured = probe === nothing ? nothing : get(probe, idx, nothing)
        usable = measured !== nothing && isfinite(measured) && measured >= 0.0
        total = total + WorkCounts(
            probe_ns = usable ? measured : 0.0,
            unknown_effectors = 1,
            in_domain = usable,
        )
    end

    return WorkCounts(
        simd_terms = total.simd_terms,
        scalar_items = total.scalar_items,
        coeff_touches = total.coeff_touches,
        coeff_table_bytes = total.coeff_table_bytes,
        simd_workspace_bytes_per_sat = total.simd_workspace_bytes_per_sat,
        queue_nodes = Float64(n_sats * queue_effectors),
        probe_ns = total.probe_ns,
        unknown_effectors = total.unknown_effectors,
        in_domain = total.in_domain &&
            model_in_cost_domain(args.environment_model.density_model),
    )
end
