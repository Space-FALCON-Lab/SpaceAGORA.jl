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

    # Phase 2 sub-diagonal: 1 seed + (L+1) rows.
    # Phase 3 longitude recurrence: 1 seed + (M+1) steps.
    # Phase 4 per degree: accumulator reset, zonal m=0 block, tail accumulate.
    simd_terms = Float64((L + 2) + (M + 2) + alf + 3 * L + active_total)

    # Phase 1 gather and phase 5 back-transform, one scalar pass each.
    scalar_items = 2.0

    # N1/N2 per ALF step; C/VR01/VR11 for each zonal term; C/S/VR01/VR11 per
    # active order. Counted as reads, not bytes: the coefficient matrices are
    # column-major and the kernel walks a row at fixed `row` over varying `j`,
    # so each read strides by (L+2) elements and the traffic per read is set by
    # cache-line granularity and table size rather than by sizeof(Float64).
    # ns_per_coeff_touch absorbs that, which is why it is calibrated rather than
    # derived.
    coeff_touches = Float64(2 * alf + 3 * L + 4 * active_total)

    return WorkCounts(
        simd_terms = simd_terms,
        scalar_items = scalar_items,
        coeff_touches = coeff_touches,
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
`probe_effector_costs!`. An effector that is in neither is counted as unknown
with zero cost, which is only correct when the caller has already decided to
abstain; `probe_effector_costs!` covers every undeclared effector, so this
should not arise in a wired-up run.
"""
function constellation_work_counts(
    args,
    n_active_sats::Int;
    probe::Union{Nothing, Dict{Int, Float64}} = nothing,
)::WorkCounts
    n_sats = max(0, n_active_sats)
    effectors = args.dynamics_model.dynamic_effectors
    total = WorkCounts()

    @inbounds for idx in eachindex(effectors)
        declared = effector_cost_terms(effectors[idx])
        if declared !== nothing
            total = total + declared
            continue
        end
        measured = probe === nothing ? nothing : get(probe, idx, nothing)
        total = total + WorkCounts(
            probe_ns = measured === nothing ? 0.0 : max(0.0, measured),
            unknown_effectors = 1,
        )
    end

    return WorkCounts(
        simd_terms = total.simd_terms,
        scalar_items = total.scalar_items,
        coeff_touches = total.coeff_touches,
        queue_nodes = Float64(n_sats * length(effectors)),
        probe_ns = total.probe_ns,
        unknown_effectors = total.unknown_effectors,
        in_domain = total.in_domain &&
            model_in_cost_domain(args.environment_model.density_model),
    )
end
