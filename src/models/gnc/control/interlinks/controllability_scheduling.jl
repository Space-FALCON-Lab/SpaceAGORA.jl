module ControllabilityScheduling

# Finite-horizon controllability of a constellation whose inputs are pairwise
# link forces, and the schedules of active links that raise it. The dynamics
# enter only as a window of linear time-varying steps (`EdgeLinearStep`), so a
# model in any state representation (Cartesian, orbital elements, KS) serves;
# the feasible-edge graphs (`feasible_edge_graph`) supply which links each step
# admits.

using LinearAlgebra
using ..InterLinkModels: FeasibleEdgeGraph, DutyCycleLimit, terminal_assignment, _pair_set, _duty_window_ok

export EdgeLinearStep, EdgeLinearModel, barycentric_transform, edge_input_matrix, zoh_discretize
export reachability_matrix, controllability_gramian
export AbstractControllabilityMetric, LargestEigenvalueMetric, SmallestEigenvalueMetric, LogDetMetric,
    DirectionEnergyMetric, controllability_metric, gramian_metrics
export AbstractEdgeScheduleStrategy, ConstantEdgeSet, GreedyEdgeSchedule, EdgeSchedule, schedule_edges,
    minimum_input_energy

const EdgePair = Tuple{Int, Int}

# 1. The linear model of a window

"""
    EdgeLinearStep(A, B, edges; columns=nothing)
    EdgeLinearStep(A, blocks::AbstractVector, edges)

One step `k` of a linear time-varying model whose inputs are link forces,
`δx_{k+1} = A_k δx_k + Σ_e B_k^e δu_{k,e}`, in any state representation.

- `A` is the `n × n` state transition of the step.
- `edges` lists the candidate links `(i, j)` (spacecraft vector indices; kept
  as `(min, max)`), and `B` (`n × m`) holds their input columns: the block
  `B[:, columns[e]]` belongs to `edges[e]`. By default each edge has one
  column, its force magnitude, in order.
- The second form takes one block per edge (a vector for one column, or an
  `n × m_e` matrix) and builds `B` and `columns` from them.

The inputs are in the units the columns are scaled to: one newton of link
force for `edge_input_matrix`. A step from another representation (orbital
elements, Kustaanheimo–Stiefel coordinates) needs only its own `A_k` and the
input columns of each candidate edge.
"""
struct EdgeLinearStep
    A::Matrix{Float64}
    B::Matrix{Float64}
    edges::Vector{EdgePair}
    columns::Vector{UnitRange{Int}}

    function EdgeLinearStep(A::AbstractMatrix, B::AbstractMatrix, edges::AbstractVector,
                            columns::AbstractVector{<:AbstractUnitRange})
        n = size(A, 1)
        size(A, 2) == n || throw(DimensionMismatch("A step's A must be square; got $(size(A))."))
        size(B, 1) == n || throw(DimensionMismatch("A step's B has $(size(B, 1)) rows for a state of $n."))
        keys = [minmax(Int(i), Int(j)) for (i, j) in edges]
        all(((i, j),) -> i != j, keys) || throw(ArgumentError("An edge must join two different spacecraft."))
        allunique(keys) || throw(ArgumentError("A step lists an edge twice."))
        length(columns) == length(keys) ||
            throw(ArgumentError("A step has $(length(columns)) column ranges for $(length(keys)) edges."))
        used = falses(size(B, 2))
        for range in columns
            (!isempty(range) && first(range) >= 1 && last(range) <= size(B, 2)) || throw(ArgumentError(
                "Column range $(range) is empty or outside B's $(size(B, 2)) columns."))
            any(used[range]) && throw(ArgumentError("Two edges share a column of B."))
            used[range] .= true
        end
        all(isfinite, A) && all(isfinite, B) || throw(ArgumentError("A step's A and B must be finite."))
        new(Matrix{Float64}(A), Matrix{Float64}(B), keys, [UnitRange{Int}(range) for range in columns])
    end
end

function EdgeLinearStep(A::AbstractMatrix, B::AbstractMatrix, edges::AbstractVector; columns=nothing)
    ranges = columns === nothing ? [e:e for e in eachindex(edges)] : columns
    columns === nothing && size(B, 2) != length(edges) && throw(DimensionMismatch(
        "B has $(size(B, 2)) columns for $(length(edges)) edges; give each edge's column range as `columns`."))
    return EdgeLinearStep(A, B, edges, ranges)
end

function EdgeLinearStep(A::AbstractMatrix, blocks::AbstractVector{<:AbstractArray}, edges::AbstractVector)
    length(blocks) == length(edges) || throw(ArgumentError("A step has $(length(blocks)) input blocks for $(length(edges)) edges."))
    all(block -> ndims(block) <= 2, blocks) || throw(ArgumentError("Each input block is a vector or a matrix."))
    widths = [size(block, 2) for block in blocks]
    stops = cumsum(widths)
    B = isempty(blocks) ? zeros(size(A, 1), 0) : reduce(hcat, [reshape(block, size(block, 1), :) for block in blocks])
    return EdgeLinearStep(A, B, edges, [(stop - width + 1):stop for (stop, width) in zip(stops, widths)])
end

"""
    EdgeLinearModel(steps, projection=I)
    EdgeLinearModel(As, Bs, edges, projection=I)

A window of `N` steps `k₀ … k₀+N−1` of a linear time-varying model with link
inputs (`EdgeLinearStep`), and the map `projection` (`r × n`) from the state at
the window's end, `δx_{k₀+N}`, to the coordinates whose controllability is
assessed. This is the scheduler's representation-agnostic interface: a
sequence of state transitions `A_k` with each candidate edge's input columns
`B_k^e`, and the output map `C` and the barycentric transformation `T_rel`
folded into `projection`.

- **The relative subspace.** Internal link forces cannot move the
  constellation's barycentre, so the useful projection drops it: for a stacked
  Cartesian state, `projection = T[relative_rows, :]` from
  `barycentric_transform`.
- **Physical outputs.** For an output `δy = C δz_rel` (ranges, an RTN
  formation state), `projection = C * T[relative_rows, :]`. A model in other
  coordinates first maps to Cartesian ones with its own output matrix at the
  window's end: `projection = T[relative_rows, :] * C_cart`.
- **Requirements.** Every step lists the same edges in the same order and has
  the state dimension of `projection`'s columns. `As`, `Bs` give one `A_k` and
  one `B_k` (one column per edge) per step.
"""
struct EdgeLinearModel
    steps::Vector{EdgeLinearStep}
    projection::Matrix{Float64}
    edges::Vector{EdgePair}

    function EdgeLinearModel(steps::AbstractVector{EdgeLinearStep}, projection::AbstractMatrix)
        isempty(steps) && throw(ArgumentError("An EdgeLinearModel needs at least one step."))
        n = size(first(steps).A, 1)
        edges = first(steps).edges
        for (k, step) in enumerate(steps)
            size(step.A, 1) == n || throw(DimensionMismatch("Step $k has a state of $(size(step.A, 1)), not $n."))
            step.edges == edges || throw(ArgumentError("Step $k lists other edges than step 1; every step lists the same."))
        end
        size(projection, 2) == n || throw(DimensionMismatch(
            "The projection has $(size(projection, 2)) columns for a state of $n."))
        all(isfinite, projection) || throw(ArgumentError("The projection must be finite."))
        new(collect(steps), Matrix{Float64}(projection), copy(edges))
    end
end

function EdgeLinearModel(steps::AbstractVector{EdgeLinearStep}, projection::UniformScaling=I)
    isempty(steps) && throw(ArgumentError("An EdgeLinearModel needs at least one step."))
    n = size(first(steps).A, 1)
    return EdgeLinearModel(steps, Matrix{Float64}(projection, n, n))
end

function EdgeLinearModel(As::AbstractVector{<:AbstractMatrix}, Bs::AbstractVector{<:AbstractMatrix}, edges::AbstractVector,
                         projection=I)
    length(As) == length(Bs) || throw(ArgumentError("$(length(As)) state matrices for $(length(Bs)) input matrices."))
    return EdgeLinearModel([EdgeLinearStep(A, B, edges) for (A, B) in zip(As, Bs)], projection)
end

"""
    barycentric_transform(masses; block_size=6, relative=:chief) -> (; T, barycenter_rows, relative_rows)

The transformation `T_rel` from a stacked state of `length(masses)` spacecraft
blocks of `block_size` entries each to barycentre and relative coordinates,
`z = T x`. Every entry of a block must add with mass, as Cartesian position and
velocity components do.

- Rows `barycenter_rows` (the first block) are the barycentre,
  `Σ m_i x_i / Σ m_i`.
- Rows `relative_rows` are the relative coordinates of spacecraft `2 … N`:
  `x_i − x_1` with `relative=:chief`, or `x_i − x_cm` with
  `relative=:barycenter`.

`T` is invertible. Masses must be finite and positive.
"""
function barycentric_transform(masses::AbstractVector{<:Real}; block_size::Integer=6, relative::Symbol=:chief)
    relative in (:chief, :barycenter) || throw(ArgumentError("relative is :chief or :barycenter; got $(repr(relative))."))
    block_size >= 1 || throw(ArgumentError("block_size must be positive; got $(block_size)."))
    m = collect(Float64, masses)
    isempty(m) && throw(ArgumentError("barycentric_transform needs at least one mass."))
    all(x -> isfinite(x) && x > 0, m) || throw(ArgumentError("Masses must be finite and positive."))
    count, b = length(m), Int(block_size)
    block(i) = (i - 1) * b + 1:i * b
    Ib = Matrix{Float64}(I, b, b)
    T = zeros(count * b, count * b)
    total = sum(m)
    for i in 1:count
        T[block(1), block(i)] .= (m[i] / total) .* Ib
    end
    for i in 2:count
        rows = block(i)
        if relative === :chief
            T[rows, block(i)] .+= Ib
            T[rows, block(1)] .-= Ib
        else
            T[rows, :] .= -T[block(1), :]
            T[rows, block(i)] .+= Ib
        end
    end
    return (; T, barycenter_rows=block(1), relative_rows=(b + 1):(count * b))
end

"""
    edge_input_matrix(graph::FeasibleEdgeGraph, masses; input_maps=:cartesian, edges=graph.pairs) -> Matrix{Float64}

The continuous input matrix `B_c` (`n × length(edges)`) of a constellation's
link forces at the geometry of `graph`: column `e` is the state rate per newton
of the link `edges[e] = (i, j)`. The link pushes its ends apart along the line
of sight, as the laser links' `link_force`: spacecraft `i` accelerates by
`−ℓ̂_ij / m_i` and `j` by `+ℓ̂_ij / m_j`, with `ℓ̂_ij` the graph's `direction`.

- `input_maps` maps each spacecraft's inertial acceleration into the state
  rate: a vector with one `n × 3` matrix per spacecraft, or `:cartesian` for
  stacked `[r_i; v_i]` blocks of six.
- The force's sign does not change a Gramian (each column enters as `g gᵀ`).
"""
function edge_input_matrix(graph::FeasibleEdgeGraph, masses::AbstractVector{<:Real}; input_maps=:cartesian,
                           edges=graph.pairs)
    count = length(masses)
    count == length(graph.terminal_count) || throw(ArgumentError(
        "$(count) masses for the $(length(graph.terminal_count)) spacecraft of the feasible-edge graph."))
    maps = input_maps === :cartesian ? [_cartesian_input_map(i, count) for i in 1:count] : collect(input_maps)
    length(maps) == count || throw(ArgumentError("$(length(maps)) input maps for $count spacecraft."))
    n = size(first(maps), 1)
    all(rate_map -> size(rate_map) == (n, 3), maps) ||
        throw(DimensionMismatch("Every input map is n × 3 for one state size n."))
    all(x -> isfinite(x) && x > 0, masses) || throw(ArgumentError("Masses must be finite and positive."))
    index = Dict(pair => e for (e, pair) in enumerate(graph.pairs))
    B = zeros(n, length(edges))
    for (column, (i, j)) in enumerate(edges)
        e = get(index, minmax(i, j), 0)
        e == 0 && throw(ArgumentError("Edge ($i, $j) is not a pair of the feasible-edge graph."))
        first_end, second_end = graph.pairs[e]
        sight = graph.direction[e]
        B[:, column] .= maps[second_end] * (sight / masses[second_end]) .- maps[first_end] * (sight / masses[first_end])
    end
    return B
end

function _cartesian_input_map(i::Int, count::Int)
    rate_map = zeros(6 * count, 3)
    rate_map[6 * (i - 1) .+ (4:6), :] .= Matrix{Float64}(I, 3, 3)
    return rate_map
end

"""
    zoh_discretize(Ac, Bc, dt) -> (A, B)

The zero-order-hold discretization of `ẋ = Ac x + Bc u` over a step of `dt`:
`A = exp(Ac dt)`, `B = ∫₀^dt exp(Ac τ) dτ Bc`, both from the exponential of
the augmented matrix `[Ac Bc; 0 0] dt`.
"""
function zoh_discretize(Ac::AbstractMatrix, Bc::AbstractVecOrMat, dt::Real)
    n = size(Ac, 1)
    size(Ac, 2) == n || throw(DimensionMismatch("Ac must be square; got $(size(Ac))."))
    size(Bc, 1) == n || throw(DimensionMismatch("Bc has $(size(Bc, 1)) rows for a state of $n."))
    isfinite(dt) && dt > 0 || throw(ArgumentError("dt must be finite and positive; got $(dt)."))
    B = reshape(Bc, n, :)
    m = size(B, 2)
    augmented = zeros(n + m, n + m)
    augmented[1:n, 1:n] .= Ac
    augmented[1:n, n+1:end] .= B
    E = exp(Float64(dt) .* augmented)
    return E[1:n, 1:n], E[1:n, n+1:end]
end

# 2. Reachability and the Gramian

# Each step's per-edge contributions to the window's projected end state,
# `P Φ(k₀+N, k+1) B_k^e`, and the free response `P Φ(k₀+N, k₀)`.
function _reachability_blocks(model::EdgeLinearModel)
    S = copy(model.projection)
    blocks = Vector{Vector{Matrix{Float64}}}(undef, length(model.steps))
    for k in reverse(eachindex(model.steps))
        step = model.steps[k]
        blocks[k] = [S * step.B[:, range] for range in step.columns]
        S = S * step.A
    end
    return blocks, S
end

# The active edge sets per step as model edge indices, checked.
function _active_indices(model::EdgeLinearModel, active)
    index = Dict(pair => e for (e, pair) in enumerate(model.edges))
    active === nothing && return [collect(eachindex(model.edges)) for _ in model.steps]
    length(active) == length(model.steps) || throw(ArgumentError(
        "The schedule has $(length(active)) steps; the model has $(length(model.steps))."))
    return map(enumerate(active)) do (k, set)
        sort!(map(collect(_pair_set(set))) do pair
            e = get(index, pair, 0)
            e == 0 && throw(ArgumentError("Step $k activates $pair, which is not an edge of the model."))
            e
        end)
    end
end

"""
    reachability_matrix(model::EdgeLinearModel; active=nothing) -> (; matrix, columns)

The projected reachability matrix of the window, `R_N(k₀) = P [Φ(k₀+N, k₀+1)
B_{k₀} … B_{k₀+N−1}]` with `Φ(k, j) = A_{k−1} ⋯ A_j`, restricted to the
active links: `active[k]` lists the edges active at step `k` (all edges at
every step when `nothing`). `columns[c]` names column `c` as `(step, edge,
input)`, the input being the column's place in the edge's block.

With the projection `P` a physical output map, column `(k, e)` is the change
of that output at the window's end per unit input on edge `e` at step `k`: the
terminal gain of that activation.
"""
function reachability_matrix(model::EdgeLinearModel; active=nothing)
    blocks, _ = _reachability_blocks(model)
    chosen = _active_indices(model, active)
    columns = Tuple{Int, EdgePair, Int}[]
    parts = Matrix{Float64}[]
    for (k, edges) in enumerate(chosen), e in edges
        push!(parts, blocks[k][e])
        append!(columns, [(k, model.edges[e], c) for c in 1:size(blocks[k][e], 2)])
    end
    matrix = isempty(parts) ? zeros(size(model.projection, 1), 0) : reduce(hcat, parts)
    return (; matrix, columns)
end

"""
    controllability_gramian(model::EdgeLinearModel; active=nothing) -> Matrix{Float64}

The finite-horizon controllability Gramian of the window, `W_N(k₀) = R_N R_Nᵀ`
(`reachability_matrix`), in the projected coordinates: with the relative
projection, the relative-subspace Gramian; with an output map, the output
Gramian `W_N^y`. `active` as for `reachability_matrix`.
"""
function controllability_gramian(model::EdgeLinearModel; active=nothing)
    R = reachability_matrix(model; active).matrix
    return Matrix(Symmetric(R * R'))
end

# 3. Metrics

"""
    AbstractControllabilityMetric

Supertype of the scalar measures `Ψ(W)` of a controllability Gramian that the
edge schedules maximize: `LargestEigenvalueMetric`, `SmallestEigenvalueMetric`,
`LogDetMetric` and `DirectionEnergyMetric`. A metric implements
`controllability_metric`. Each built-in metric is nondecreasing as links are
added, since every activation adds a positive semidefinite term to `W`.
"""
abstract type AbstractControllabilityMetric end

"""
    LargestEigenvalueMetric()

`λ₁(W)`: the input energy reach of the easiest direction, whose eigenvector is
the dominant deformation mode.
"""
struct LargestEigenvalueMetric <: AbstractControllabilityMetric end

"""
    SmallestEigenvalueMetric()

`λ_min(W)`: the reach of the hardest direction. It stays zero until the active
links make every projected direction reachable.
"""
struct SmallestEigenvalueMetric <: AbstractControllabilityMetric end

"""
    LogDetMetric(; epsilon)

`log det(W + ε I)`: the log volume of the reachable ellipsoid, with the
regularizer `ε > 0` (finite, in the units of `W`).
"""
struct LogDetMetric <: AbstractControllabilityMetric
    epsilon::Float64

    function LogDetMetric(epsilon::Real)
        isfinite(epsilon) && epsilon > 0 || throw(ArgumentError("epsilon must be finite and positive; got $(epsilon)."))
        new(Float64(epsilon))
    end
end

LogDetMetric(; epsilon::Real) = LogDetMetric(epsilon)

"""
    DirectionEnergyMetric(direction)

`dᵀ W d`: the reach along a mission-chosen deformation `d` (a finite, nonzero
vector in the projected coordinates, used as given).
"""
struct DirectionEnergyMetric <: AbstractControllabilityMetric
    direction::Vector{Float64}

    function DirectionEnergyMetric(direction::AbstractVector{<:Real})
        d = collect(Float64, direction)
        all(isfinite, d) && any(!iszero, d) || throw(ArgumentError("The direction must be finite and nonzero."))
        new(d)
    end
end

"""
    controllability_metric(metric::AbstractControllabilityMetric, W) -> Float64

The value of `metric` for the symmetric positive semidefinite Gramian `W`.
"""
function controllability_metric end

controllability_metric(::LargestEigenvalueMetric, W::AbstractMatrix) = last(eigvals(Symmetric(Matrix{Float64}(W))))
controllability_metric(::SmallestEigenvalueMetric, W::AbstractMatrix) = first(eigvals(Symmetric(Matrix{Float64}(W))))
# From the eigenvalues, with the roundoff below zero of a semidefinite `W` taken as zero.
controllability_metric(metric::LogDetMetric, W::AbstractMatrix) =
    sum(value -> log(max(value, 0.0) + metric.epsilon), eigvals(Symmetric(Matrix{Float64}(W))))

function controllability_metric(metric::DirectionEnergyMetric, W::AbstractMatrix)
    length(metric.direction) == size(W, 1) || throw(DimensionMismatch(
        "The direction has $(length(metric.direction)) entries for a Gramian of $(size(W, 1))."))
    return dot(metric.direction, Symmetric(Matrix{Float64}(W)), metric.direction)
end

"""
    gramian_metrics(W; epsilon, direction=nothing) -> NamedTuple

Every metric of the Gramian `W` at once, with its eigenstructure:
`largest_eigenvalue` (`λ₁`), `smallest_eigenvalue` (`λ_min`), `log_det`
(`log det(W + εI)` with `epsilon`), `direction_energy` (`dᵀWd`, `NaN` without a
`direction`), `eigenvalues` (ascending), `modes` (the matching unit
eigenvectors, by column) and `dominant_mode` (`v₁`, the eigenvector of `λ₁`).
"""
function gramian_metrics(W::AbstractMatrix; epsilon::Real, direction=nothing)
    decomposition = eigen(Symmetric(Matrix{Float64}(W)))
    values, modes = decomposition.values, decomposition.vectors
    energy = direction === nothing ? NaN : controllability_metric(DirectionEnergyMetric(direction), W)
    return (; largest_eigenvalue=last(values), smallest_eigenvalue=first(values),
        log_det=controllability_metric(LogDetMetric(epsilon), W), direction_energy=energy, eigenvalues=values,
        modes, dominant_mode=modes[:, end])
end

# The difference below which two values of `metric` count as equal: a relative
# share of the Gramian's mean eigenvalue (or, for the log volume, of one neper).
_tie_tolerance(::LogDetMetric, W) = 1e-10
_tie_tolerance(metric::DirectionEnergyMetric, W) = 1e-10 * dot(metric.direction, metric.direction) * tr(W) / size(W, 1)
_tie_tolerance(::AbstractControllabilityMetric, W) = 1e-10 * tr(W) / size(W, 1)

# 4. Schedules

"""
    AbstractEdgeScheduleStrategy

Supertype of the ways `schedule_edges` picks the active links of a window:
`ConstantEdgeSet` and `GreedyEdgeSchedule`.
"""
abstract type AbstractEdgeScheduleStrategy end

"""
    ConstantEdgeSet(; edges=nothing, metric=nothing, tie_break=nothing)

One link set kept over the whole window: at each step, its edges that are
admissible there (feasible, a free terminal at each end, within the duty
cycle), in the set's order.

- With `edges`, that set (at most `max_active` edges).
- Otherwise the set is built from the candidate edges by `metric`: each round
  adds the edge whose constant use gives the largest metric, until `max_active`
  edges; ties (within the tolerance of `GreedyEdgeSchedule`) go to
  `tie_break`, a second metric, and then to the earlier edge.
"""
struct ConstantEdgeSet{M, T} <: AbstractEdgeScheduleStrategy
    edges::Union{Nothing, Vector{EdgePair}}
    metric::M
    tie_break::T

    function ConstantEdgeSet(edges, metric::M, tie_break::T) where {M, T}
        (edges === nothing) == (metric === nothing) && throw(ArgumentError(
            "A ConstantEdgeSet takes its edges or a metric to choose them, one of the two."))
        (M <: Union{Nothing, AbstractControllabilityMetric} && T <: Union{Nothing, AbstractControllabilityMetric}) ||
            throw(ArgumentError("metric and tie_break are controllability metrics or nothing."))
        set = edges === nothing ? nothing : [minmax(Int(i), Int(j)) for (i, j) in edges]
        set === nothing || allunique(set) || throw(ArgumentError("A ConstantEdgeSet lists an edge twice."))
        new{M, T}(set, metric, tie_break)
    end
end

ConstantEdgeSet(; edges=nothing, metric=nothing, tie_break=nothing) = ConstantEdgeSet(edges, metric, tie_break)

"""
    GreedyEdgeSchedule(metric; tie_break=nothing)

A schedule built one activation at a time over the whole window: each round
adds the `(step, edge)` activation, among those still admissible, that gives
the largest `metric` of the Gramian; ties (within a relative `1e-10` of the
Gramian's mean eigenvalue, or `1e-10` for `LogDetMetric`) go to `tie_break`, a
second metric, and then to the earlier step and edge. It fills every admissible
slot: the built-in metrics never fall as links are added, so the metric decides
which links are active, not how many.

- `GreedyEdgeSchedule(LargestEigenvalueMetric())` is the greedy `λ₁` schedule.
- `GreedyEdgeSchedule(LogDetMetric(; epsilon))` is the greedy volume schedule.
  `log det(W + εI)` is monotone and submodular in the activations, the property
  greedy approximation bounds rest on; no bound is claimed here under terminal
  and duty-cycle limits.
- `GreedyEdgeSchedule(SmallestEigenvalueMetric(); tie_break=LogDetMetric(; epsilon))`
  is the greedy `λ_min` schedule: `λ_min` stays zero until every projected
  direction is reachable, so the tie break chooses until then.
"""
struct GreedyEdgeSchedule{M <: AbstractControllabilityMetric, T} <: AbstractEdgeScheduleStrategy
    metric::M
    tie_break::T

    function GreedyEdgeSchedule(metric::M, tie_break::T) where {M <: AbstractControllabilityMetric, T}
        T <: Union{Nothing, AbstractControllabilityMetric} ||
            throw(ArgumentError("tie_break is a controllability metric or nothing."))
        new{M, T}(metric, tie_break)
    end
end

GreedyEdgeSchedule(metric::AbstractControllabilityMetric; tie_break=nothing) = GreedyEdgeSchedule(metric, tie_break)

"""
    EdgeSchedule

The active links of a window, from `schedule_edges`.

- `active[k]`: the edges active at step `k`, in the model's edge order (the
  `active` argument of `reachability_matrix` and `controllability_gramian`).
- `gramian`: the window's projected Gramian under this schedule.
- `terminals[k]`: with feasible-edge graphs, each active edge's terminals
  `(at i, at j)` at step `k` (`terminal_assignment`); empty otherwise.
"""
struct EdgeSchedule
    active::Vector{Vector{EdgePair}}
    gramian::Matrix{Float64}
    terminals::Vector{Dict{EdgePair, NTuple{2, Int}}}
end

# What a schedule may do at each step: the window's contributions, which edges
# each step admits, and the limits on a step's set.
struct _SchedulePlan
    rows::Int
    edges::Vector{EdgePair}
    blocks::Vector{Vector{Matrix{Float64}}}
    allowed::Matrix{Bool}
    graphs::Union{Nothing, Vector{FeasibleEdgeGraph}}
    max_active::Int
    duty::Union{Nothing, DutyCycleLimit}
    previous::Vector{Set{EdgePair}}
end

"""
    schedule_edges(model::EdgeLinearModel, strategy; max_active, graphs=nothing, feasible=nothing,
                   duty_cycle=nothing, previous=()) -> EdgeSchedule

The active links `σ_k` of each step of the window `model`, chosen by `strategy`
(`ConstantEdgeSet`, `GreedyEdgeSchedule`) to make the Gramian large, subject
to `σ_k ⊆ E_k` and `|σ_k| ≤ max_active`.

- `graphs` gives the feasible-edge graph `G_k` of each step
  (`feasible_edge_graph`): its line of sight, range and field-of-regard flags
  give `E_k`, and a step's set must have a `terminal_assignment`. Each model
  edge must be a pair of every graph.
- Or `feasible`, an `edges × steps` Boolean matrix, gives `E_k` alone, with no
  terminal limit; with neither, every edge is feasible at every step.
- `duty_cycle` (by default the graphs' `DutyCycleLimit`) limits the schedule
  over the window, counting the steps `previous` (active sets, oldest first)
  just before it.
"""
function schedule_edges(model::EdgeLinearModel, strategy::AbstractEdgeScheduleStrategy; max_active::Integer,
                        graphs=nothing, feasible=nothing, duty_cycle=nothing, previous=())
    plan = _schedule_plan(model, max_active, graphs, feasible, duty_cycle, previous)
    active = _choose(strategy, plan)
    ordered = [sort(set; by=edge -> findfirst(==(edge), plan.edges)) for set in active]
    terminals = graphs === nothing ? [Dict{EdgePair, NTuple{2, Int}}() for _ in active] :
        [Dict{EdgePair, NTuple{2, Int}}(terminal_assignment(graphs[k], set)) for (k, set) in enumerate(ordered)]
    return EdgeSchedule(ordered, _gramian(plan, ordered), terminals)
end

function _schedule_plan(model::EdgeLinearModel, max_active::Integer, graphs, feasible, duty_cycle, previous)
    max_active >= 0 || throw(ArgumentError("max_active must be nonnegative; got $(max_active)."))
    graphs === nothing || feasible === nothing ||
        throw(ArgumentError("Give the feasible-edge graphs or a feasibility matrix, not both."))
    N, edges = length(model.steps), model.edges
    allowed = fill(true, length(edges), N)
    if graphs !== nothing
        length(graphs) == N || throw(ArgumentError("$(length(graphs)) feasible-edge graphs for $N steps."))
        for (k, graph) in enumerate(graphs)
            index = Dict(pair => e for (e, pair) in enumerate(graph.pairs))
            for (e, edge) in enumerate(edges)
                g = get(index, edge, 0)
                g == 0 && throw(ArgumentError("Edge $edge of the model is not a pair of the feasible-edge graph of step $k."))
                allowed[e, k] = graph.line_of_sight[g] && graph.in_range[g] && all(graph.field_of_regard[g])
            end
        end
    elseif feasible !== nothing
        size(feasible) == (length(edges), N) || throw(DimensionMismatch(
            "The feasibility matrix is $(size(feasible)); the model has $(length(edges)) edges and $N steps."))
        allowed .= feasible
    end
    duty = duty_cycle !== nothing ? duty_cycle : graphs === nothing ? nothing : first(graphs).duty_limit
    blocks, _ = _reachability_blocks(model)
    return _SchedulePlan(size(model.projection, 1), copy(edges), blocks, allowed, graphs === nothing ? nothing : collect(graphs), Int(max_active),
        duty, [_pair_set(set) for set in previous])
end

# Whether edge `e` can join the edges `active[k]` already active at step `k`.
function _admits(plan::_SchedulePlan, active::Vector{Vector{EdgePair}}, k::Int, e::Int)
    edge = plan.edges[e]
    plan.allowed[e, k] && !(edge in active[k]) && length(active[k]) < plan.max_active || return false
    candidate = vcat(active[k], [edge])
    plan.graphs === nothing || terminal_assignment(plan.graphs[k], candidate) !== nothing || return false
    plan.duty === nothing && return true
    steps = vcat(plan.previous, [Set(set) for set in active])
    steps[length(plan.previous) + k] = Set(candidate)
    return _duty_window_ok(plan.duty, steps, length(plan.previous) + k)
end

function _gramian(plan::_SchedulePlan, active)
    W = zeros(plan.rows, plan.rows)
    index = Dict(edge => e for (e, edge) in enumerate(plan.edges))
    for (k, set) in enumerate(active), edge in set
        G = plan.blocks[k][index[edge]]
        W .+= G * G'
    end
    return Matrix(Symmetric(W))
end

# The value of `metric` and `tie_break` for `W`, with the tolerance of a tie.
function _score(metric, tie_break, W)
    value = controllability_metric(metric, W)
    second = tie_break === nothing ? 0.0 : controllability_metric(tie_break, W)
    return (value, second, _tie_tolerance(metric, W))
end

_beats(score, best) = best === nothing || score[1] > best[1] + max(score[3], best[3]) ||
    (abs(score[1] - best[1]) <= max(score[3], best[3]) && score[2] > best[2])

function _choose(strategy::GreedyEdgeSchedule, plan::_SchedulePlan)
    N = length(plan.blocks)
    active = [EdgePair[] for _ in 1:N]
    W = zeros(plan.rows, plan.rows)
    while true
        best, choice = nothing, (0, 0)
        for k in 1:N, e in eachindex(plan.edges)
            _admits(plan, active, k, e) || continue
            G = plan.blocks[k][e]
            score = _score(strategy.metric, strategy.tie_break, Symmetric(W + G * G'))
            _beats(score, best) && ((best, choice) = (score, (k, e)))
        end
        best === nothing && return active
        k, e = choice
        push!(active[k], plan.edges[e])
        W .+= plan.blocks[k][e] * plan.blocks[k][e]'
    end
end

function _choose(strategy::ConstantEdgeSet, plan::_SchedulePlan)
    if strategy.edges !== nothing
        length(strategy.edges) <= plan.max_active || throw(ArgumentError(
            "The constant set has $(length(strategy.edges)) edges; max_active is $(plan.max_active)."))
        for edge in strategy.edges
            edge in plan.edges || throw(ArgumentError("Edge $edge of the constant set is not an edge of the model."))
        end
        return _constant_active(plan, strategy.edges)
    end
    set = EdgePair[]
    while length(set) < plan.max_active
        best, choice = nothing, nothing
        for edge in plan.edges
            edge in set && continue
            score = _score(strategy.metric, strategy.tie_break, _gramian(plan, _constant_active(plan, vcat(set, [edge]))))
            _beats(score, best) && ((best, choice) = (score, edge))
        end
        choice === nothing && break
        push!(set, choice)
    end
    return _constant_active(plan, set)
end

# The active sets of a constant edge set: at each step in turn, its admissible
# edges in the set's order.
function _constant_active(plan::_SchedulePlan, set::Vector{EdgePair})
    active = [EdgePair[] for _ in plan.blocks]
    index = Dict(edge => e for (e, edge) in enumerate(plan.edges))
    for k in eachindex(active), edge in set
        _admits(plan, active, k, index[edge]) && push!(active[k], edge)
    end
    return active
end

# 5. Input energy

"""
    minimum_input_energy(model::EdgeLinearModel, schedule, target; initial_state=nothing, rtol=1e-10, reach_tol=1e-8)
        -> (; energy, inputs, residual, reachable)

The least input energy `Σ |δu|²` that moves the window's projected end state
to `target` with the links of `schedule` (an `EdgeSchedule`, or active sets as
for `reachability_matrix`), from `δx_{k₀} = 0` or from `initial_state`, whose
free response the inputs then correct.

- `energy` is `yᵀ W⁺ y`, the minimum-norm solution of `R u = y` for the
  required change `y`, with singular values below `rtol` times the largest
  dropped.
- `inputs[k]` maps each edge active at step `k` to its input vector.
- `residual` is `|R u − y| / |y|` (zero for `y = 0`), and `reachable` whether it
  is at most `reach_tol`: a target outside the reachable subspace gets the
  least-squares inputs, `reachable = false`.
"""
function minimum_input_energy(model::EdgeLinearModel, schedule, target::AbstractVector{<:Real}; initial_state=nothing,
                              rtol::Real=1e-10, reach_tol::Real=1e-8)
    active = schedule isa EdgeSchedule ? schedule.active : schedule
    length(target) == size(model.projection, 1) || throw(DimensionMismatch(
        "The target has $(length(target)) entries for $(size(model.projection, 1)) projected coordinates."))
    reach = reachability_matrix(model; active)
    required = collect(Float64, target)
    if initial_state !== nothing
        _, free = _reachability_blocks(model)
        length(initial_state) == size(free, 2) || throw(DimensionMismatch(
            "The initial state has $(length(initial_state)) entries for a state of $(size(free, 2))."))
        required .-= free * collect(Float64, initial_state)
    end
    R = reach.matrix
    u = isempty(R) ? zeros(size(R, 2)) : pinv(R; rtol=Float64(rtol)) * required
    scale = norm(required)
    residual = scale == 0 ? 0.0 : norm(R * u - required) / scale
    inputs = [Dict{EdgePair, Vector{Float64}}() for _ in model.steps]
    for (c, (k, edge, _)) in enumerate(reach.columns)
        push!(get!(inputs[k], edge, Float64[]), u[c])
    end
    return (; energy=dot(u, u), inputs, residual, reachable=residual <= reach_tol)
end

end
