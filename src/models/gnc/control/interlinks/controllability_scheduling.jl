module ControllabilityScheduling

# Finite-horizon controllability of a constellation whose inputs are pairwise
# link forces, and the schedules of active links that raise it. The dynamics
# enter only as a window of linear time-varying steps (`EdgeLinearStep`), so a
# model in any state representation (Cartesian, orbital elements, KS) serves;
# the feasible-edge graphs (`feasible_edge_graph`) supply which links each step
# admits.

using LinearAlgebra
using StaticArrays: SVector
using ..EffectorSampling: StateSample
using ..InterLinkModels: FeasibleEdgeGraph, DutyCycleLimit, terminal_assignment, _pair_set, _duty_window_ok

export EdgeLinearStep, EdgeLinearModel, barycentric_transform, edge_input_matrix, zoh_discretize
export link_force_jacobian, range_output_matrix
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
inputs (`EdgeLinearStep`), and the map `projection` (`r × n`, `r ≥ 1`) from the
state at the window's end, `δx_{k₀+N}`, to the coordinates whose
controllability is assessed. This is the scheduler's representation-agnostic
interface: a sequence of state transitions `A_k` with each candidate edge's
input columns `B_k^e`, and the output map `C` and the barycentric
transformation `T_rel` folded into `projection`.

- **The relative subspace.** The link forces have no barycenter component at
  any instant (the continuous `B_c`), so for a stacked Cartesian state
  `projection = T[relative_rows, :]` from `barycentric_transform` drops a mode
  they cannot drive directly. Over a window in orbit, unequal gravity couples
  the barycenter weakly to the relative motion, so a discretized `B_k` has
  small, nonzero barycenter rows.
- **Rigid rotation.** The relative rows still hold the formation's rigid
  rotation, which line-of-centers forces cannot drive in free space and drive
  only weakly in orbit (they leave the total angular momentum about the
  planet's center unchanged). For `λ_min` and `log det`, prefer a
  rotation-free output such as `range_output_matrix`.
- **Units.** Rows in different units (metres and metres per second) spread the
  Gramian's eigenvalues by the units alone; scale them to one unit, for
  example the velocity or rate rows by `1 / n` for the mean motion `n`.
- **Physical outputs.** For an output `δy = C δz_rel` (ranges, an RTN
  formation state), `projection = C * T[relative_rows, :]`. A model in other
  coordinates first maps to Cartesian ones with its own output matrix at the
  window's end: `projection = T[relative_rows, :] * C_cart`.
- **Requirements.** Every step lists the same edges in the same order and has
  the state dimension of `projection`'s columns. `As`, `Bs` give one `A_k` and
  one `B_k` (one column per edge) per step. In a model whose spacecraft keep
  their own time variables (Kustaanheimo–Stiefel fictitious time), each step
  must span one common physical-time interval for every spacecraft, so that
  both ends of a link are forced together.
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
        size(projection, 1) >= 1 || throw(ArgumentError("The projection needs at least one row."))
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
blocks of `block_size` entries each to barycenter and relative coordinates,
`z = T x`. Every entry of a block must add with mass, as Cartesian position and
velocity components do.

- Rows `barycenter_rows` (the first block) are the barycenter,
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
of the link `edges[e] = (i, j)`.

- **Sign.** A positive input pushes the two ends apart along the line of
  sight, as the laser links' `link_force` does: spacecraft `i` accelerates by
  `−ℓ̂_ij / m_i` and `j` by `+ℓ̂_ij / m_j`, with `ℓ̂_ij` the graph's
  `direction`. A negative input pulls them together, which a laser cannot do:
  a model linearized about zero link force can use only nonnegative inputs.
  Giving an edge as `(j, i)` does not change its column.
- **Input maps.** `input_maps` maps each spacecraft's inertial acceleration
  into the state rate: a vector with one `n × 3` matrix per spacecraft, or
  `:cartesian` for stacked `[r_i; v_i]` blocks of six.
- The sign does not change a Gramian (each column enters as `g gᵀ`), but it
  does change the inputs of `minimum_input_energy`.
"""
function edge_input_matrix(graph::FeasibleEdgeGraph, masses::AbstractVector{<:Real}; input_maps=:cartesian,
                           edges=graph.pairs)
    count = _check_masses(graph, masses)
    maps = input_maps === :cartesian ? [_cartesian_input_map(i, count) for i in 1:count] : collect(input_maps)
    length(maps) == count || throw(ArgumentError("$(length(maps)) input maps for $count spacecraft."))
    n = size(first(maps), 1)
    all(rate_map -> size(rate_map) == (n, 3), maps) ||
        throw(DimensionMismatch("Every input map is n × 3 for one state size n."))
    B = zeros(n, length(edges))
    for (column, edge) in enumerate(edges)
        e = _graph_pair(graph, edge)
        first_end, second_end = graph.pairs[e]
        sight = graph.direction[e]
        B[:, column] .= maps[second_end] * (sight / masses[second_end]) .- maps[first_end] * (sight / masses[first_end])
    end
    return B
end

"""
    link_force_jacobian(graph::FeasibleEdgeGraph, masses, forces; edges=graph.pairs) -> Matrix{Float64}

The state Jacobian, for stacked Cartesian states `[r_i; v_i]`, of nominal link
forces `forces[e] ≥ 0` (N) held on `edges` at the geometry of `graph`. A model
linearized about a positive nominal force adds it to its continuous `Ac`;
the signed inputs about that force are then perturbations within
`−F* ≤ δu ≤ F_max − F*`.

Each force pushes its ends apart along `ℓ̂_ij`, which turns with the relative
position: the force on `j` changes by `(F*/ρ)(I − ℓ̂ℓ̂ᵀ)(δr_j − δr_i)` and the
force on `i` by its negative. The force's magnitude is held fixed; a beam force
that falls with range adds a radial term that is not included.
"""
function link_force_jacobian(graph::FeasibleEdgeGraph, masses::AbstractVector{<:Real}, forces::AbstractVector{<:Real};
                             edges=graph.pairs)
    count = _check_masses(graph, masses)
    length(forces) == length(edges) || throw(ArgumentError("$(length(forces)) forces for $(length(edges)) edges."))
    all(f -> isfinite(f) && f >= 0, forces) || throw(ArgumentError("The nominal forces must be finite and nonnegative."))
    J = zeros(6 * count, 6 * count)
    position(s) = 6 * (s - 1) .+ (1:3)
    velocity(s) = 6 * (s - 1) .+ (4:6)
    for (force, edge) in zip(forces, edges)
        force == 0 && continue
        e = _graph_pair(graph, edge)
        i, j = graph.pairs[e]
        distance = graph.range_m[e]
        distance > 0 || throw(ArgumentError("Edge $(graph.pairs[e]) has no line of centers: its ends coincide."))
        sight = graph.direction[e]
        stiffness = (force / distance) .* (Matrix{Float64}(I, 3, 3) - sight * sight')
        J[velocity(j), position(j)] .+= stiffness ./ masses[j]
        J[velocity(j), position(i)] .-= stiffness ./ masses[j]
        J[velocity(i), position(i)] .+= stiffness ./ masses[i]
        J[velocity(i), position(j)] .-= stiffness ./ masses[i]
    end
    return J
end

"""
    range_output_matrix(samples, pairs; rates=true, rate_scale=1.0) -> Matrix{Float64}
    range_output_matrix(state, pairs; rates=true, rate_scale=1.0)

An output map for stacked Cartesian states `[r_i; v_i]`: the first-order change
of the range `ρ_ij` of each pair in `pairs`, then, with `rates`, of each range
rate `ρ̇_ij` times `rate_scale`, at the given states (one `StateSample` per
spacecraft, or a run state with spacecraft `i` at `state.sc[i]`).

These outputs ignore the formation's rigid motion: a translation, or a
rotation of every position and velocity together, changes none of them. Used as
the `projection` of an `EdgeLinearModel`, they keep `λ_min` and `log det` on the
formation's shape. `rate_scale = 1 / n`, for the mean motion `n`, puts the
rates in metres like the ranges.
"""
function range_output_matrix(samples::AbstractVector{<:StateSample}, pairs; rates::Bool=true, rate_scale::Real=1.0)
    return _range_output(samples, length(samples), pairs, rates, Float64(rate_scale))
end

function range_output_matrix(state, pairs; rates::Bool=true, rate_scale::Real=1.0)
    count = length(state.sc)
    samples = [(pos_ii=SVector{3, Float64}(state.sc[i].pos), vel_ii=SVector{3, Float64}(state.sc[i].vel)) for i in 1:count]
    return _range_output(samples, count, pairs, rates, Float64(rate_scale))
end

function _range_output(samples, count::Int, pairs, rates::Bool, rate_scale::Float64)
    isfinite(rate_scale) && rate_scale > 0 || throw(ArgumentError("rate_scale must be finite and positive; got $(rate_scale)."))
    keys = [minmax(Int(i), Int(j)) for (i, j) in pairs]
    isempty(keys) && throw(ArgumentError("range_output_matrix needs at least one pair."))
    allunique(keys) || throw(ArgumentError("A pair is listed twice."))
    C = zeros(rates ? 2 * length(keys) : length(keys), 6 * count)
    for (row, (i, j)) in enumerate(keys)
        (1 <= i && j <= count && i != j) || throw(ArgumentError("Pair ($i, $j) names a spacecraft outside 1:$count."))
        separation = samples[j].pos_ii - samples[i].pos_ii
        distance = norm(separation)
        distance > 0 || throw(ArgumentError("Pair ($i, $j) has no range: its ends coincide."))
        sight = separation / distance
        position(s) = 6 * (s - 1) .+ (1:3)
        velocity(s) = 6 * (s - 1) .+ (4:6)
        C[row, position(j)] .= sight
        C[row, position(i)] .= -sight
        rates || continue
        closing = samples[j].vel_ii - samples[i].vel_ii
        turning = (closing - dot(sight, closing) * sight) / distance
        rate_row = length(keys) + row
        C[rate_row, velocity(j)] .= rate_scale .* sight
        C[rate_row, velocity(i)] .= -rate_scale .* sight
        C[rate_row, position(j)] .= rate_scale .* turning
        C[rate_row, position(i)] .= -rate_scale .* turning
    end
    return C
end

# The number of spacecraft of `graph`, checked against `masses`.
function _check_masses(graph::FeasibleEdgeGraph, masses)
    count = length(masses)
    count == length(graph.terminal_count) || throw(ArgumentError(
        "$(count) masses for the $(length(graph.terminal_count)) spacecraft of the feasible-edge graph."))
    all(x -> isfinite(x) && x > 0, masses) || throw(ArgumentError("Masses must be finite and positive."))
    return count
end

# The index of `edge` among `graph`'s pairs.
function _graph_pair(graph::FeasibleEdgeGraph, (i, j))
    e = findfirst(==(minmax(Int(i), Int(j))), graph.pairs)
    e === nothing && throw(ArgumentError("Edge ($i, $j) is not a pair of the feasible-edge graph."))
    return e
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
added, since every activation adds a positive semidefinite term to `W`. The
schedules count two values as equal within the roundoff of a symmetric
eigensolver (`eps · n · λ_max`, with the trace for `λ_max`), scaled to the
metric; a metric of one's own gets that eigenvalue tolerance.
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
links make every projected direction reachable, so it needs a projection on
which that can happen: one in common units, without the barycenter or the
rigid rotation (`range_output_matrix`).
"""
struct SmallestEigenvalueMetric <: AbstractControllabilityMetric end

"""
    LogDetMetric(; epsilon)

`log det(W + ε I)`: the log volume of the reachable ellipsoid, with the
regularizer `ε > 0` (finite, in the units of `W`). Directions whose eigenvalue
lies below `ε` count little, so choose `ε` against the spectrum of a projection
in common units.
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

# The difference below which two values of a metric count as equal, for a
# Gramian of size `n` and trace `trace`: the roundoff of a symmetric
# eigensolver, about eps · n · λ_max (the trace bounds λ_max), carried into the
# metric, plus the rounding of the value itself.
const _TIE_ULPS = 16
_tie_tolerance(::Nothing, n, trace, value) = 0.0
_tie_tolerance(::AbstractControllabilityMetric, n, trace, value) =
    _TIE_ULPS * eps() * (n * max(trace, 0.0) + abs(value))
_tie_tolerance(metric::DirectionEnergyMetric, n, trace, value) =
    _TIE_ULPS * eps() * (n * dot(metric.direction, metric.direction) * max(trace, 0.0) + abs(value))
_tie_tolerance(metric::LogDetMetric, n, trace, value) =
    _TIE_ULPS * eps() * (n * (max(trace, 0.0) / metric.epsilon + 1.0) + abs(value))

# The metric of `W + G Gᵀ` for the candidate blocks `G` of one greedy round:
# from a factorization of `W + εI` for the log volume, from `dᵀWd` for a
# direction, and by eigenvalues otherwise.
_round(::Nothing, W) = nothing
_round(metric::AbstractControllabilityMetric, W) = (metric, W)
_evaluate(::Nothing, G) = 0.0
_evaluate((metric, W)::Tuple, G) = controllability_metric(metric, Symmetric(W + G * G'))

struct _LogDetRound
    base::Float64
    factor::Cholesky{Float64, Matrix{Float64}}
end

function _round(metric::LogDetMetric, W)
    factor = cholesky(Symmetric(W + metric.epsilon * I); check=false)
    issuccess(factor) || return (metric, W)
    return _LogDetRound(controllability_metric(metric, W), factor)
end

function _evaluate(round::_LogDetRound, G)
    H = round.factor.L \ G
    return round.base + logdet(cholesky(Symmetric(Matrix{Float64}(I, size(H, 2), size(H, 2)) + H' * H)))
end

_round(metric::DirectionEnergyMetric, W) = (controllability_metric(metric, W), metric.direction)
_evaluate((base, direction)::Tuple{Float64, Vector{Float64}}, G) = base + sum(abs2, G' * direction)

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
  edges. Ties (within the roundoff tolerance of `GreedyEdgeSchedule`) go to
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
the largest `metric` of the Gramian.

- **Ties.** Two values within the roundoff of a symmetric eigensolver (about
  `eps · n · λ_max`, carried into the metric) are equal. A tie goes to
  `tie_break`, a second metric compared the same way, and then to the earlier
  step and edge. `EdgeSchedule.tied_rounds` counts the rounds in which every
  candidate tied on `metric`, so that the tie break or the order chose.
- **Filling.** It fills every admissible slot: the built-in metrics never fall
  as links are added, so the metric decides which links are active, not how
  many. When no limit binds, every strategy gives the same all-feasible
  schedule; strategies differ only under a binding `max_active`, terminal or
  duty-cycle limit.
- `GreedyEdgeSchedule(LargestEigenvalueMetric())` is the greedy `λ₁` schedule.
- `GreedyEdgeSchedule(LogDetMetric(; epsilon))` is the greedy volume schedule.
  `log det(W + εI)` is monotone and submodular in the activations, the property
  greedy approximation bounds rest on; no bound is claimed here under terminal
  and duty-cycle limits.
- `GreedyEdgeSchedule(SmallestEigenvalueMetric(); tie_break=LogDetMetric(; epsilon))`
  is the greedy `λ_min` schedule. `λ_min` stays zero until every projected
  direction is reachable, so the tie break chooses until then; on a projection
  that keeps directions the links cannot reach, it chooses throughout.
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
- `tied_rounds`: the selection rounds in which every admissible candidate tied
  on the strategy's metric, so that its tie break or the candidate order chose
  (0 for a given constant set).
"""
struct EdgeSchedule
    active::Vector{Vector{EdgePair}}
    gramian::Matrix{Float64}
    terminals::Vector{Dict{EdgePair, NTuple{2, Int}}}
    tied_rounds::Int
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
- `duty_cycle` limits the schedule over the window, counting the steps
  `previous` (active sets, oldest first) just before it. By default it is the
  graphs' `DutyCycleLimit`, which must then be the same in every graph.
"""
function schedule_edges(model::EdgeLinearModel, strategy::AbstractEdgeScheduleStrategy; max_active::Integer,
                        graphs=nothing, feasible=nothing, duty_cycle=nothing, previous=())
    plan = _schedule_plan(model, max_active, graphs, feasible, duty_cycle, previous)
    active, tied = _choose(strategy, plan)
    ordered = [sort(set; by=edge -> findfirst(==(edge), plan.edges)) for set in active]
    terminals = graphs === nothing ? [Dict{EdgePair, NTuple{2, Int}}() for _ in active] :
        [Dict{EdgePair, NTuple{2, Int}}(terminal_assignment(graphs[k], set)) for (k, set) in enumerate(ordered)]
    return EdgeSchedule(ordered, _gramian(plan, ordered), terminals, tied)
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
    duty = duty_cycle
    if duty === nothing && graphs !== nothing
        limits = unique([graph.duty_limit for graph in graphs])
        length(limits) == 1 || throw(ArgumentError(
            "The feasible-edge graphs carry different duty-cycle limits; give the window's limit as duty_cycle."))
        duty = only(limits)
    end
    blocks, _ = _reachability_blocks(model)
    return _SchedulePlan(size(model.projection, 1), copy(edges), blocks, allowed,
        graphs === nothing ? nothing : collect(graphs), Int(max_active), duty, [_pair_set(set) for set in previous])
end

# The duty-cycle history of the schedule `active`: the steps before the window,
# then the window's active sets.
_history(plan::_SchedulePlan, active) = vcat(plan.previous, [Set{EdgePair}(set) for set in active])

# Whether edge `e` can join the edges `active[k]` already active at step `k`,
# with `history` the schedule's duty-cycle history.
function _admits(plan::_SchedulePlan, active::Vector{Vector{EdgePair}}, history, k::Int, e::Int)
    edge = plan.edges[e]
    plan.allowed[e, k] && !(edge in active[k]) && length(active[k]) < plan.max_active || return false
    plan.graphs === nothing || terminal_assignment(plan.graphs[k], vcat(active[k], [edge])) !== nothing || return false
    plan.duty === nothing && return true
    slot = length(plan.previous) + k
    push!(history[slot], edge)
    within = _duty_window_ok(plan.duty, history, slot)
    delete!(history[slot], edge)
    return within
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

# A candidate's metric and tie-break values, each with its tie tolerance.
function _score(metric, tie_break, value, second, n, trace)
    return (; value, tolerance=_tie_tolerance(metric, n, trace, value), second,
        second_tolerance=_tie_tolerance(tie_break, n, trace, second))
end

_score(metric, tie_break, W) =
    _score(metric, tie_break, controllability_metric(metric, W),
        tie_break === nothing ? 0.0 : controllability_metric(tie_break, W), size(W, 1), tr(W))

_same(a, b) = abs(a.value - b.value) <= max(a.tolerance, b.tolerance)

function _beats(score, best)
    best === nothing && return true
    _same(score, best) || return score.value > best.value
    return score.second > best.second + max(score.second_tolerance, best.second_tolerance)
end

# A round's best score and choice, and whether every candidate tied on the metric.
mutable struct _RoundRecord{C}
    best::Any
    choice::C
    lowest::Any
    candidates::Int
end

function _record!(round::_RoundRecord, score, choice)
    round.candidates += 1
    _beats(score, round.best) && (round.best = score; round.choice = choice)
    (round.lowest === nothing || score.value < round.lowest.value) && (round.lowest = score)
    return round
end

_all_tied(round::_RoundRecord) = round.candidates >= 2 && _same(round.best, round.lowest)

function _choose(strategy::GreedyEdgeSchedule, plan::_SchedulePlan)
    N = length(plan.blocks)
    active = [EdgePair[] for _ in 1:N]
    history = _history(plan, active)
    W = zeros(plan.rows, plan.rows)
    tied = 0
    while true
        primary, secondary = _round(strategy.metric, W), _round(strategy.tie_break, W)
        trace = tr(W)
        round = _RoundRecord(nothing, (0, 0), nothing, 0)
        for k in 1:N, e in eachindex(plan.edges)
            _admits(plan, active, history, k, e) || continue
            G = plan.blocks[k][e]
            _record!(round, _score(strategy.metric, strategy.tie_break, _evaluate(primary, G), _evaluate(secondary, G),
                plan.rows, trace + sum(abs2, G)), (k, e))
        end
        round.best === nothing && return active, tied
        _all_tied(round) && (tied += 1)
        k, e = round.choice
        push!(active[k], plan.edges[e])
        push!(history[length(plan.previous) + k], plan.edges[e])
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
        return _constant_active(plan, strategy.edges), 0
    end
    set = EdgePair[]
    tied = 0
    while length(set) < plan.max_active
        round = _RoundRecord(nothing, EdgePair((0, 0)), nothing, 0)
        for edge in plan.edges
            edge in set && continue
            W = _gramian(plan, _constant_active(plan, vcat(set, [edge])))
            _record!(round, _score(strategy.metric, strategy.tie_break, W), edge)
        end
        round.best === nothing && break
        _all_tied(round) && (tied += 1)
        push!(set, round.choice)
    end
    return _constant_active(plan, set), tied
end

# The active sets of a constant edge set: at each step in turn, its admissible
# edges in the set's order.
function _constant_active(plan::_SchedulePlan, set::Vector{EdgePair})
    active = [EdgePair[] for _ in plan.blocks]
    history = _history(plan, active)
    index = Dict(edge => e for (e, edge) in enumerate(plan.edges))
    for k in eachindex(active), edge in set
        if _admits(plan, active, history, k, index[edge])
            push!(active[k], edge)
            push!(history[length(plan.previous) + k], edge)
        end
    end
    return active
end

# 5. Input energy

"""
    minimum_input_energy(model::EdgeLinearModel, schedule, target; initial_state=nothing, rtol=1e-10, reach_tol=1e-8)
        -> (; energy, inputs, residual, reachable, most_negative_input, nonnegative)

The least input energy `Σ |δu|²` that moves the window's projected end state
to `target` with the links of `schedule` (an `EdgeSchedule`, or active sets as
for `reachability_matrix`), from `δx_{k₀} = 0` or from `initial_state`, whose
free response the inputs then correct.

- **Energy and inputs.** `energy` is `yᵀ W⁺ y`, the minimum-norm solution of
  `R u = y` for the required change `y`, with singular values below `rtol`
  times the largest dropped. `inputs[k]` maps each edge active at step `k` to
  its input vector.
- **Reach.** `residual` is `|R u − y| / |y|` (zero for `y = 0`), and
  `reachable` whether it is at most `reach_tol`: a target outside the reachable
  subspace gets the least-squares inputs, `reachable = false`.
- **Signs.** The inputs are signed, in the sign of the model's columns: for
  `edge_input_matrix`, a positive input pushes the ends apart and a negative
  one pulls them together. `most_negative_input` is the smallest input (0 when
  none is negative), and `nonnegative` whether every input is at least
  `−√eps` times the largest magnitude. About zero link force, a laser that
  only pushes realizes the solution only when `nonnegative` holds; otherwise
  the energy is that of a two-way actuator, a lower bound on the push-only one.
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
    lowest = min(minimum(u; init=0.0), 0.0)
    return (; energy=dot(u, u), inputs, residual, reachable=residual <= reach_tol, most_negative_input=lowest,
        nonnegative=lowest >= -sqrt(eps()) * maximum(abs, u; init=0.0))
end

end
