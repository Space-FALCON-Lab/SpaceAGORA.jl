using Test
using SpaceAGORA
using SparseArrays
using ComponentArrays

const SE = SpaceAGORA.SimulationEngine
const SM = SpaceAGORA.SimulationModel

# Stand-ins for the shapes the gate has to tell apart. The real coupling
# effectors (RPOMPCControlModel, the RPO guidance models) carry chaser_idx and
# target_idx and read both satellites' states inside the RHS; a per-satellite
# effector carries at most one satellite index and reads only that satellite.
struct _CouplingEffectorStub
    chaser_idx::Int
    target_idx::Int
end

struct _SingleSatEffectorStub
    chaser_idx::Int
end

struct _PlainEffectorStub end

@testset "cross-satellite coupling detection" begin
    @test SE._effector_couples_satellites(_CouplingEffectorStub(1, 2))
    @test SE._effector_couples_satellites(_CouplingEffectorStub(2, 1))

    # Same satellite named twice is not coupling: the derivative still depends
    # on that one block only.
    @test !SE._effector_couples_satellites(_CouplingEffectorStub(3, 3))

    # One index alone names the satellite the effector acts on, not a second one.
    @test !SE._effector_couples_satellites(_SingleSatEffectorStub(1))
    @test !SE._effector_couples_satellites(_PlainEffectorStub())
    @test !SE._effector_couples_satellites(SM.InverseSquaredJ2GravityModel())

    # The real coupling effectors must be recognised, so the gate cannot rot
    # away from the types it exists to catch. Both default to chaser 1 /
    # target 2 and read both satellites' states inside the RHS.
    @test SE._effector_couples_satellites(SM.RPOMPCControlModel())
    @test SE._effector_couples_satellites(SM.RPOGuidanceModel())
end

@testset "coupling detection scans every effector surface" begin
    plain = (dynamic_effectors=(_PlainEffectorStub(),),)
    empty_surface = (guidance_effectors=(),)

    uncoupled = (
        dynamics_model=(dynamic_effectors=(SM.InverseSquaredJ2GravityModel(),),),
        guidance_model=(guidance_effectors=(),),
        navigation_model=(navigation_effectors=(),),
        control_model=(control_effectors=(_SingleSatEffectorStub(1),),),
    )
    @test !SE._config_has_cross_satellite_coupling(uncoupled)

    # Coupling has to be found on each of the four surfaces, not just controls.
    for surface in (
        (:dynamics_model, :dynamic_effectors),
        (:guidance_model, :guidance_effectors),
        (:navigation_model, :navigation_effectors),
        (:control_model, :control_effectors),
    )
        model_name, effectors_name = surface
        args = merge(
            uncoupled,
            NamedTuple{(model_name,)}((
                NamedTuple{(effectors_name,)}(((_CouplingEffectorStub(1, 2),),)),
            )),
        )
        @test SE._config_has_cross_satellite_coupling(args)
    end

    # A config missing a surface entirely must not error.
    @test !SE._config_has_cross_satellite_coupling((dynamics_model=plain,))
    @test !SE._config_has_cross_satellite_coupling((guidance_model=empty_surface,))
    @test !SE._config_has_cross_satellite_coupling(NamedTuple())
end

@testset "block-diagonal prototype structure" begin
    u0 = ComponentVector(sc=[(pos=zeros(3), vel=zeros(3)) for _ in 1:4])
    J = SE._build_block_diagonal_jac_prototype(u0)
    n = length(u0)

    @test size(J) == (n, n)
    @test nnz(J) == 4 * 36          # four 6x6 dense blocks

    # Every stored entry must lie inside one satellite's own block, and every
    # in-block entry must be present. This is the assertion the gate protects:
    # an effector reading another satellite needs an entry that is structurally
    # absent here.
    rows_, cols_, _ = findnz(J)
    @test all(((r, c),) -> (r - 1) ÷ 6 == (c - 1) ÷ 6, zip(rows_, cols_))
    for blk in 0:3, j in 1:6, k in 1:6
        @test J[blk * 6 + k, blk * 6 + j] != 0.0
    end
end

@testset "split problems carry the prototype on the inner component function" begin
    # A SplitFunction-level prototype is silently ignored by the coloring, so
    # the prototype must land on the ODEFunction wrapping the component.
    J = sparse([1, 2], [1, 2], [1.0, 1.0], 2, 2)
    f(du, u, p, t) = (du .= u; nothing)

    wrapped = SE._split_component_function(f, J)
    @test wrapped !== f
    @test hasproperty(wrapped, :jac_prototype)
    @test wrapped.jac_prototype === J
    @test SE._has_sparse_jac_prototype(wrapped)

    # No prototype means the bare function, so the explicit paths keep their
    # existing dispatch and allocate nothing extra.
    @test SE._split_component_function(f, nothing) === f
    @test !SE._has_sparse_jac_prototype(f)
end

@testset "inactive satellites shrink to their diagonal" begin
    u0 = ComponentVector(sc=[(pos=zeros(3), vel=zeros(3)) for _ in 1:4])
    full = SE._build_block_diagonal_jac_prototype(u0)
    @test nnz(full) == 4 * 36

    # All active must match the no-argument form exactly.
    @test SE._build_block_diagonal_jac_prototype(u0, Bool[1, 1, 1, 1]) == full

    # Two inactive: their blocks collapse to 6 diagonal entries each.
    part = SE._build_block_diagonal_jac_prototype(u0, Bool[1, 0, 1, 0])
    @test nnz(part) == 2 * 36 + 2 * 6

    # The retained structure must be a subset of the full pattern, never new
    # entries, and every active block must survive intact.
    rows_p, cols_p, _ = findnz(part)
    @test all(((r, c),) -> full[r, c] != 0.0, zip(rows_p, cols_p))
    for blk in (0, 2), j in 1:6, k in 1:6
        @test part[blk * 6 + k, blk * 6 + j] != 0.0
    end
    # An inactive block keeps only its diagonal.
    for blk in (1, 3), j in 1:6, k in 1:6
        expected_nonzero = (j == k)
        @test (part[blk * 6 + k, blk * 6 + j] != 0.0) == expected_nonzero
    end
end

@testset "direct-CSC path agrees with the general path" begin
    # The contiguous fast path and the index-list fallback must build the same
    # matrix; otherwise the optimisation changes the Jacobian it describes.
    for n_sats in (1, 3, 8)
        u0 = ComponentVector(
            sc=[(pos=zeros(3), vel=zeros(3), mass=0.0, heat_loads=zeros(2)) for _ in 1:n_sats])
        fast = SE._build_block_diagonal_jac_prototype(u0)
        @test fast isa SparseMatrixCSC{Float64, Int}
        @test size(fast) == (length(u0), length(u0))
        # Rebuild the same pattern from first principles.
        b = 9
        expected = spzeros(length(u0), length(u0))
        for s in 0:(n_sats - 1), j in 1:b, k in 1:b
            expected[s * b + k, s * b + j] = 1.0
        end
        @test fast == expected
    end
end
