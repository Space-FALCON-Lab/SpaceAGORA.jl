module LinearToleranceTests
using Test
using SpaceAGORA
using OrdinaryDiffEq
using ComponentArrays
using LinearSolve
using SparseArrays
using LinearAlgebra
const SE = SpaceAGORA.SimulationEngine

@testset "Implicit integration preserves componentwise error tolerances" begin
    u0 = ComponentVector(x=1.0, y=2.0)
    rtol = ComponentVector(x=1e-9, y=1e-10)
    atol = ComponentVector(x=1e-10, y=1e-11)
    rhs!(du, u, p, t) = (du .= -u; nothing)
    for sparse_jac in (false, true)
        f = sparse_jac ? ODEFunction(rhs!; jac_prototype=spdiagm(ones(2))) : rhs!
        prob = ODEProblem(f, u0, (0.0, 1.0))
        alg = Rodas5P(autodiff=SE.AutoFiniteDiff(),
            linsolve=SE._sparse_linsolve_or_default(sparse_jac))
        integ = init(prob, alg; reltol=rtol, abstol=atol)
        @test integ.opts.reltol == rtol
        @test integ.opts.abstol == atol
        sol = solve!(integ)
        @test string(sol.retcode) == "Success"
        @test sol.u[end] ≈ u0 .* exp(-1) rtol=1e-8
        @test integ.opts.reltol == rtol
        @test integ.opts.abstol == atol
    end
end

@testset "Linear tolerance adapter retains backend and cache reuse" begin
    A = [3.0 1.0; 1.0 4.0]
    b = [2.0, 3.0]
    for sparse_jac in (false, true), componentwise in (false, true)
        matrix = sparse_jac ? sparse(A) : A
        tol = componentwise ? [1e-9, 1e-10] : 1e-10
        alg = SE._sparse_linsolve_or_default(sparse_jac)
        prob = LinearProblem(matrix, b)
        cache = init(prob, alg; reltol=tol, abstol=tol)
        direct = init(prob, alg.algorithm; reltol=1e-10, abstol=1e-10)
        @test typeof(cache.alg) === typeof(direct.alg)
        @test cache.reltol == cache.abstol == 1e-10
        @test solve!(cache).u ≈ A \ b rtol=1e-12
        next_A = 2 .* matrix
        next_b = [4.0, -1.0]
        # The normal in-place factorization may overwrite the reinit matrix.
        expected = next_A \ next_b
        LinearSolve.reinit!(cache; A=next_A, b=next_b)
        @test solve!(cache).u ≈ expected rtol=1e-12
    end
end
end
