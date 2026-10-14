module VectorEventDirectionTests
using Test
using OrdinaryDiffEq
using SpaceAGORA
const CB = SpaceAGORA.SimulationModel.SimulationCallbacks

@testset "Directional events under the solver's simultaneous-event API" begin
    # Independent linear trajectories cross zero at t=1 (two directions) and
    # t=2 (one direction). A solver step spans both events; root finding must
    # dispatch every simultaneous crossing once, to the appropriate handler.
    function rhs!(du, u, p, t)
        du .= (1.0, -1.0, 1.0)
    end
    function condition!(out, u, t, integrator)
        out .= u
    end
    for directions in (:both, :up, :down)
        observed = Tuple{Float64, Int, Symbol}[]
        up! = directions === :down ? nothing : (i, n) -> push!(observed, (i.t, n, :up))
        down! = directions === :up ? nothing : (i, n) -> push!(observed, (i.t, n, :down))
        callback = CB._directional_vector_callback(condition!, up!, down!, 3)
        prob = ODEProblem(rhs!, [-1.0, 1.0, -2.0], (0.0, 3.0))
        sol = solve(prob, Tsit5(); callback, dt=2.5, abstol=1e-12, reltol=1e-12)
        expected = directions === :both ? [(1.0, 1, :up), (1.0, 2, :down), (2.0, 3, :up)] :
            directions === :up ? [(1.0, 1, :up), (2.0, 3, :up)] : [(1.0, 2, :down)]
        @test string(sol.retcode) == "Success"
        @test length(observed) == length(expected)
        @test [(i, d) for (_, i, d) in observed] == [(i, d) for (_, i, d) in expected]
        @test first.(observed) ≈ first.(expected) atol=1e-12 rtol=0
        @test sol.u[end] ≈ [2.0, -2.0, 1.0] atol=1e-12 rtol=0
    end

    # A direction-filtered event must not invalidate unchanged solver caches.
    for ignored_direction in (-1, 1)
        calls = Ref(0)
        effect!(integrator, idx) = (calls[] += 1)
        up! = ignored_direction < 0 ? effect! : nothing
        down! = ignored_direction > 0 ? effect! : nothing
        cb = CB._directional_vector_callback(condition!, up!, down!, 3)
        prob = ODEProblem(rhs!, [-1.0, 1.0, -2.0], (0.0, 3.0))
        integrator = init(prob, Tsit5(); callback=cb)
        before = copy(integrator.u)
        integrator.derivative_discontinuity = true
        cb.affect!(integrator, Int8[ignored_direction, 0, 0])
        @test calls[] == 0
        @test integrator.u == before
        @test !integrator.derivative_discontinuity
    end
end
end
