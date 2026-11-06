module AutoStiffSwitchTests
using Test, SpaceAGORA, OrdinaryDiffEq
const SE = SpaceAGORA.SimulationEngine

@testset "auto-stiff switch from the Jacobian count" begin
    vdp!(du, u, mu, t) = (du[1] = u[2]; du[2] = mu * (1 - u[1]^2) * u[2] - u[1])
    alg = AutoTsit5(Rodas5P())
    for (mu, tf, expect) in ((1.0e4, 300.0, true), (0.1, 5.0, false))
        prob = ODEProblem(vdp!, [2.0, 0.0], (0.0, tf), mu)
        full = solve(prob, alg; save_everystep=true)
        lean = solve(prob, alg; save_everystep=false, save_on=false)
        @test SE._auto_stiff_switched(full) == expect
        @test SE._auto_stiff_switched_by_jacobians(lean) == expect
        @test lean.u[end] == full.u[end]
    end
end
end
