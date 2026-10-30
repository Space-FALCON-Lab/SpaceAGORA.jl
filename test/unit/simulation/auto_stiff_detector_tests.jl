module AutoStiffDetectorTests
# The AutoTsit5 stiffness detector must not read a state component that does not
# change (a constant mass, an idle heat load) as stiffness. OrdinaryDiffEq's own
# Tsit5 estimate gives 0/0 = NaN for such a component, which AutoSwitch treats as
# "stiff"; the solver then left Tsit5 after a handful of steps on a smooth orbit.
# Genuine stiffness must still be detected.
using Test
using SpaceAGORA
const SE = SpaceAGORA.SimulationEngine

# Keplerian orbit with two components whose derivative is identically zero.
function orbit_rhs!(du, u, p, t)
    r = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
    k = -p / r^3
    du[1] = u[4]; du[2] = u[5]; du[3] = u[6]
    du[4] = k * u[1]; du[5] = k * u[2]; du[6] = k * u[3]
    du[7] = 0.0
    du[8] = 0.0
    return nothing
end

# Van der Pol with a large mu (stiff) plus the same kind of constant component.
function vdp_rhs!(du, u, p, t)
    du[1] = u[2]
    du[2] = p * ((1 - u[1]^2) * u[2]) - u[1]
    du[3] = 0.0
    return nothing
end

auto_alg() = SE._auto_tsit5(SE.Rodas5P(autodiff=SE.AutoFiniteDiff()))

@testset "AutoTsit5 stiffness detector" begin
    @testset "smooth orbit with constant components stays on Tsit5" begin
        mu = 3.986004418e14
        r0 = 7.0e6
        u0 = [r0, 0.0, 0.0, 0.0, sqrt(mu / r0), 0.0, 500.0, 0.0]
        prob = SE.ODEProblem(orbit_rhs!, u0, (0.0, 5000.0), mu)
        sol = SE.solve(prob, auto_alg(); reltol=1e-9, abstol=1e-9)
        @test sol.retcode == SE.ReturnCode.Success
        @test all(==(1), sol.alg_choice)
        @test sol.stats.njacs == 0
    end

    @testset "stiff problem still switches to Rodas5P" begin
        prob = SE.ODEProblem(vdp_rhs!, [2.0, 0.0, 1.0], (0.0, 3000.0), 1000.0)
        sol = SE.solve(prob, auto_alg(); reltol=1e-8, abstol=1e-8)
        @test sol.retcode == SE.ReturnCode.Success
        @test 2 in sol.alg_choice
        @test sol.stats.njacs > 0
    end

    @testset "legacy detector is selectable" begin
        @test SE._auto_tsit5(SE.Rodas5P()).choice_function.nonstiffalg isa SE._Tsit5StiffnessProbe
        withenv("SPACEAGORA_AUTO_STIFF_DETECTOR" => "legacy") do
            @test SE._auto_tsit5(SE.Rodas5P()).choice_function.nonstiffalg isa SE.Tsit5
        end
    end
end

end
