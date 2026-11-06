module ColoredFDJacobianTests

# The implicit solvers' Jacobian (`ColoredFDJacobian`) must equal column-by-column dense
# finite differences with the same per-component step, on the same-satellite entries,
# while perturbing fewer colors than the block size and never perturbing heat_loads.

using Test
using StaticArrays
using ComponentArrays
using SparseArrays
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SE_CFD = SpaceAGORA.SimulationEngine
const EM_CFD = SpaceAGORA.SimulationModel.EnvironmentModels

function cfd_spacecraft(planet, id::Int, n::Int)
    root = Link(root=true, m=500.0, ref_area=12.0)
    links = Link[root]
    for k in 1:2
        push!(links, Link(root=false, m=8.0, ref_area=3.0,
            r=MVector{3, Float64}(1.8 * cos(pi * (k - 1)), 1.8 * sin(pi * (k - 1)), 0.4)))
    end
    ic = InitialCondition(ra=planet.Rp_e + 260e3 + 2e3 * id, rp=planet.Rp_e + 150e3 + 1e3 * id,
        i=53.0, ω=24.0, Ω=10.0, ν=(360.0 * (id - 1)) / n)
    return SpacecraftModel(Joint[], links, root, true, sum(l.m for l in links), 0.0, root.inertia, 0, 0, ic, id)
end

function cfd_args(planet, n_sats::Int)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false,
            normalize=false, save_csv=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, keplerian=true,
            number_of_orbits=1, mission_time=10.0, orientation_sim=false, num_steps_to_save=10, data_rate=1.0),
        environment_model=EnvironmentModel(planet=planet, EI=120.0,
            density_model=EM_CFD.ExponentialAtmosphereModel(planet),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false, ephemerides_model=SimpleEphemeridesModel()),
        dynamics_model=DynamicsModel([cfd_spacecraft(planet, i, n_sats) for i in 1:n_sats],
            (InverseSquaredJ2GravityModel(), AerodynamicCoefficientfM())),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0),
    )
end

# Dense forward FD, one column at a time, FiniteDiff's default per-component step.
function cfd_dense_fd(f!, x0)
    n = length(x0)
    relstep = sqrt(eps(Float64))
    f0 = zeros(n); f!(f0, x0)
    J = zeros(n, n)
    for j in 1:n
        h = max(relstep * abs(x0[j]), relstep)
        x = copy(x0); x[j] += h
        f1 = zeros(n); f!(f1, x)
        J[:, j] = (f1 .- f0) ./ h
    end
    return J
end

@testset "colored per-component FD Jacobian" begin
    planet = Earth()
    n_sats = 3
    args = cfd_args(planet, n_sats)
    u = SE_CFD.build_initial_conditions(args)
    p = ODEParams(n_sats=n_sats, args=args)
    SE_CFD._initialize_save_cache_buffers!(p)
    SE_CFD._initialize_heat_rate_buffers!(p)
    t = 0.0
    ax = ComponentArrays.getaxes(u)
    x0 = Vector{Float64}(u)
    f!(dx, x) = SE_CFD.spacecraft_dynamics!(ComponentArray(dx, ax), ComponentArray(x, ax), p, t)

    proto = SE_CFD._build_block_diagonal_jac_prototype(u, p.is_active)
    jf = SE_CFD.ColoredFDJacobian(u, proto)
    heat = SE_CFD._heat_load_mask(u)
    block = length(u) ÷ n_sats

    @test count(heat) == 3 * 3   # three links per satellite
    @test all(jf.colors[heat] .== 0)
    @test all(jf.colors[.!heat] .> 0)
    @test maximum(jf.colors) == block - 3   # fewer colors than the clique size of the block

    Jd = cfd_dense_fd(f!, x0)
    # The RHS never reads heat_loads: the dense columns are exactly zero. If a term
    # starts to read it, this fails and the exclusion in ColoredFDJacobian is wrong.
    @test all(iszero, Jd[:, heat])

    J = copy(proto)
    jf(J, x0, p, t)
    mask = Matrix(proto) .!= 0
    Jm = Matrix(J)
    scale = maximum(abs, Jd)
    @test maximum(abs, (Jm .- Jd)[mask]) <= 1e-9 * scale
    @test all(iszero, Jm[:, heat])
    # heat_loads ROWS keep their entries (d heat_rate / d state is nonzero).
    @test any(!iszero, Jm[heat, .!heat])
    @test jf.rhs_calls[] == 1 + maximum(jf.colors)

    # single satellite: dense J, one color per non-heat column
    args1 = cfd_args(planet, 1)
    u1 = SE_CFD.build_initial_conditions(args1)
    jf1 = SE_CFD.ColoredFDJacobian(u1, nothing)
    @test maximum(jf1.colors) == length(u1) - 3

    # wiring: implicit full-state modes get jac, explicit ones do not
    for (mode, expect) in ((:rodas5p, true), (:auto_stiff, true), (:tsit5, false))
        prob = SE_CFD._build_typed_solver_problem(u, (0.0, 1.0), p, nothing, mode, proto)
        @test (prob.f.jac isa SE_CFD.ColoredFDJacobian) == expect
    end
    p1 = ODEParams(n_sats=1, args=args1)
    prob1 = SE_CFD._build_typed_solver_problem(u1, (0.0, 1.0), p1, nothing, :rodas5p, nothing)
    @test prob1.f.jac isa SE_CFD.ColoredFDJacobian
end

end # module
