module ControlReadoutTests
using Test, StaticArrays, SpaceAGORA
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const CH = SM.ControlHooks
const ZERO3 = SVector(0.0, 0.0, 0.0)

# These tests call the public-to-the-engine entry through real ODEParams and
# SimulationConfiguration, including its abstractly stored ControlModel field.
# No propagation, SPICE data, atmosphere data or callbacks are needed.
function readout_params(effectors::Tuple; n_sats::Int=32)
    planet = SM.Earth()
    spacecraft = SM.SpacecraftModel[]
    for i in 1:n_sats
        root = SM.Link(root=true, m=500.0, ref_area=12.0)
        ic = SM.InitialCondition(ra=planet.Rp_e + 550e3, rp=planet.Rp_e + 550e3,
            i=53.0, ω=0.0, Ω=10.0, ν=360.0 * (i - 1) / n_sats)
        push!(spacecraft, SM.SpacecraftModel(SM.Joint[], [root], root, true,
            500.0, 0.0, root.inertia, 0, 0, ic, i))
    end
    args = SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime,
            keplerian=true, number_of_orbits=1, mission_time=1.0,
            orientation_sim=true, num_steps_to_save=1),
        environment_model=SM.EnvironmentModel(planet=planet, EI=120.0,
            density_model=SM.NoAtmosphereModel(),
            ephemerides_model=SM.SimpleEphemeridesModel(),
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false),
        dynamics_model=SM.DynamicsModel(spacecraft, (SM.InverseSquaredGravityModel(),)),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=effectors,
            control_rates=fill(1.0, length(effectors))),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
    return SM.ODEParams(n_sats=n_sats, args=args)
end

struct ReadoutProbe{ID,F,T,W} <: SM.AbstractControlEffectorModel
    force::F
    torque::T
    wheel_torque::W
    mass_rate::Float64
    sat_idx::Int
    calls::Vector{Tuple{Symbol,Symbol,Int}}
end
function probe(id::Symbol, force, torque, wheel_torque, mass_rate, calls; sat_idx=999)
    return ReadoutProbe{id,typeof(force),typeof(torque),typeof(wheel_torque)}(
        force, torque, wheel_torque, mass_rate, sat_idx, calls)
end
function CH.calcControlForceTorque(m::ReadoutProbe{ID}, u::AbstractVector,
        p::SM.ODEParams, i::Int64, t::Float64) where {ID}
    push!(m.calls, (ID, :force, i))
    return m.force, m.torque
end
function CH.calcControlMassFlowRate(m::ReadoutProbe{ID}, u::AbstractVector,
        p::SM.ODEParams, i::Int64, t::Float64) where {ID}
    push!(m.calls, (ID, :mass, i))
    return m.mass_rate
end
function CH.calcReactionWheelTorque(m::ReadoutProbe{ID}, u::AbstractVector,
        p::SM.ODEParams, i::Int64, t::Float64) where {ID}
    push!(m.calls, (ID, :wheel, i))
    return m.wheel_torque
end

# ControlModel accepts arbitrary extension types. Only the force hook is
# implemented here, so generic mass-flow and reaction-wheel fallbacks matter.
struct GlobalForceOnly
    calls::Vector{Tuple{Symbol,Symbol,Int}}
end
function CH.calcControlForceTorque(m::GlobalForceOnly, u::AbstractVector, p::SM.ODEParams,
        i::Int64, t::Float64)
    push!(m.calls, (:global, :force, i))
    return SVector(0.0, 0.0, Float64(i)), ZERO3
end

function readout(p, i; force=ZERO3, torque=ZERO3, wheel=ZERO3)
    f, tq, rw = MVector{3,Float64}(force), MVector{3,Float64}(torque), MVector{3,Float64}(wheel)
    rate = SE._accumulate_control_effectors!(f, tq, rw, Float64[], p, i, 2.0, false)
    return SVector{3,Float64}(f), SVector{3,Float64}(tq), SVector{3,Float64}(rw), rate
end

function held_managers(n)
    return ntuple(n) do i
        SM.MagneticMomentumManagerModel(sat_idx=i,
            held_torque_body=SVector(i / 1024.0, i / 512.0, -i / 256.0))
    end
end

@testset "control readout retains held torque ownership without updating control state" begin
    managers = held_managers(32)
    p = readout_params(managers)
    f0, tq0, rw0 = SVector(2.0,-3.0,5.0), SVector(-2.0,3.0,-5.0), SVector(7.0,11.0,13.0)
    for i in 1:32
        f, tq, rw, rate = readout(p, i; force=f0, torque=tq0, wheel=rw0)
        @test f == f0
        @test tq == tq0 + managers[i].held_torque_body
        @test rw == rw0
        @test rate === 0.0
    end
    @test all(m -> m.ticks == 0 && !m.initialized && isnan(m.last_update_s), managers)
    @test all(m -> m.h_wheels == ZERO3, managers)
end

@testset "heterogeneous and global control readout preserves hook order and vector shapes" begin
    calls = Tuple{Symbol,Symbol,Int}[]
    # The sat_idx field is deliberately unrelated to the evaluated spacecraft:
    # field naming is not a contract permitting generic ownership filtering.
    vector_model = probe(:vector, [1.0,2.0,3.0], SVector(4.0,5.0,6.0), nothing, -2.0, calls)
    static_model = probe(:static, SVector(-4.0,5.0,-6.0), [-2.0,4.0,5.0],
        [3.0,-4.0,5.0], -3.0, calls; sat_idx=1)
    p = readout_params((vector_model, GlobalForceOnly(calls), static_model); n_sats=2)
    for i in 1:2
        empty!(calls)
        f, tq, rw, rate = readout(p, i; force=SVector(10.0,20.0,30.0),
            torque=SVector(1.0,2.0,3.0), wheel=SVector(7.0,8.0,9.0))
        @test f == SVector(7.0,27.0,27.0 + i)
        @test tq == SVector(3.0,11.0,14.0)
        @test rw == SVector(10.0,4.0,14.0)
        @test rate === -5.0
        @test calls == [(:vector,:force,i), (:vector,:mass,i), (:vector,:wheel,i),
            (:global,:force,i), (:static,:force,i), (:static,:mass,i), (:static,:wheel,i)]
    end
end

@testset "control readout sums left to right and filters nonfinite mass rates" begin
    calls = Tuple{Symbol,Symbol,Int}[]
    ids = (:large, :cancel, :unit, :nan, :positive_inf, :negative_inf)
    rates = (1.0e16, -1.0e16, 1.0, NaN, Inf, -Inf)
    effectors = ntuple(length(ids)) do i
        # Force and wheel sums also expose accidental regrouping. Nonfinite
        # mass rates must not suppress the same controller's finite wrench.
        x = i <= 3 ? rates[i] : 2.0
        probe(ids[i], SVector(x,0.0,0.0), SVector(0.0,x,0.0),
            SVector(0.0,0.0,x), rates[i], calls)
    end
    p = readout_params(effectors; n_sats=1)
    f, tq, rw, rate = readout(p, 1)
    @test rate === 1.0
    @test f == SVector(7.0,0.0,0.0)
    @test tq == SVector(0.0,7.0,0.0)
    @test rw == SVector(0.0,0.0,7.0)
    @test calls == [(id,hook,1) for id in ids for hook in (:force,:mass,:wheel)]
end

@testset "empty control tuple preserves existing accumulators" begin
    p = readout_params((); n_sats=1)
    initial = (SVector(1.0,2.0,3.0), SVector(4.0,5.0,6.0), SVector(7.0,8.0,9.0))
    f, tq, rw, rate = readout(p, 1; force=initial[1], torque=initial[2], wheel=initial[3])
    @test (f,tq,rw) == initial
    @test rate === 0.0
end

# Warmed, allocation-only regression through the real type-erasing config
# boundary. Keep buffers/state outside the measurement and observe a checksum.
# The tuple is still obtained through an abstract ControlModel field, so this
# entry retains one runtime dispatch boundary per spacecraft. ODEParams is a
# large immutable aggregate; boxing it at that boundary can allocate even when
# the homogeneous inner loop infers. Julia 1.12.1 measured 2,096 B/call for this
# fixture after warmup. A 4 KiB/call budget allows that measured residual while
# catching a return to the original per-controller hook-result boxing. This is
# an allocation regression bound, not a zero-allocation or wall-time claim.
@noinline function readout_batch!(p, state, forces, torques, wheels, repetitions)
    checksum = 0.0
    for k in 1:repetitions
        fill!(forces, 0.0); fill!(torques, 0.0); fill!(wheels, 0.0)
        i = mod1(k, p.n_sats)
        rate = SE._accumulate_control_effectors!(forces, torques, wheels, state, p, i, 2.0, false)
        checksum += forces[1] + torques[1] + wheels[1] + rate
    end
    return checksum
end
function readout_allocated(p, state, forces, torques, wheels, repetitions, checksum)
    return @allocated checksum[] = readout_batch!(p, state, forces, torques, wheels, repetitions)
end

@testset "warmed 32-manager control readout has bounded allocation" begin
    p = readout_params(held_managers(32))
    state = Float64[]
    f, tq, rw = MVector{3,Float64}(ZERO3), MVector{3,Float64}(ZERO3), MVector{3,Float64}(ZERO3)
    checksum = Ref(0.0)
    repetitions = 32
    readout_batch!(p, state, f, tq, rw, repetitions)
    readout_allocated(p, state, f, tq, rw, repetitions, checksum)
    bytes = readout_allocated(p, state, f, tq, rw, repetitions, checksum)
    @test checksum[] == sum(1:32) / 1024.0
    @test bytes <= 4096 * repetitions
end
end # module ControlReadoutTests
