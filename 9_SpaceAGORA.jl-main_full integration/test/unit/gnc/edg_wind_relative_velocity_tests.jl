module EDGWindRelativeVelocityTests
# The E-EDG targeting and control environment and the aerodynamic wrench must
# agree on the atmosphere-relative velocity.
#
# Convention: a density model's wind is the air's own velocity relative to the
# rotating planet, in local east/north/up components (GRAM's ewWind, nsWind and
# verticalWind; GRAM prints ewWind as "Eastward Wind"). `vel_pp` is the
# spacecraft's velocity relative to the rotating planet in planet-fixed axes.
# The velocity the flow sees is therefore v_rel = vel_pp - w, with w the wind
# rotated into planet-fixed axes.
using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const EM = SM.EnvironmentModels
const CB = SM.SimulationCallbacks
const CH = SM.ControlHooks

# A constant density, temperature and wind everywhere.
struct ConstantWindModel <: SM.AbstractDensityModel
    wind_enu::SVector{3, Float64}
end
EM.getDensity(m::ConstantWindModel, h::Float64, lat::Float64, lon::Float64,
    t::Float64, wind::Bool) = (1.0e-9, 700.0, m.wind_enu)
EM.getDensity(m::ConstantWindModel, h::Float64, lat::Float64, lon::Float64,
    t::Float64, wind::Bool, p) = EM.getDensity(m, h, lat, lon, t, wind)

function config(density_model; planet=SM.Earth())
    root = Link(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(ra=planet.Rp_e + 400.0e3, rp=planet.Rp_e + 150.0e3,
        i=53.0, ω=20.0, Ω=40.0, ν=-10.0)
    sc = SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, 1)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=MissionConfiguration(MissionTime, true, 1, 20.0, false, 20, 2.0),
        environment_model=EnvironmentModel(planet=planet, EI=600.0,
            density_model=density_model,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=true, ephemerides_model=SimpleEphemeridesModel()),
        dynamics_model=DynamicsModel([sc], (InverseSquaredGravityModel(), AerodynamicCoefficientfM())),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
    )
end

# Local east/north/up unit vectors in planet-fixed axes, written out here
# rather than taken from `latlongtoNED`, so that the test pins the convention.
function enu_axes(lat, lon)
    e = SVector(-sin(lon), cos(lon), 0.0)
    n = SVector(-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat))
    u = SVector(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
    return e, n, u
end

function setup(model)
    args = config(model)
    p = SM.ODEParams(n_sats=1, args=args)
    p.shared_buffers.et_start[] = SM.ephemerides_time_seconds(args.initial_time,
        args.environment_model.ephemerides_model)
    p.shared_buffers.callback_env_config[] = CB._snapshot_callback_env_config(args)
    p.shared_buffers.current_time[] = 0.0
    u = SE.build_initial_conditions(args)
    pos = SVector{3, Float64}(u.sc[1].pos)
    vel = SVector{3, Float64}(u.sc[1].vel)
    planet = args.environment_model.planet
    pos_pp, vel_pp = SM.FrameTransforms.r_intor_p!(pos, vel, planet, p.shared_buffers.et_start[],
        args.environment_model.ephemerides_model)
    lla = SM.FrameTransforms.rtolatlong(pos_pp, planet)
    return (; args, p, u, pos, vel, pos_pp, vel_pp, lat=lla[2], lon=lla[3])
end

# The wind rotated into planet-fixed axes.
function wind_pp(s, wind_enu)
    e, n, u = enu_axes(s.lat, s.lon)
    return wind_enu[1] * e + wind_enu[2] * n + wind_enu[3] * u
end

@testset "E-EDG relative velocity is v - w, as in aerodynamics" begin
    wind_enu = SVector(60.0, -35.0, 4.0)
    s = setup(ConstantWindModel(wind_enu))
    expected = s.vel_pp - wind_pp(s, wind_enu)

    withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "off", "SPACEAGORA_GRAM_PROCESS_POOL" => "off") do
        env = CH._edg_targeting_prediction_environment(s.p, s.pos, s.vel, 0.0)
        @test env.vel_pp ≈ s.vel_pp rtol = 1e-12
        @test env.vel_pp_rw ≈ expected rtol = 1e-12
        @test env.speed ≈ norm(expected) rtol = 1e-12

        control_env = CH._edg_environment_state(s.u, s.p, 0.0, 1)
        @test control_env.speed ≈ norm(expected) rtol = 1e-12
        @test control_env.dynamic_pressure ≈ 0.5 * 1.0e-9 * norm(expected)^2 rtol = 1e-12

        # The aerodynamic wrench's drag opposes the same relative velocity.
        cb = CB.get_density_callback(1, s.args.dynamics_model.dynamic_effectors, s.args)
        cb.affect!((p=s.p, u=s.u, t=0.0))
        SM.calcForceTorque(AerodynamicCoefficientfM(), s.u.sc[1], s.p, 1)
        drag_ii = s.p.save_cache.drag_cache[1]
        @test norm(drag_ii) > 0.0
        drag_hat_pp = env.l_pi * (drag_ii / norm(drag_ii))
        @test drag_hat_pp ≈ -env.vel_pp_rw / env.speed atol = 1e-10
    end
end

@testset "a tailwind lowers the E-EDG airspeed" begin
    # Air moving with the spacecraft at 40% of its planet-relative velocity.
    s0 = setup(ConstantWindModel(SVector(0.0, 0.0, 0.0)))
    e, n, u = enu_axes(s0.lat, s0.lon)
    w = 0.4 * s0.vel_pp
    s = setup(ConstantWindModel(SVector(dot(w, e), dot(w, n), dot(w, u))))
    env = CH._edg_targeting_prediction_environment(s.p, s.pos, s.vel, 0.0)
    @test env.speed ≈ 0.6 * norm(s.vel_pp) rtol = 1e-9
    @test CH._edg_environment_state(s.u, s.p, 0.0, 1).speed ≈ 0.6 * norm(s.vel_pp) rtol = 1e-9
end

end # module
