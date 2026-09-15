using Test
using SpaceAGORA
using StaticArrays
using LinearAlgebra
import SpaceAGORA.TelemetryVerification: make_example_config

const SM_TL = SpaceAGORA.SimulationModel
const CH_TL = SM_TL.ControlHooks

# The Apollo LM's actuators: the descent engine on the body axis plus four RCS
# quads of four jets each, the same layout scripts/dev/viewer_demos/apollo11_landing.jl flies.
const TL_DPS_THRUST_N = 45_040.0
const TL_RCS_THRUST_N = 445.0
const TL_QUAD_ARM_M = 1.65
const TL_QUAD_Z_M = -0.3

function tl_lunar_module_link()
    dps = SM_TL.Thruster(max_thrust=TL_DPS_THRUST_N, location=MVector{3, Float64}(0.0, 0.0, 1.5),
        direction=MVector{3, Float64}(0.0, 0.0, 1.0), Isp=311.0)
    rcs = SM_TL.Thruster[]
    for sx in (-1.0, 1.0), sy in (-1.0, 1.0)
        quad = MVector{3, Float64}(TL_QUAD_ARM_M * sx, TL_QUAD_ARM_M * sy, TL_QUAD_Z_M)
        for d in ((sx, 0.0, 0.0), (0.0, sy, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, -1.0))
            push!(rcs, SM_TL.Thruster(max_thrust=TL_RCS_THRUST_N, location=copy(quad),
                direction=MVector{3, Float64}(d...), Isp=290.0))
        end
    end
    return SM_TL.Link(root=true, m=6_900.0, dims=MVector{3, Float64}(4.2, 4.2, 7.0), ref_area=1.0, thrusters=[dps; rcs])
end

function tl_lunar_module()
    root = tl_lunar_module_link()
    ic = SM_TL.CartesianInitialCondition(SVector(1.8e6, 0.0, 0.0), SVector(0.0, 1.6e3, 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
    return SM_TL.SpacecraftModel(links=[root], root=root, prop_mass=8_200.0,
        inertia_tensor=SMatrix{3, 3, Float64}(2.2e4, 0, 0, 0, 2.4e4, 0, 0, 0, 2.0e4), initial_condition=ic, id=1)
end

# Which quad jets face which body axis, in the order the layout lists them
# (the descent engine is thruster 1, so jet j is thruster j + 1).
tl_jet_axis(j::Int) = ((j - 1) % 4) + 1   # 1 = ±x, 2 = ±y, 3 = +z, 4 = -z

struct TlPlainEffector <: SM_TL.AbstractControlEffectorModel end

@testset "ThrusterLevels" begin
    @testset "hook default is nothing" begin
        @test CH_TL.control_thruster_levels(TlPlainEffector(), 1) === nothing
        @test CH_TL.control_thruster_levels(nothing, 3) === nothing
        @test SpaceAGORA.control_thruster_levels(TlPlainEffector(), 1) === nothing
        # a vehicle without thrusters has no levels to report
        empty_root = SM_TL.Link(root=true, m=10.0, dims=MVector{3, Float64}(1.0, 1.0, 1.0), ref_area=1.0)
        empty_sc = SM_TL.SpacecraftModel(links=[empty_root], root=empty_root,
            initial_condition=SM_TL.CartesianInitialCondition(SVector(1.8e6, 0.0, 0.0), SVector(0.0, 1.6e3, 0.0)), id=1)
        @test CH_TL.descent_thruster_layout(empty_sc) === nothing
    end

    @testset "layout and least-norm RCS allocation" begin
        sc = tl_lunar_module()
        layout = CH_TL.descent_thruster_layout(sc)
        @test layout !== nothing
        @test layout.engine == 1                       # the highest-rated thruster is the descent engine
        @test layout.engine_max_thrust_n == TL_DPS_THRUST_N
        @test length(layout.jets) == 16
        @test layout.jets == collect(2:17)
        A = layout.torque_arms_nm
        @test size(A) == (3, 16)
        # A jet's arm is location x (-max_thrust * exhaust direction): a +x jet on the
        # (+,+) quad pushes -x and yaws the vehicle.
        @test A[:, 1] ≈ cross(SVector(-TL_QUAD_ARM_M, -TL_QUAD_ARM_M, TL_QUAD_Z_M), SVector(TL_RCS_THRUST_N, 0.0, 0.0))

        # The least-norm allocation reproduces a pure torque about each body axis.
        for (axis, τ) in ((:roll, SVector(120.0, 0.0, 0.0)), (:pitch, SVector(0.0, 120.0, 0.0)), (:yaw, SVector(0.0, 0.0, 120.0)))
            u = layout.torque_to_levels * Vector(τ)
            @test A * u ≈ Vector(τ) atol=1e-8
            # least-norm: the solution lies in the row space of A
            @test norm(u - A' * ((A * A') \ Vector(τ))) < 1e-8
            # and it fires only the jets whose arms reach that axis
            idle = axis === :roll ? 1 : axis === :pitch ? 2 : 3      # yaw idles both z-facing sets
            for (j, _) in enumerate(layout.jets)
                a = tl_jet_axis(j)
                reaches = axis === :yaw ? a <= 2 : a != idle
                if reaches
                    @test abs(u[j]) > 1e-9
                else
                    @test abs(u[j]) < 1e-12
                end
            end
        end
    end

    @testset "levels: engine throttle and clipped jets" begin
        sc = tl_lunar_module()
        layout = CH_TL.descent_thruster_layout(sc)
        levels = zeros(17)
        # the descent engine's level is its actual thrust over its rating
        CH_TL.descent_thruster_levels!(levels, layout, 0.4 * TL_DPS_THRUST_N, SVector(0.0, 0.0, 0.0))
        @test levels[1] ≈ 0.4
        @test all(==(0.0), levels[2:end])
        CH_TL.descent_thruster_levels!(levels, layout, 3.0 * TL_DPS_THRUST_N, SVector(0.0, 0.0, 0.0))
        @test levels[1] == 1.0                          # clipped at full thrust
        CH_TL.descent_thruster_levels!(levels, layout, 0.0, SVector(0.0, 0.0, 0.0))
        @test levels[1] == 0.0

        # a commanded yaw torque fires the jets with a positive allocation and
        # leaves the ones the least-norm solution pushes negative idle
        τ = SVector(0.0, 0.0, 400.0)
        raw = layout.torque_to_levels * Vector(τ)
        CH_TL.descent_thruster_levels!(levels, layout, 0.5 * TL_DPS_THRUST_N, τ)
        @test levels[1] ≈ 0.5
        for j in eachindex(layout.jets)
            @test levels[layout.jets[j]] ≈ clamp(raw[j], 0.0, 1.0)
        end
        @test any(>(0.0), levels[2:end])
        @test all(l -> 0.0 <= l <= 1.0, levels)
    end

    @testset "descent control reports levels and the run saves them" begin
        planet = SM_TL.Moon()
        ephem = SM_TL.SimpleEphemeridesModel()
        sc = tl_lunar_module()
        braking, approach = apollo11_descent_targets()
        gcfg = ApolloDescentConfig(site_lat_deg=0.674, site_lon_deg=23.473, approach_azimuth_deg=270.0,
            braking=braking, approach=approach)
        state = ApolloDescentState(1)
        control = ApolloDescentControlModel(ApolloDescentControlConfig(), gcfg, state, NoTerrainModel())
        # before the first control cycle the hook still answers: no thruster is firing
        @test CH_TL.control_thruster_levels(control, 1) == Float64[]
        @test CH_TL.control_thruster_levels(control, 2) === nothing

        base = make_example_config(planet=planet, spacecraft=sc, mission_time=60.0,
            initial_time=SM_TL.InitialTime(year=1969, month=7, day=20, hour=20, minute=10, second=0.0),
            dynamic_effectors=(SM_TL.InverseSquaredGravityModel(),), density_model=SM_TL.NoAtmosphereModel(),
            ephemerides_model=ephem, orientation_sim=true, keplerian=false, EI_km=1.0, verbose=false, results=false)
        args = SM_TL.SimulationConfiguration(file_paths=base.file_paths, simulation_settings=base.simulation_settings,
            mission_configuration=base.mission_configuration, environment_model=base.environment_model,
            dynamics_model=base.dynamics_model, guidance_model=base.guidance_model, navigation_model=base.navigation_model,
            control_model=SM_TL.ControlModel(control_effectors=(control,), control_rates=[0.05]),
            initial_time=base.initial_time, integration_tolerances=base.integration_tolerances)

        counts = SM_TL.SimulationCallbacks.thruster_level_counts(args)
        @test counts == [17]
        fields = SM_TL.SimulationCallbacks.default_save_fields(args)
        @test :thruster_level in [f.name for f in fields]
        field = fields[findfirst(f -> f.name == :thruster_level, fields)]
        @test field.per_satellite && field.column_prefix == "thruster_level"

        # without a reporting effector the field is absent
        plain = SM_TL.SimulationConfiguration(file_paths=args.file_paths, simulation_settings=args.simulation_settings,
            mission_configuration=args.mission_configuration, environment_model=args.environment_model,
            dynamics_model=args.dynamics_model, guidance_model=args.guidance_model, navigation_model=args.navigation_model,
            control_model=SM_TL.ControlModel(control_effectors=(TlPlainEffector(),), control_rates=[1.0]),
            initial_time=args.initial_time, integration_tolerances=args.integration_tolerances)
        @test SM_TL.SimulationCallbacks.thruster_level_counts(plain) == [0]
        @test !(:thruster_level in [f.name for f in SM_TL.SimulationCallbacks.default_save_fields(plain)])
    end
end
