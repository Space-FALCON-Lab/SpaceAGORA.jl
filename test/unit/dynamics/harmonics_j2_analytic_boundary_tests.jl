module HarmonicsJ2AnalyticBoundaryTests

using Test, Logging, LinearAlgebra, StaticArrays, SpaceAGORA
import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft
const SM = SpaceAGORA.SimulationModel
const PE = SM.DynamicEffectors.PerturbationEffectors
const ES = SM.EffectorSampling
const ZERO = SVector(0.0,0.0,0.0)
const EYE = SMatrix{3,3,Float64,9}(I)
const TURN = @SMatrix [0.0 0.0 1.0; 0.0 1.0 0.0; -1.0 0.0 0.0]
const GM = 2.0
const RADIUS = 1.0
const J2 = 0.01

function model_at(directory, planet; include_central=true, j2=J2)
    path = joinpath(directory, "j2-$(j2).csv")
    write(path, "degree,order,C,S\n0,0,1.0,0.0\n1,0,0.0,0.0\n2,0,$(-j2/sqrt(5.0)),0.0\n")
    return PE.GravitationalHarmonicsModel(2,0,path,planet; coefficient_normalization=:full,
        reference_radius_m=RADIUS,gm_m3s2=GM,include_central)
end

# Independent potential-energy oracle, not the Pines recurrence or another
# SpaceAGORA acceleration helper. Differentiate Φ with a symmetric BigFloat
# difference at 256 bits; h=1e-20 makes truncation/cancellation negligible beside
# Float64 rounding for these order-one synthetic distances.
function analytic_acceleration(position; central=true)
    return setprecision(BigFloat,256) do
        p = collect(BigFloat.(position))
        h = big"1e-20"
        potential(q) = begin
            r = sqrt(sum(abs2,q))
            (central ? -BigFloat(GM)/r : zero(r)) +
                BigFloat(GM)*BigFloat(J2)*BigFloat(RADIUS)^2/(2r^3)*(3q[3]^2/r^2-1)
        end
        values = ntuple(3) do axis
            plus, minus = copy(p), copy(p)
            plus[axis] += h; minus[axis] -= h
            Float64(-(potential(plus)-potential(minus))/(2h))
        end
        SVector{3,Float64}(values)
    end
end

function sampled(position, planet; rotation=EYE,mass=7.0)
    p = SVector{3,Float64}(position)
    state = ES.StateSample(p,ZERO,mass)
    frame = ES.PlanetFrameSample(rotation,rotation*p,ZERO,norm(p)-planet.Rp_e,0.0,0.0)
    return state,ES.EnvironmentSample(planet;planet_frame=frame)
end

function parameters(planet, model)
    craft = make_three_body_spacecraft(bus_dims=(1.0,1.0,1.0), panel_dims=(0.01,1.0,0.5),
        bus_mass=100.0,panel_mass_each=1.0,panel_offset_y=1.0,
        ic=SM.InitialCondition(ra=7.0e6,rp=7.0e6,i=0.0,ω=0.0,Ω=0.0,ν=0.0),prop_mass=0.0,id=1)
    args = make_example_config(; planet,spacecraft=craft,mission_time=0.1,initial_time=SM.InitialTime(year=2025,month=1,day=1,hour=0,minute=0,second=0.0),
        dynamic_effectors=(model,),density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(),orientation_sim=false,
        keplerian=true,verbose=false,results=false)
    return SM.ODEParams(n_sats=1,args=args)
end

@testset "J2 harmonic wrenches and relative gravity match the independent potential" begin
    planet = make_no_gram_planet(:earth)
    mktempdir() do directory
        for central in (false,true)
            model = model_at(directory,planet;include_central=central)
            @test ES.environment_requirements(model).planet_frame
            @test ES.gravity_backbone_structure(model) == :position_only_static_gravity
            for position in (SVector(3.0,0.0,0.0), SVector(0.0,0.0,3.0), SVector(2.3,-1.2,0.7)), rotation in (EYE,TURN)
                state,environment = sampled(position,planet;rotation)
                expected = rotation' * analytic_acceleration(rotation*position;central)
                force,torque = ES.wrench(model,state,environment,0.0)
                @test force/state.mass_kg ≈ expected rtol=3e-13 atol=2e-16
                @test torque === ZERO
                @test ES.gravity_backbone_acceleration_ii(model,state,environment,0.0) ≈ expected rtol=3e-13 atol=2e-16
            end
            base = SVector(3.0,0.4,-0.2)
            displacement = SVector(0.2,-0.03,0.01)
            far = base+displacement
            xb,eb = sampled(base,planet;mass=11.0)
            xf,ef = sampled(far,planet;mass=29.0)
            expected = analytic_acceleration(far;central)-analytic_acceleration(base;central)
            relative = ES.gravity_backbone_relative_acceleration_ii(model,xb,displacement,xf,eb,ef,0.0)
            @test relative ≈ expected rtol=3e-12 atol=2e-16
            @test ES.gravity_backbone_relative_acceleration_ii(model,xb,ZERO,xb,eb,eb,0.0) == ZERO
            unavailable = ES.EnvironmentSample(planet)
            @test_throws ArgumentError ES.wrench(model,xb,unavailable,0.0)
            @test_throws ArgumentError ES.gravity_backbone_relative_acceleration_ii(model,xb,displacement,xf,unavailable,ef,0.0)
            @test_throws ArgumentError ES.gravity_backbone_relative_acceleration_ii(model,xb,displacement,xf,eb,unavailable,0.0)
        end
    end
end

@testset "Optional J2 diagnostics preserve force and bound their logging" begin
    planet = make_no_gram_planet(:earth)
    mktempdir() do directory
        model = model_at(directory,planet)
        parameter = parameters(planet,model)
        position = SVector(2.3,-1.2,0.7)
        state = [position...,0.0,0.0,0.0,7.0]
        old_flag,old_count = PE._DEBUG_COMPARE_J2[],PE._J2_COMPARE_LOG_COUNT[]
        try
            PE._DEBUG_COMPARE_J2[] = false
            baseline = PE._harmonics_calcforcetorque_with_lpi(model,state,parameter,Int64(1),TURN)
            PE._DEBUG_COMPARE_J2[] = true
            PE._J2_COMPARE_LOG_COUNT[] = 0
            debugged = @test_logs (:info,"J2 generic-vs-analytic comparison") PE._harmonics_calcforcetorque_with_lpi(model,state,parameter,Int64(1),TURN)
            @test debugged == baseline
            @test debugged[1]/7.0 ≈ TURN' * analytic_acceleration(TURN*position) rtol=3e-13 atol=2e-16
            @test debugged[2] === ZERO
            @test PE._J2_COMPARE_LOG_COUNT[] == 1
            PE._J2_COMPARE_LOG_COUNT[] = 20
            silent = @test_logs min_level=Logging.Info PE._harmonics_calcforcetorque_with_lpi(model,state,parameter,Int64(1),TURN)
            @test silent == baseline
            @test PE._J2_COMPARE_LOG_COUNT[] == 21

            spherical = model_at(directory,planet;j2=0.0)
            count_before = PE._J2_COMPARE_LOG_COUNT[]
            central = @test_logs min_level=Logging.Info PE._harmonics_calcforcetorque_with_lpi(spherical,state,parameter,Int64(1),TURN)
            @test central[1]/7.0 ≈ -GM*position/norm(position)^3
            @test PE._J2_COMPARE_LOG_COUNT[] == count_before
        finally
            PE._DEBUG_COMPARE_J2[] = old_flag
            PE._J2_COMPARE_LOG_COUNT[] = old_count
        end
        @test PE._DEBUG_COMPARE_J2[] == old_flag
        @test PE._J2_COMPARE_LOG_COUNT[] == old_count
    end
end

end # module HarmonicsJ2AnalyticBoundaryTests
