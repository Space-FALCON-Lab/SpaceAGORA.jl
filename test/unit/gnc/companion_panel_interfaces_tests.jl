module CompanionPanelInterfacesTests
using Test, SpaceAGORA, StaticArrays, LinearAlgebra
const SM = SpaceAGORA.SimulationModel

@testset "Analytic density slope and clamp conventions" begin
    rho(m,h)=getDensity(m,Float64(h),0.,0.,0.,false,(args=(environment_model=(planet=make_no_gram_planet(:earth),),),))[1]
    model=PolynomialFitAtmosphereModel([-0.001,-0.1,-2.];
        valid_min_altitude_m=100e3,valid_max_altitude_m=200e3)
    for h in (120e3,150e3,180e3)
        numerical=(rho(model,h+0.1)-rho(model,h-0.1))/0.2
        @test density_altitude_derivative(model,h) ≈ numerical rtol=2e-8
    end
    for h in (0.,99_999.,200_001.,300e3)
        @test density_altitude_derivative(model,h) == 0.
    end
    @test density_altitude_derivative(model,100e3) ≈ rho(model,100e3)*(-0.3)/1e3
    collapsed=PolynomialFitAtmosphereModel([-0.1,0.];
        valid_min_altitude_m=100e3,valid_max_altitude_m=100e3)
    @test rho(collapsed,0.) == rho(collapsed,100e3) == rho(collapsed,200e3)
    @test density_altitude_derivative(collapsed,100e3) == 0.
    saturated=PolynomialFitAtmosphereModel([10.,1e6])
    @test density_altitude_derivative(saturated,150e3) == 0.
    low=PolynomialFitAtmosphereModel([-10.,-1e6])
    @test density_altitude_derivative(low,150e3) == 0.
    piecewise=PiecewiseExponentialAtmosphereModel([0.,100.,200.],[1.,0.5],[20.,50.])
    @test density_altitude_derivative(piecewise,100.) == -0.5/50.
    @test density_altitude_derivative(piecewise,prevfloat(100.)) ≈ -exp(-5)/20.
    exponential=ExponentialAtmosphereModel(1.,0.,20.;valid_min_altitude_m=0.,valid_max_altitude_m=100.)
    @test density_altitude_derivative(exponential,200.) ≈ -exp(-10)/20.
    @test density_altitude_derivative(NoAtmosphereModel(),100.) == 0.
end

function vehicle()
    make_three_body_spacecraft(bus_dims=(1.,2.,3.),panel_dims=(0.01,2.,1.),
        bus_mass=50.,panel_mass_each=5.,panel_offset_y=2.,
        ic=InitialCondition(ra=7e6,rp=7e6,i=0.,Ω=0.,ω=0.,ν=0.))
end
@testset "Public panel actuation and commanded-area drag" begin
    sc=vehicle(); panel=SolarPanelAngleOfAttackControlModel((2,3))
    planet=make_no_gram_planet(:earth)
    position=SVector(planet.Rp_e+100e3,0.,0.); velocity=SVector(0.,7500.,0.)
    # At lat=lon=0 east is +y, north +z and up +x.
    env=EnvironmentSample(planet;planet_frame=PlanetFrameSample(SMatrix{3,3,Float64}(I),
        position,velocity,100e3,0.,0.),atmosphere=AtmosphereSample(1e-8,200.,SVector(100.,20.,3.)))
    sample=StateSample(position,velocity,60.;spacecraft=sc)
    model=AerodynamicCommandedAreaDragModel()
    relative=SVector(-3.,7400.,-20.)
    for alpha in (0.,pi/2)
        apply_panel_angle_of_attack!(panel,sc,alpha,0.)
        force,torque=wrench(model,sample,env,0.)
        area=sc.links[1].ref_area+sum(l.ref_area for l in sc.links[2:3])*sin(alpha)
        expected=-0.5*1e-8*2.2*area*norm(relative)*relative
        @test force ≈ expected rtol=1e-14
        @test torque == zeros(3)
        @test dot(force,relative)<0
    end
    @test environment_requirements(model).atmosphere
    @test environment_requirements(model).planet_frame
    @test solver_partition(model) == :implicit
    @test_throws ArgumentError apply_panel_angle_of_attack!(SolarPanelAngleOfAttackControlModel((1,)),sc,0.,0.)
    @test_throws ArgumentError apply_panel_angle_of_attack!(SolarPanelAngleOfAttackControlModel((4,)),sc,0.,0.)
    hinge=PanelHinge(axis=SVector(0.,1.,0.),min_angle_rad=-1.,max_angle_rad=1.,max_rate_rad_s=0.1,angle_rad=0.2)
    add_actuator!(sc,sc.links[2],hinge)
    apply_panel_angle_of_attack!(panel,sc,0.,0.)
    apply_panel_angle_of_attack!(panel,sc,0.,2.)
    @test hinge.angle_rad ≈ 0.4
    @test sc.links[2].α ≈ pi/2-0.4
    apply_panel_angle_of_attack!(panel,sc,0.,100.)
    @test hinge.angle_rad == 1.
    @test sc.links[2].α ≈ pi/2-1.
    @test_throws ArgumentError apply_panel_angle_of_attack!(panel,sc,0.,NaN)
end

@testset "Public orbital and geodetic conversion" begin
    planet=make_no_gram_planet(:earth)
    for (position,latitude) in (([planet.Rp_e+100e3,0.,0.],0.),([0.,0.,planet.Rp_p+100e3],pi/2))
        h,lat,lon=planet_fixed_to_geodetic(planet,position)
        @test h ≈ 100e3 atol=1e-6
        @test lat ≈ latitude atol=1e-12
    end
    phi=0.7; lambda=0.3; h=123e3
    e2=1-(planet.Rp_p/planet.Rp_e)^2; N=planet.Rp_e/sqrt(1-e2*sin(phi)^2)
    position=[(N+h)*cos(phi)*cos(lambda),(N+h)*cos(phi)*sin(lambda),(N*(1-e2)+h)*sin(phi)]
    actual=planet_fixed_to_geodetic(planet,position)
    @test actual[1] ≈ h atol=1e-5
    # Latitude budget corresponds to 1 mm on the reference surface; the
    # existing ellipsoid conversion is an approximation, not an exact inverse.
    @test actual[2] ≈ phi atol=1e-3/planet.Rp_e
    @test actual[3] ≈ lambda atol=1e-11
    ic=InitialCondition(ra=7e6,rp=7e6,i=0.,Ω=0.,ω=0.,ν=0.)
    r,v=orbital_elements_to_cartesian(ic,planet)
    @test r ≈ [7e6,0.,0.]
    @test v ≈ [0.,sqrt(planet.μ/7e6),0.]
end
end
