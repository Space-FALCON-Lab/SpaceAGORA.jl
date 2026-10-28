# Normal package execution; caller activates integration/environment at process startup.
using SpaceAGORA, LinearAlgebra
const S=SpaceAGORA.SimulationModel
const SE=SpaceAGORA.SimulationEngine
const MU=3.986004418e14
const RE=6378136.3
const SPIN=7.292115e-5
const J2=0.00108262668
const SMA=7e6
const ECC=.01
const INC=.7
const NODE=.3
const PERI=.2
const INITIAL_MASS=500.
const RHO_REF=1e-11
const H_REF=621863.
const SCALE_HEIGHT=60000.
const CD_AREA=22.
const TEMP=200.
function orbit_initial()
    # Active rotation from perifocal frame to inertial axes.
    rz(x)=[cos(x) -sin(x) 0.;sin(x) cos(x) 0.;0. 0. 1.]
    rx(x)=[1. 0. 0.;0. cos(x) -sin(x);0. sin(x) cos(x)]
    Q=rz(NODE)*rx(INC)*rz(PERI)
    return Q*[SMA*(1-ECC),0.,0.],Q*[0.,sqrt(MU*(1+ECC)/(SMA*(1-ECC))),0.]
end
function baseline_config(case::Symbol,mode::Symbol;duration=86400.,rate=60.,rtol=1e-10,atol=1e-11,dtmax=10.)
    hasj2=case in (:J2_only,:J2_plus_drag)
    hasdrag=case in (:drag_only,:J2_plus_drag)
    case in (:two_body,:J2_only,:drag_only,:J2_plus_drag) || error("Unknown case")
    planet=S.Earth(μ=MU,Rp_e=RE,Rp_p=RE,Rp_m=RE,ω=S.SVector(0.,0.,SPIN),J2=J2,R=287.1,γ=1.4)
    root=S.Link(root=true,m=INITIAL_MASS,ref_area=10.,α=pi/2,β=0.,θ=0.)
    r,v=orbit_initial()
    sc=S.SpacecraftModel(joints=S.Joint[],links=[root],root=root,prop_mass=0.,inertia_tensor=root.inertia,
        initial_condition=S.CartesianInitialCondition(r,v),id=1)
    grav=hasj2 ? S.InverseSquaredJ2GravityModel() : S.InverseSquaredGravityModel()
    effects=hasdrag ? (grav,S.AerodynamicCoefficientConstant()) : (grav,)
    density=hasdrag ? S.ExponentialAtmosphereModel(RHO_REF,H_REF,SCALE_HEIGHT;temperature_k=TEMP,
        valid_min_altitude_m=0.,valid_max_altitude_m=1e6) : S.NoAtmosphereModel()
    return S.SimulationConfiguration(
        simulation_settings=S.SimulationSettings(results=false,verbose=false,generate_plots=false,
            results_directory=mktempdir(),normalize=false),
        mission_configuration=S.MissionConfiguration(mission_type=S.MissionTime,keplerian=true,
            number_of_orbits=1,mission_time=duration,orientation_sim=false,num_steps_to_save=100000,data_rate=rate),
        environment_model=S.EnvironmentModel(planet=planet,EI=1000.,density_model=density,
            ephemerides_model=S.SimpleEphemeridesModel(),thermal_model=S.MaxwellianHeat(
                thermal_accomodation_factor=1.,planet=planet,thermal_contact=false),topography=false,wind=false),
        dynamics_model=S.DynamicsModel([sc],effects),
        guidance_model=S.GuidanceModel(guidance_effectors=(),guidance_rates=Float64[]),
        navigation_model=S.NavigationModel(navigation_effectors=(),navigation_rates=Float64[]),
        control_model=S.ControlModel(control_effectors=(),control_rates=Float64[]),
        initial_time=S.InitialTime(year=2020,month=1,day=1,hour=0,minute=0,second=0.),
        integration_tolerances=S.IntegrationTolerances(reltol=rtol,abstol=atol,reltol_orbit=rtol,abstol_orbit=atol,
            reltol_atmosphere=rtol,abstol_atmosphere=atol,reltol_mass=rtol,abstol_mass=atol,
            reltol_heat_load=rtol,abstol_heat_load=atol,dt_max=dtmax,dt_max_orbit=dtmax,dt_max_atmosphere=dtmax),
        solver_config=S.SolverConfig(solver_mode=mode))
end
# Independent Cartesian equations. No KS transformations or SpaceAGORA force calls.
function reference_rhs!(dy,y,case,t)
    r=@view y[1:3];v=@view y[4:6]
    radius=norm(r);zratio=r[3]^2/radius^2
    dy[1:3].=v
    dy[4:6].=-MU.*r./radius^3
    if case in (:J2_only,:J2_plus_drag)
        k=1.5*J2*MU*RE^2/radius^5
        dy[4]+=k*r[1]*(5zratio-1)
        dy[5]+=k*r[2]*(5zratio-1)
        dy[6]+=k*r[3]*(5zratio-3)
    end
    if case in (:drag_only,:J2_plus_drag)
        rho=RHO_REF*exp(-(radius-RE-H_REF)/SCALE_HEIGHT)
        relative=v-cross([0.,0.,SPIN],r)
        dy[4:6].-=.5*rho*CD_AREA/INITIAL_MASS*norm(relative).*relative
    end
    return nothing
end
function reference_solution(case;duration=86400.,rtol=1e-12,atol=1e-12,dtmax=20.)
    r,v=orbit_initial()
    p=SE.ODEProblem(reference_rhs!,vcat(r,v),(0.,duration),case)
    return SE.solve(p,SE.DP8();reltol=rtol,abstol=atol,dtmax=dtmax,save_everystep=true,dense=true)
end
function analytic_kepler(t)
    n=sqrt(MU/SMA^3);M=n*t;E=M
    for _ in 1:16
        E-=(E-ECC*sin(E)-M)/(1-ECC*cos(E))
    end
    r0,v0=orbit_initial()
    phat=r0/norm(r0);qhat=v0/norm(v0)
    r=SMA*(cos(E)-ECC)*phat+SMA*sqrt(1-ECC^2)*sin(E)*qhat
    v=n*SMA/(1-ECC*cos(E))*(-sin(E)*phat+sqrt(1-ECC^2)*cos(E)*qhat)
    return vcat(r,v)
end
function reference_heat(y)
    r=@view y[1:3];v=@view y[4:6]
    rho=RHO_REF*exp(-(norm(r)-RE-H_REF)/SCALE_HEIGHT)
    speed=norm(v-cross([0.,0.,SPIN],r))
    gas=287.1;thermal_speed=sqrt(2gas*TEMP);ratio=speed/thermal_speed
    # Simplification at alpha=pi/2, accommodation=1, wall=ambient T.
    # This is the independently expressed RETURNED heat flux; unused recovery
    # temperature intermediates in SpaceAGORA do not enter this quadrature.
    erfratio=getfield(parentmodule(S.MaxwellianHeat),:erf)(ratio)
    return 1e-4*rho*gas*TEMP*sqrt(gas*TEMP/(2pi))*(ratio^2*exp(-ratio^2)+
        sqrt(pi)*ratio*(ratio^2+.5)*(1+erfratio))
end
function full_reference_rhs!(dy,y,case,t)
    reference_rhs!(dy,y,case,t)
    dy[7]=0.
    dy[8]=case in (:drag_only,:J2_plus_drag) ? reference_heat(y) : 0.
    return nothing
end
function full_reference_solution(case;duration=86400.,rtol=1e-12,atol=1e-12,dtmax=20.)
    r,v=orbit_initial()
    p=SE.ODEProblem(full_reference_rhs!,vcat(r,v,INITIAL_MASS,0.),(0.,duration),case)
    return SE.solve(p,SE.DP8();reltol=rtol,abstol=atol,dtmax=dtmax,save_everystep=true,dense=true)
end
