# Shared core-only two-body fixture; no tests or cross-suite module dependency.
function fixture(; h=2.0^-12, duration=.01, position=SVector(7e6,0.,0.))
    SM=SpaceAGORA.SimulationModel
    bus = Link(root=true,m=20.,dims=MVector(1.,.8,.6))
    node = CompliantTopologyNode(:bumper; mass_kg=1.,inertia_body_kg_m2=.004,
        position=SVector(.12,0.,0.),velocity=SVector(-.25,0.,0.))
    edge = CompliantTopologyEdge(:mount,0,1; parent_point_body=Z, child_point_body=SVector(-.11,0.,0.),
        k_translation_n_m=200.,c_translation_n_s_m=.2,k_rotation_n_m_rad=2.,c_rotation_n_m_s_rad=.02)
    attachment = CompliantAttachment(model=build_compliant_topology([node],[edge]),link=bus,mount_point=SVector(.5,0.,0.))
    sc = SpacecraftModel(links=[bus],root=bus,prop_mass=0.,inertia_tensor=bus.inertia,
        initial_condition=CartesianInitialCondition(position,SVector(.25/21,0.,0.)),attachments=[attachment])
    base = make_example_config(planet=SM.make_no_gram_planet(:earth),spacecraft=sc,
        mission_time=duration,initial_time=InitialTime(year=2020,month=1,day=1),
        dynamic_effectors=(),ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=true,keplerian=false,verbose=false,results=false,
        solver_config=SolverConfig(solver_mode=:tsit5))
    return SM.SimConfig._with_configuration(base;
        mission_configuration=SM.SimConfig._with_fields(base.mission_configuration; data_rate=duration/2),
        integration_tolerances=IntegrationTolerances(dt_max_orbit=h,dt_max_atmosphere=h,
            reltol_orbit=1e-10,abstol_orbit=1e-11,reltol_quaternion=1e-10,abstol_quaternion=1e-11))
end
