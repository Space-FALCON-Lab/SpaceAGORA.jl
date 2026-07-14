function make_spacecraft(;
    ra_alt_m::Float64,
    rp_alt_m::Float64,
    i_deg::Float64=35.0,
    ω_deg::Float64=40.0,
    Ω_deg::Float64=10.0,
    ν_deg::Float64=175.0,
    orientation_state::Union{Nothing, Tuple{SVector{4, Float64}, SVector{3, Float64}}}=nothing
)
    root = Link{0}(root=true, m=500.0, ref_area=12.0)
    panel = Link{0}(root=false, m=30.0, ref_area=6.0, r=MVector{3, Float64}(0.0, 1.2, 0.0))

    if isnothing(orientation_state)
        ic = InitialCondition(
            ra=EARTH.Rp_e + ra_alt_m,
            rp=EARTH.Rp_e + rp_alt_m,
            i=i_deg,
            ω=ω_deg,
            Ω=Ω_deg,
            ν=ν_deg
        )
    else
        q0, w0 = orientation_state
        ra = EARTH.Rp_e + ra_alt_m
        rp = EARTH.Rp_e + rp_alt_m
        a = (ra + rp) / 2.0
        e = (ra - rp) / (ra + rp)
        ic = InitialCondition(a, e, i_deg, ω_deg, Ω_deg, ν_deg, q0, w0)
    end

    dry_mass = root.m + panel.m
    return SpacecraftModel(
        Joint[],
        [root, panel],
        root,
        true,
        dry_mass,
        0.0,
        root.inertia,
        0,
        0,
        ic,
        1
    )
end

function make_single_link_spacecraft(;
    ra_alt_m::Float64,
    rp_alt_m::Float64,
    i_deg::Float64=35.0,
    ω_deg::Float64=40.0,
    Ω_deg::Float64=10.0,
    ν_deg::Float64=175.0
)
    root = Link{0}(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(
        ra=EARTH.Rp_e + ra_alt_m,
        rp=EARTH.Rp_e + rp_alt_m,
        i=i_deg,
        ω=ω_deg,
        Ω=Ω_deg,
        ν=ν_deg
    )

    return SpacecraftModel(
        Joint[],
        [root],
        root,
        true,
        root.m,
        0.0,
        root.inertia,
        0,
        0,
        ic,
        1
    )
end

function make_base_thruster_model(;
    thrust::Real,
    direction::Real=0.0,
    Δv::Real=0.0,
    start_burn_time::Real=0.0,
    stop_burn_time::Real=0.0,
    Isp::Real=300.0
)
    return BaseThrusterModel(
        thrust=[Float64(thrust)],
        direction=[Float64(direction)],
        Δv=[Float64(Δv)],
        start_burn_time=[Float64(start_burn_time)],
        stop_burn_time=[Float64(stop_burn_time)],
        Isp=[Float64(Isp)]
    )
end

function make_agora_earth_spacecraft()
    main_bus = Link{0}(root=true, m=620.0, ref_area=2.05 * 2.8)
    left_panel = Link{0}(root=false, m=10.0, ref_area=5.7 * 1.0 / 2.0, r=MVector{3, Float64}(0.0, -2.05 / 2.0 - 5.7 / 4.0, 0.0))
    right_panel = Link{0}(root=false, m=10.0, ref_area=5.7 * 1.0 / 2.0, r=MVector{3, Float64}(0.0, 2.05 / 2.0 + 5.7 / 4.0, 0.0))
    ic = InitialCondition(
        ra=56_378.7978559e3,
        rp=EARTH.Rp_e + 200_590.0,
        i=89.876,
        ω=75.505,
        Ω=104.115,
        ν=175.0
    )

    return SpacecraftModel(
        Joint[],
        [main_bus, left_panel, right_panel],
        main_bus,
        true,
        main_bus.m + left_panel.m + right_panel.m,
        200.0,
        main_bus.inertia,
        0,
        0,
        ic,
        1
    )
end
