using Test
include("../../Output_6_Close_Form_Solution/oracle_closed_forms.jl")

@testset "Fixed-offset sparse-ring RTN prediction" begin
    p = Params()
    b0 = r_h(p)*deg2rad(1)
    run(inc; kw...) = leapfrog(100, b0; p, inc_deg=inc, u0_deg=1, kw...)

    # With the same trailing gap, reflecting the plane interchanges the two
    # normal impulses, without changing range eligibility or in-plane dynamics.
    up, down, flat = run(1), run(-1), run(0)
    @test up.dvNp > 0 > up.dvNm
    @test up.dvNp ≈ -down.dvNm && up.dvNm ≈ -down.dvNp
    @test up.dvT == down.dvT && up.dvR == down.dvR > 0
    @test flat.dvNp == flat.dvNm == 0
    @test up.links > 0 && flat.links > 0

    # A short first contact near peak latitude has a known initial sight line.
    seconds = 0.1
    t0 = leapfrog(1, b0; p, inc_deg=0.5, u0_deg=90,
        n_orbits=seconds*mean_motion(r_h(p), p)/(2pi))
    z = r_h(p)*sind(0.5)
    impulse = seconds*a_OCL(p)/hypot(b0, z)
    @test t0.dvT ≈ impulse*b0 rtol=1e-7
    @test t0.dvNp ≈ impulse*z rtol=1e-7
    @test t0.dvNm == 0 && t0.links == 0

    # Integrating da/dt = 2 A_T/n(a) gives this independent impulse check.
    rh = r_h(p)
    @test up.dvT ≈ sqrt(p.mu/rh) - sqrt(p.mu/(rh + 1000*up.da_km)) atol=1e-7
    finer = run(1; dt=0.5)
    @test all(f -> abs(up[f]-finer[f]) < 3e-4, (:dvR, :dvT, :dvNp, :dvNm))

    # Inside the along-track limit but outside the 3D range near peak latitude.
    off = leapfrog(100, 180e3; p, inc_deg=1, u0_deg=90, n_orbits=0.01)
    @test all(iszero, values(off))
    @test all(iszero, values(leapfrog(100, p.L_max+1; p, inc_deg=1)))
    @test_throws ArgumentError leapfrog(250, b0; p, inc_deg=1)
    @test_throws ArgumentError run(1; dt=0)
    @test propertynames(leapfrog(1, b0; p)) == (:dvT, :da_km, :links)
end
