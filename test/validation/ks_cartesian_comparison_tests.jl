using Test, LinearAlgebra
include(joinpath(@__DIR__,"..","..","examples","ks_cartesian_comparison","comparison.jl"))
const NC=NativeCoordinateComparison
const CS=NC.CS
const KS=NC.KS
BLAS.set_num_threads(1)

@testset "Native KS/Cartesian comparison" begin
    for e in (.01,.3)
        c=NC.ten_initial(e);y=NC.to_ks(c.x,c.planet)
        @test NC.to_cart(y) ≈ c.x rtol=1e-12
        dx=zeros(60);dy=zeros(100)
        NC.cart_rhs!(dx,vec(c.x),c.planet,0.)
        NC.ks_rhs!(dy,y,c.planet,0.)
        rho=sum(abs2,y[1:4])
        @test dy[10:10:end] ≈ fill(rho,10)
        for j in 1:10
            m=KS.ks_physical_tangent(y[10j-9:10j],c.planet)
            direct=CS.cartesian_model(c.x[1:3,j],c.x[4:6,j],c.planet)
            @test m.H*dy[10j-9:10j]/rho ≈ dx[6j-5:6j] rtol=1e-11 atol=1e-10
            @test m.A_cartesian ≈ direct.A_cartesian rtol=1e-11 atol=1e-12
            @test m.B_cartesian ≈ direct.B_cartesian atol=1e-12
        end
        for alg in (NC.Tsit5(),NC.DP5())
            sc=NC.solve_nonlinear(c,:cartesian,100.,alg;tol=1e-12)
            sk=NC.solve_nonlinear(c,:ks_target_s,100.,alg;tol=1e-12)
            for (sol,route) in ((sc,:cartesian),(sk,:ks_target_s))
                times=[0.,50.,100.]
                old=cat([NC.physical_sample(sol,route,t) for t in times]...;dims=3)
                @test NC.sample_trajectory(sol,route,times,10)==old
            end
            @test maximum(abs.(sk.u[end][10:10:end].-100.))<1e-6
            for t in (0.,50.,100.)
                @test maximum(abs.(NC.physical_sample(sc,:cartesian,t)-NC.physical_sample(sk,:ks_target_s,t)))<1e-3
            end
        end
    end
    case=NC.walker27_constellation();times=collect(range(0.,.02case.period_s;length=6))
    for objective in (:a,:i,:omega,:combined)
        a=NC.plan_laser_schedule(case,times;objective)
        b=NC.plan_laser_schedule(case,times;objective,reference=a.reference,linearization=:cartesian)
        @test a.command==b.command
        @test a.predicted≈b.predicted rtol=1e-9 atol=1e-8
        native=NC.inertial_predictions(case,times,a)
        @test all(isfinite,native.ks_native_s)
        @test maximum(abs.(native.ks_clock_perturbation_s))<1e-12
        cart=NC.cart_replay(case,times,a.command)
        @test length(NC.cart_elements(cart,last(times),case))==5
        ks=NC.replay_laser_schedule(case,times,a.command)
        @test maximum(abs.(reshape(cart(last(times)),6,:)-NC.to_cart(ks(last(times)))))<1e-3
    end
end
