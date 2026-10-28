include("lifecycle-fixtures.jl")
include("lifecycle-reference.jl")
using Test, TOML
coincident=length(ARGS)>0 && ARGS[1]=="coincident"
base=length(ARGS)>=2 ? abspath(ARGS[2]) : joinpath(@__DIR__,"lifecycle-attempt1")
prefix=coincident ? "coincident" : "main"
normal=deserialize(joinpath(base,"$(prefix)-dp8","dp8-direct-result.bin"))
ks=deserialize(joinpath(base,"$(prefix)-ks","ks-direct-result.bin"))
duration=coincident ? 442. : 6000.
times=collect(0.:17.:duration)
last(times)==duration || push!(times,duration)
println("LIFECYCLE_REFERENCE_BEGIN ",prefix);flush(stdout)
coarse=lifecycle_reference(;coincident,rtol=1e-12,atol=1e-12,dtmax=10.)
reference=lifecycle_reference(;coincident,rtol=2.5e-13,atol=2.5e-13,dtmax=5.)
function compare(values,truth)
    position=0.;velocity=0.;mass=0.;heat=0.;heatbudget=0.
    for t in times
        a=values(t);b=truth(t)
        position=max(position,norm(a[1:3]-b[1:3]));velocity=max(velocity,norm(a[4:6]-b[4:6]))
        mass=max(mass,abs(a[7]-b[7]));heat=max(heat,abs(a[8]-b[8]))
        heatbudget=max(heatbudget,abs(a[8]-b[8])/max(1e-6,1e-5*abs(b[8])))
    end
    return Dict("position_m"=>position,"velocity_m_s"=>velocity,"mass_kg"=>mass,"heat_J_cm2"=>heat,"heat_budget_ratio"=>heatbudget)
end
referror=compare(coarse.at,reference.at)
records=Dict{String,Any}()
@testset "Direct lifecycle $prefix" begin
    @test referror["position_m"]<=.01
    @test referror["velocity_m_s"]<=1e-5
    @test referror["mass_kg"]<=1e-9
    @test referror["heat_budget_ratio"]<=.01
    @test propertynames(normal.table)==propertynames(ks.table)
    for result in (normal,ks)
        @test result.table.time==times
        @test result.ticks==collect(13.:13.:prevfloat(duration))
        @test result.root_alias
        @test all(clock->clock[2]==0.,result.clocks)
        @test [last(x) for x in result.actions]==(coincident ? [:impulse] : [:area15,:impulse,:area8])
        planned=coincident ? [221.] : [611.25,1879.125,3599.75]
        @test maximum(abs.([first(x) for x in result.actions]-planned))<=1e-8
        lookup=Dict(t=>i for (i,t) in enumerate(result.table.time))
        at(t)=lifecycle_table_state(result.table,lookup[t])
        m=compare(at,reference.at)
        @test m["position_m"]<=1.
        @test m["velocity_m_s"]<=.001
        @test m["mass_kg"]<=1e-7
        @test m["heat_budget_ratio"]<=1.
        @test maximum(abs(at(t)[7]-(coincident ? 500. : lifecycle_mass_oracle(t))) for t in times)<=1e-7
        @test length(result.roots)==length(reference.roots)==(coincident ? 0 : 1)
        for (found,expected) in zip(result.roots,reference.roots)
            @test abs(found[1]-expected)<=1e-4
            @test abs(found[4])<=1e-4  # stored[t,x,y,z,...]
            @test found[7]<0.        # stored vz; descending crossing
        end
        if coincident
            order=[label for (t,label) in result.trace if t==221.]
            @test order==[:impulse,:control,:output]
            @test result.table.lifecycle_held_force_N[lookup[221.]]≈.02*sin(221/1000.) atol=1e-14 rtol=0.
        end
        # Force-output columns are separate diagnostics, not post-hoc gates.
        heat_rate_diagnostic=maximum(abs(result.table.sc1_heat_rate[i]-reference_heat(reference.at(t))) for (i,t) in enumerate(times))
        records[string(result.mode)]=Dict("errors"=>m,"heat_rate_output_max_difference_W_cm2"=>heat_rate_diagnostic,
            "ticks"=>length(result.ticks),"rows"=>size(result.table,1),"roots"=>[r[1] for r in result.roots])
        println("LIFECYCLE_DIRECT_RESULT ",result.mode," ",m," force_output_heat_rate_diagnostic=",heat_rate_diagnostic);flush(stdout)
    end
end
open(joinpath(base,"$(prefix)-direct-verdict.toml"),"w") do io
    TOML.print(io,Dict("case"=>prefix,"status"=>"passed direct lifecycle checks; checkpoint/resume separate", "reference_refinement"=>referror,"methods"=>records))
end
