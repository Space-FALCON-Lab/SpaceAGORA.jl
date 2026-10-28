include("lifecycle-fixtures.jl")
include("lifecycle-reference.jl")
using Test,TOML
coincident=length(ARGS)>0 && ARGS[1]=="coincident"
prefix=coincident ? "coincident" : "main"
base1=length(ARGS)>=2 ? abspath(ARGS[2]) : joinpath(@__DIR__,"lifecycle-attempt1")
base2=length(ARGS)>=3 ? abspath(ARGS[3]) : joinpath(@__DIR__,"lifecycle-attempt2")
direct=deserialize(joinpath(base1,"$(prefix)-ks","ks-direct-result.bin"))
resume=deserialize(joinpath(base2,"$(prefix)-resume-ks","ks-resume-result.bin"))
audit=deserialize(joinpath(base2,"$(prefix)-resume-ks","ks-interrupted-audit.bin"))
checkpoint=coincident ? nothing : deserialize(joinpath(base2,"main-checkpoint-ks","ks-checkpoint-result.bin"))
ref=lifecycle_reference(;coincident)
start=coincident ? 221. : 997.3
duration=coincident ? 442. : 6000.
times=collect(0.:17.:duration);last(times)==duration || push!(times,duration)
normal_checkpoint=coincident ? nothing : deserialize(joinpath(base2,"main-checkpoint-dp8","dp8-checkpoint-result.bin"))
records=Dict{String,Any}()
@testset "Checkpoint/resume lifecycle $prefix" begin
    @test audit.checkpoint_time==start
    initial=first(resume.initial_audits)
    @test initial.epoch==start
    @test initial.root_alias && resume.root_alias
    @test initial.area==(coincident ? 10. : 15.)
    expected_tick=coincident ? 221. : 988.
    @test initial.last_tick==expected_tick
    @test initial.held_force_N≈.02*sin(expected_tick/1000.) atol=1e-14 rtol=0.
    @test initial.burn_active==(coincident ? Bool[] : [true])
    @test [last(e) for e in initial.committed_actions]==(coincident ? [:impulse] : [:area15])
    if coincident
        model=audit.checkpoint_snapshot.runtime_state.models.control[1]
        @test count(==(221.),model.ticks)==1
        @test model.actions==[(221.,:impulse)]
        order=[label for (t,label) in model.trace if t==221.]
        @test count(==(:impulse),order)==1 && count(==(:control),order)==1
        @test findfirst(==(:impulse),order)<findfirst(==(:control),order)<findfirst(==(:output),order)
        @test model.held_force_N≈.02*sin(221/1000.) atol=1e-14 rtol=0.
    end
    for (name,result) in (("resume",resume),("checkpoint",checkpoint),("normal_checkpoint",normal_checkpoint))
        result===nothing && continue
        @test result.ticks==direct.ticks
        @test result.actions==direct.actions
        @test all(clock->clock[2]==0.,result.clocks)
        @test length(result.roots)==length(direct.roots)
        for (a,b) in zip(result.roots,direct.roots)
            @test abs(a[1]-b[1])<=1e-4
        end
        dp=norm(result.final[1:3]-direct.final[1:3]);dv=norm(result.final[4:6]-direct.final[4:6])
        @test dp<=.1
        @test dv<=.0001
        @test abs(result.final[7]-direct.final[7])<=1e-7
        @test abs(result.final[8]-direct.final[8])<=max(1e-6,1e-5*abs(direct.final[8]))
        @test propertynames(result.table)==propertynames(direct.table)
        @test length(unique(result.table.time))==length(result.table.time)
        # Current upstream restores a valid checkpoint output prefix.
        # This is the normal results-table contract, not an accuracy relaxation.
        @test first(result.table.time)==0.
        if name=="resume"
            @test start in result.table.time
            @test all(t -> t in result.table.time, filter(<(start),times))
        end
        expected=times
        missing=setdiff(expected,result.table.time)
        @test isempty(missing)
        pe=0.;ve=0.;me=0.;qe=0.
        for (i,t) in enumerate(result.table.time)
            u=lifecycle_table_state(result.table,i);v=ref.at(t)
            pe=max(pe,norm(u[1:3]-v[1:3]));ve=max(ve,norm(u[4:6]-v[4:6]))
            me=max(me,abs(u[7]-v[7]));qe=max(qe,abs(u[8]-v[8])/max(1e-6,1e-5*abs(v[8])))
        end
        @test pe<=1.
        @test ve<=.001
        @test me<=1e-7
        @test qe<=1.
        records[name]=Dict("final_difference_position_m"=>dp,"final_difference_velocity_m_s"=>dv,
            "max_independent_position_error_m"=>pe,"max_independent_velocity_error_m_s"=>ve,
            "max_mass_error_kg"=>me,"max_heat_budget_ratio"=>qe,"missing_output_epochs"=>missing,
            "rows"=>size(result.table,1),"ticks"=>length(result.ticks))
        println("LIFECYCLE_CHECKPOINT_RESULT ",prefix," ",name," ",records[name]);flush(stdout)
    end
end
open(joinpath(base2,"$(prefix)-checkpoint-verdict.toml"),"w") do io
    TOML.print(io,Dict("case"=>prefix,"status"=>"passed frozen checkpoint/resume checks","methods"=>records))
end
