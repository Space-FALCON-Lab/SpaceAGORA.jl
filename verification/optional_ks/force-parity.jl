ENV["SPACEAGORA_DENSITY_FREEZE_PER_STEP"]="false"
ENV["SPACEAGORA_GRAM_TRACK_CACHE_MODE"]="off"
include("baseline-fixtures.jl")
using Test
function run_case(case)
    function factory(args)
        function init(cb,u,t,int)
            for scale in (1.,1.0001,.9999,1.)
                trial=copy(u);trial.sc[1].pos .*=scale
                # Deliberately revisit the identical physical epoch after a
                # different state, which time-only sample caches cannot key.
                got=zero(trial);SE._trial_cartesian_rhs!(got,trial,int.p,Float64(t))
                expected=zeros(8);full_reference_rhs!(expected,Vector(trial),case,t)
                accelerror=norm(got[4:6]-expected[4:6]);heaterror=abs(got[8]-expected[8])
                println("FORCE_PARITY ",case," scale=",scale," accel_error_m_s2=",accelerror,
                    " heat_error_W_cm2=",heaterror," heat=",got[8]);flush(stdout)
                @test norm(got[1:3]-expected[1:3])<=1e-11
                @test accelerror<=1e-11
                @test heaterror<=1e-12
                @test got[7]==0.
            end
            SE._trial_cartesian_rhs!(zero(u),u,int.p,Float64(t))
        end
        SE.DiscreteCallback((u,t,int)->false,int->nothing;initialize=init,save_positions=(false,false))
    end
    sol=run_simulation(baseline_config(case,:dp8;duration=1.,rate=.1);return_solution=true,extra_callback_factories=(factory,))
    @test sol.t[end]==1.
end
@testset "Actual configured Cartesian RHS vs independent physical equations" begin
    for case in (:two_body,:J2_only,:drag_only,:J2_plus_drag)
        run_case(case)
    end
end
