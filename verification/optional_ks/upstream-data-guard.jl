include("baseline-fixtures.jl")
using KSPropagator, Test
@testset "New upstream mounted data remain outside initial KS scope" begin
    for field in (:data_stores, :data_sources)
        cfg=baseline_config(:two_body,:ks;duration=2.)
        sc=only(cfg.dynamics_model.spacecraft)
        item=field===:data_stores ? S.DataStore(capacity_bits=1e6) : S.DataSource(rate_bps=1.)
        S.add_data!(sc,sc.root,item)
        err=try SE._validate_ks_config!(cfg);nothing catch e;e end
        println("MOUNTED_DATA_GUARD field=",field," error=",err);flush(stdout)
        @test err isa ArgumentError
        @test err isa ArgumentError && occursin(string(field),sprint(showerror,err))
    end
end
