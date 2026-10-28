for (i,name) in enumerate(("upstream-data-guard.jl","upstream-scope-guard.jl","force-parity.jl","optional-guards.jl","dense-history.jl"))
    println("BEGIN_ADAPTER_CHECK ",name);flush(stdout)
    symbol=Symbol("AdapterCheck",i)
    path=joinpath(@__DIR__,name)
    Core.eval(Main, :(module $symbol
        include($path)
    end))
end
println("ADAPTER_CURRENT_SOURCE_CHECKS_SUCCESS")
