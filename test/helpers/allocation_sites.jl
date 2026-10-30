# Temporary diagnostic: record every allocation of one call with the allocation
# profiler and print the distinct sites (type, size, stack), so a CI log shows
# where an "allocation-free" right-hand side allocates on the CI platform.
# Profile is loaded only when the first report runs, after the test's own
# measurement, so that loading it cannot change which code the measurement sees.
const _PROFILE = Ref{Module}()
function _profile_module()
    isassigned(_PROFILE) && return _PROFILE[]
    _PROFILE[] = Base.require(Base.PkgId(Base.UUID("9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"), "Profile"))
    return _PROFILE[]
end

function _report_allocation_sites(f::F, label::AbstractString) where {F}
    Profile = _profile_module()
    f()                                   # the closure itself is compiled before the profile
    Base.invokelatest(Profile.Allocs.clear)
    Base.invokelatest(Profile.Allocs.start; sample_rate=1)
    try
        f()
    finally
        Base.invokelatest(Profile.Allocs.stop)
    end
    results = Base.invokelatest(Profile.Allocs.fetch)
    sites = Dict{Tuple{Any, Int, String}, Int}()
    for a in results.allocs
        frames = String[]
        for fr in a.stacktrace
            name = string(fr.func)
            file = basename(string(fr.file))
            (startswith(name, "_report_allocation_sites") || startswith(name, "jl_") || file in ("Allocs.jl", "boot.jl")) && continue
            push!(frames, string(name, fr.inlined ? "[inl]" : "", " @ ", file, ":", fr.line))
            length(frames) >= 18 && break
        end
        key = (a.type, a.size, join(frames, "\n        <- "))
        sites[key] = get(sites, key, 0) + 1
    end
    total = sum(a.size for a in results.allocs; init=0)
    println("---- allocation sites: ", label, " (", length(results.allocs), " allocations, ", total, " bytes, pkgimages=",
        Base.JLOptions().use_pkgimages, ", threads=", Threads.nthreads(), ")")
    for ((T, sz, fr), n) in sort(collect(sites); by = x -> (-x[2] * x[1][2], x[1][3]))
        println("  ", n, " x ", sz, " B  ", T, "\n        ", fr)
    end
    flush(stdout)
    return nothing
end
