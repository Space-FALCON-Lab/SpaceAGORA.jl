using Test

# Execute the real constructor/loading functions with an instrumented SPICE
# boundary. No external kernels or package initialization are needed: the test
# observes the ordering that determines binary-PCK precedence.
module PCKOrderFixture
    const _SPICE_BODY_LOCK = ReentrantLock()
    const _FURNISHED_KERNELS = Set{String}()
    const _SPICE_KERNEL_POOL_GENERATION = Ref{UInt}(0)
    const _PCK_OVERRIDE_PATHS = Ref{Union{Nothing, Vector{String}}}(nothing)
    const _PCK_APPLIED_KERNEL_COUNT = Ref(-1)
    const loaded = String[]
    const calls = String[]
    const _EARTH_CACHE = Dict{Tuple{String, String}, Any}()
    const _MOON_CACHE = Dict{Tuple{String, String}, Any}()
    const _MARS_CACHE = Dict{Tuple{String, String}, Any}()
    const _VENUS_CACHE = Dict{Tuple{String, String}, Any}()
    const _TITAN_CACHE = Dict{Tuple{String, String}, Any}()
    Earth(; kwargs...) = :earth
    Moon(; kwargs...) = :moon
    Mars(; kwargs...) = :mars
    Venus(; kwargs...) = :venus
    Titan(; kwargs...) = :titan
    _furnsh_planetary_kernel(_) = nothing
    _gravity_constants_kernel_if_available(dir) = _furnsh_required(dir, "pck/gm_de440.tpc")
    _furnsh_mars_pck(dir) = _furnsh_required(dir, "pck/pck00011.tpc")
    _furnsh_mars_system_kernel(dir) = _furnsh_required(dir, "spk/satellites/mar099.bsp")
    _spice_backed_planet_kwargs(_) = (;)
    function furnsh(path)
        push!(calls, path)
        push!(loaded, path)
    end
    module SPICE
        function unload(path)
            filter!(!=(path), parentmodule(@__MODULE__).loaded)
        end
    end
    unload(path) = SPICE.unload(path)
    source = Meta.parse(read(joinpath(@__DIR__, "..", "..", "..", "src", "models", "environment", "ephemerides", "planets.jl"), String))
    wanted = Set([:_reset_furnished_kernels!, :_furnsh_once, :_furnsh_required, :_furnsh_first_existing_if_available,
                  :_furnsh_first_existing, :_check_pck_override_policy!, :_furnsh_pck_overrides, :Earth, :Moon, :Mars, :Venus, :Titan])
    found = Set{Symbol}()
    for expr in source.args[3].args
        unwrapped = expr
        if unwrapped isa Expr && unwrapped.head == :macrocall && unwrapped.args[1] == Symbol("@inline")
            unwrapped = unwrapped.args[end]
        end
        unwrapped isa Expr && unwrapped.head == :function || continue
        signature = unwrapped.args[1]
        while signature isa Expr && signature.head in (:(::), :where)
            signature = signature.args[1]
        end
        signature isa Expr && signature.head == :call || continue
        name = signature.args[1]
        name in wanted || continue
        Core.eval(@__MODULE__, expr)
        push!(found, name)
    end
    found == wanted || error("constructor fixture missed production functions: $(setdiff(wanted, found))")
    function reset!()
        _reset_furnished_kernels!()
        empty!(loaded)
        empty!(calls)
    end
end

@testset "kernel pool generations follow changes and explicit reset" begin
    mktempdir() do dir
        PCKOrderFixture.reset!()
        first_kernel = joinpath(dir, "first.tpc")
        next_kernel = joinpath(dir, "next.tpc")
        override = joinpath(dir, "override.tpc")
        foreach(touch, (first_kernel, next_kernel, override))
        generation = PCKOrderFixture._SPICE_KERNEL_POOL_GENERATION[]

        PCKOrderFixture._furnsh_once(first_kernel)
        @test PCKOrderFixture._SPICE_KERNEL_POOL_GENERATION[] == generation + 1
        PCKOrderFixture._furnsh_once(first_kernel)
        @test PCKOrderFixture._SPICE_KERNEL_POOL_GENERATION[] == generation + 1
        @test PCKOrderFixture.loaded == [first_kernel]

        withenv("SPACEAGORA_SPICE_PCK_OVERRIDES" => override) do
            PCKOrderFixture._furnsh_pck_overrides(dir)
            @test PCKOrderFixture._SPICE_KERNEL_POOL_GENERATION[] == generation + 2
            PCKOrderFixture._furnsh_pck_overrides(dir)
            @test PCKOrderFixture._SPICE_KERNEL_POOL_GENERATION[] == generation + 2

            PCKOrderFixture._furnsh_once(next_kernel)
            @test PCKOrderFixture._SPICE_KERNEL_POOL_GENERATION[] == generation + 3
            # Reapplying the override changes the pool twice: unload, then load.
            PCKOrderFixture._furnsh_pck_overrides(dir)
            @test PCKOrderFixture._SPICE_KERNEL_POOL_GENERATION[] == generation + 5
            @test PCKOrderFixture.loaded == [first_kernel, next_kernel, override]
            @test length(PCKOrderFixture._FURNISHED_KERNELS) == 3
        end

        caches = (PCKOrderFixture._EARTH_CACHE, PCKOrderFixture._MOON_CACHE,
                  PCKOrderFixture._MARS_CACHE, PCKOrderFixture._VENUS_CACHE,
                  PCKOrderFixture._TITAN_CACHE)
        for cache in caches
            cache[("", dir)] = :cached_planet
        end
        PCKOrderFixture._reset_furnished_kernels!()
        @test PCKOrderFixture._SPICE_KERNEL_POOL_GENERATION[] == generation + 6
        @test isempty(PCKOrderFixture._FURNISHED_KERNELS)
        @test PCKOrderFixture._PCK_OVERRIDE_PATHS[] === nothing
        @test PCKOrderFixture._PCK_APPLIED_KERNEL_COUNT[] == -1
        @test all(isempty, caches)
        # An explicit reset also invalidates readers when the path set was empty.
        PCKOrderFixture._reset_furnished_kernels!()
        @test PCKOrderFixture._SPICE_KERNEL_POOL_GENERATION[] == generation + 7
    end
end

@testset "PCK overrides follow all defaults across constructor orders" begin
    mktempdir() do dir
        defaults = ["pck/pck00011.tpc", "lsk/naif0012.tls",
                    "pck/earth_latest_high_prec.bpc", "tf/earth_assoc_itrf93.tf",
                    "spk/satellites/SPICELunaCurrentKernel.bpc", "tf/SPICELunaFrameKernel.tf",
                    "pck/gm_de440.tpc", "spk/satellites/mar099.bsp", "spk/satellites/sat441.bsp"]
        overrides = ["replacement-a.bpc", "replacement-b.bpc"]
        for relative in [defaults; overrides]
            path = joinpath(dir, relative)
            mkpath(dirname(path))
            touch(path)
        end
        for order in ((:Earth, :Moon), (:Moon, :Earth), (:Mars, :Moon, :Earth),
                      (:Venus, :Earth, :Moon), (:Titan, :Moon, :Earth))
            PCKOrderFixture.reset!()
            withenv("SPACEAGORA_SPICE_PCK_OVERRIDES" => join(overrides, ",")) do
                for body in order
                    constructor = getfield(PCKOrderFixture, body)
                    constructor("", dir)
                    @test PCKOrderFixture.loaded[end-1:end] == joinpath.(dir, overrides)
                    @test count(==(joinpath(dir, overrides[1])), PCKOrderFixture.loaded) == 1
                    before = copy(PCKOrderFixture.calls)
                    constructor("", dir)
                    @test PCKOrderFixture.calls == before
                end
            end
        end
        PCKOrderFixture.reset!()
        withenv("SPACEAGORA_SPICE_PCK_OVERRIDES" => nothing) do
            for body in (:Earth, :Moon, :Mars, :Venus, :Titan)
                getfield(PCKOrderFixture, body)("", dir)
            end
            @test !any(path -> basename(path) in overrides, PCKOrderFixture.loaded)
        end
        before = copy(PCKOrderFixture.loaded)
        withenv("SPACEAGORA_SPICE_PCK_OVERRIDES" => "missing.bpc") do
            @test_throws ArgumentError PCKOrderFixture._furnsh_pck_overrides(dir)
        end
        @test PCKOrderFixture.loaded == before
    end
end
