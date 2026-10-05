using Test

const WORKSPACE_ROOT = normpath(joinpath(@__DIR__, ".."))

function static_path(expression, file)
    expression isa String && return expression
    expression isa Expr || return nothing
    if expression.head === :macrocall && expression.args[1] === Symbol("@__DIR__")
        return dirname(file)
    elseif expression.head === :macrocall && expression.args[1] === Symbol("@__FILE__")
        return file
    elseif expression.head === :call && expression.args[1] in (:joinpath, :normpath, :abspath, :dirname)
        arguments = [static_path(argument, file) for argument in expression.args[2:end]]
        any(isnothing, arguments) && return nothing
        operation = expression.args[1]
        operation === :joinpath && return joinpath(arguments...)
        operation === :normpath && return normpath(only(arguments))
        operation === :dirname && return dirname(only(arguments))
        operation === :abspath && return abspath(dirname(file), arguments...)
    end
    return nothing
end

function inspect_syntax(expression, file, includes, errors)
    expression isa Expr || return nothing
    if expression.head in (:error, :incomplete)
        push!(errors, (file, expression))
    elseif expression.head === :call && expression.args[1] === :include && length(expression.args) == 2
        target = static_path(expression.args[2], file)
        target === nothing || push!(includes, (file, normpath(joinpath(dirname(file), target))))
    end
    for argument in expression.args
        inspect_syntax(argument, file, includes, errors)
    end
    return nothing
end

@testset "Repository Julia syntax and static include paths" begin
    files = split(read(`git -C $WORKSPACE_ROOT ls-files -z -- '*.jl'`, String), '\0'; keepempty=false)
    push!(files, relpath(@__FILE__, WORKSPACE_ROOT))
    sort!(unique!(files))
    includes = Tuple{String, String}[]
    errors = Tuple{String, Expr}[]
    for relative in files
        file = joinpath(WORKSPACE_ROOT, relative)
        inspect_syntax(Meta.parseall(read(file, String); filename=file), file, includes, errors)
    end
    @test isempty(errors)
    for (file, target) in includes
        if !isfile(target) && occursin("/data/GRAMSuite.jl/", target)
            @test_skip isfile(target)
        else
            @test isfile(target)
            isfile(target) || @error "Missing include" source=file target=target
        end
    end
    projects = ("2_SpaceAGORA.jl", "3_SpaceAGORA.jl-main_blank", "4_SpaceAGORA.jl-main_v1",
        "5_SpaceAGORA.jl-main_v2", "6_SpaceAGORA.jl-main_v3")
    for project in projects
        contracts = joinpath(WORKSPACE_ROOT, project, "test", "contracts")
        pr = Set(match.captures[1] for match in eachmatch(r"\"([^\"]+_gate\.jl)\"", read(joinpath(contracts, "pr_runtests.jl"), String)))
        nightly = Set(match.captures[1] for match in eachmatch(r"\"([^\"]+_gate\.jl)\"", read(joinpath(contracts, "nightly_runtests.jl"), String)))
        @test issubset(nightly, pr)
        builder = joinpath(WORKSPACE_ROOT, project, "data", "GRAMSuite.jl", "GRAM Suite 2.0",
            "simulation", "GRAM", "build_offline_static_grids.jl")
        if isfile(builder)
            @test isfile(builder)
        else
            @info "Optional GRAM builder is unavailable; external-asset validation skipped" project=project
            @test_skip isfile(builder)
        end
    end
    println("Checked $(length(files)) Julia files and $(length(includes)) static includes.")
end
