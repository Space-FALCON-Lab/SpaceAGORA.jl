#!/usr/bin/env julia

function usage(io::IO=stdout)
    println(io, "Usage: julia --project=. scripts/ensure_gram_native.jl [--clean] [GRAM_ROOT]")
    println(io)
    println(io, "Ensure a native GRAM shared library exists for the current operating system.")
    println(io)
    println(io, "Arguments:")
    println(io, "  --clean       force a clean rebuild")
    println(io, "  GRAM_ROOT     optional GRAM Suite root; defaults to vendored GRAM tree")
end

function expected_lib_ext()
    if Sys.isapple()
        return "dylib"
    elseif Sys.iswindows()
        return "dll"
    end
    return "so"
end

function normalize_gram_root(path::AbstractString)
    expanded = expanduser(path)
    isdir(expanded) || error("GRAM root not found: $(abspath(expanded))")
    return realpath(expanded)
end

function build_command(build_helper::AbstractString, clean::Bool, gram_root::AbstractString)
    args = clean ? ["--clean", gram_root] : [gram_root]
    if Sys.iswindows()
        cmd_helper = replace(build_helper, '/' => '\\')
        return Cmd(["cmd", "/C", cmd_helper, args...])
    end
    return Cmd(["bash", build_helper, args...])
end

function main(args::Vector{String})
    clean = false
    gram_root_arg = nothing

    for arg in args
        if arg == "--clean"
            clean = true
        elseif arg == "-h" || arg == "--help"
            usage()
            return 0
        elseif isnothing(gram_root_arg)
            gram_root_arg = arg
        else
            usage(stderr)
            println(stderr)
            println(stderr, "Only one GRAM_ROOT argument is allowed.")
            return 1
        end
    end

    repo_root = normpath(joinpath(@__DIR__, ".."))
    default_gram_root = joinpath(repo_root, "data", "GRAMSuite.jl", "GRAM Suite 2.0")
    gram_root = normalize_gram_root(something(gram_root_arg, get(ENV, "GRAM_ROOT", default_gram_root)))

    build_helper = Sys.iswindows() ?
        joinpath(gram_root, "simulation", "GRAM", "build_gram.cmd") :
        joinpath(gram_root, "simulation", "GRAM", "build_gram.sh")
    isfile(build_helper) || error("GRAM build helper not found: $build_helper")

    expected_lib = joinpath(gram_root, "Build", "lib", "libGRAM.$(expected_lib_ext())")

    println("GRAM root:        ", gram_root)
    println("Expected library: ", expected_lib)

    if isfile(expected_lib) && !clean
        println("Native GRAM library already present for this host.")
        return 0
    end

    println("Native GRAM library missing for this host. Building now...")
    run(build_command(build_helper, clean, gram_root))

    isfile(expected_lib) || error("Build completed but expected native library is still missing: $expected_lib")

    println("Native GRAM library ready:")
    println("  ", expected_lib)
    return 0
end

try
    exit(main(copy(ARGS)))
catch err
    println(stderr, err)
    exit(1)
end
