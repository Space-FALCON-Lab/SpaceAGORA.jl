# build_sysimage.jl
#
# Builds a custom Julia system image that has all SpaceAGORA.jl packages
# pre-compiled into it.  Subsequent Julia sessions launched with
#
#   julia --project=. --sysimage SpaceAGORA.so <script.jl>
#
# start in seconds rather than minutes because they skip the per-package
# compilation step entirely.
#
# Usage (from the repository root):
#
#   julia --project=build build/build_sysimage.jl
#
# The resulting sysimage is written to SpaceAGORA.so in the repository root.
# That file is listed in .gitignore and should NOT be committed to version
# control.  Rebuild it whenever the package set in Project.toml
# changes.

using PackageCompiler

# All packages from the root Project.toml that should be baked into the
# sysimage.  PackageCompiler needs them listed by symbol.
const PACKAGES = [
    :Arrow,
    :AssociatedLegendrePolynomials,
    :AstroTime,
    :BenchmarkTools,
    :CSV,
    :ControlSystems,
    :ControlSystemsBase,
    :DataFrames,
    :DateFormats,
    :DiffEqBase,
    :DiffEqCallbacks,
    :DifferentialEquations,
    :Distributions,
    :Enzyme,
    :FiniteDiff,
    :ForwardDiff,
    :GRAMSuite,
    :Interpolations,
    :LaTeXStrings,
    :LegendrePolynomials,
    :LoopVectorization,
    :ManifoldsBase,
    :MatrixEquations,
    :ModelingToolkit,
    :NLsolve,
    :NonlinearSolve,
    :OrdinaryDiffEq,
    :PlotlyJS,
    :Plots,
    :PreallocationTools,
    :Profile,
    :PythonCall,
    :Quaternions,
    :Reexport,
    :Roots,
    :Rotations,
    :SPICE,
    :SatelliteToolbox,
    :SatelliteToolboxAtmosphericModels,
    :SatelliteToolboxGeomagneticField,
    :SatelliteToolboxTransformations,
    :SimpleDiffEq,
    :SparseArrays,
    :SpecialFunctions,
    :StaticArrays,
]

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))

# The sysimage is written next to this build script's parent (i.e. repo root).
const SYSIMAGE_PATH = joinpath(PROJECT_ROOT, "SpaceAGORA.so")

# Optional: provide a precompile script so PackageCompiler can record which
# methods are actually called and include their native code.
const PRECOMPILE_SCRIPT = let
    override = get(ENV, "SPACEAGORA_SYSIMAGE_PRECOMPILE_FILE", "")
    default_script = joinpath(PROJECT_ROOT, "src", "examples", "Earth.jl")

    if !isempty(override)
        override
    elseif get(ENV, "CI", "false") == "true"
        nothing
    elseif isfile(default_script)
        default_script
    else
        nothing
    end
end

println("Building sysimage → $(abspath(SYSIMAGE_PATH))")
println("This will take several minutes the first time …")
if PRECOMPILE_SCRIPT === nothing
    println("Precompile execution trace: skipped")
else
    println("Precompile execution trace: $(abspath(PRECOMPILE_SCRIPT))")
end

kwargs = Dict{Symbol,Any}(
    :sysimage_path => SYSIMAGE_PATH
    # add other kwargs you already use...
    # Use the Earth example as a precompile execution trace so that the
    # methods hot-called during a real run are compiled into the image.
    # CI skips this trace by default to avoid depending on native runtime data.
    precompile_execution_file = PRECOMPILE_SCRIPT,
    project = PROJECT_ROOT,
)

if precompile_execution_file !== nothing
    kwargs[:precompile_execution_file] = precompile_execution_file  # String or Vector{String}
end

PackageCompiler.create_sysimage(packages; kwargs...)

println("Done.  Launch Julia with:")
println("  julia --project=. --sysimage SpaceAGORA.so <script.jl>")
