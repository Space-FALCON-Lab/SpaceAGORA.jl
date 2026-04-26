# build_sysimage.jl
#
# Builds a custom Julia system image that has all SpaceAGORA.jl packages
# pre-compiled into it.  Subsequent Julia sessions launched with
#
#   julia --project=.SpaceAGORA --sysimage SpaceAGORA.so <script.jl>
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
# control.  Rebuild it whenever the package set in .SpaceAGORA/Project.toml
# changes.

using PackageCompiler

# All packages from .SpaceAGORA/Project.toml that should be baked into the
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

# The sysimage is written next to this build script's parent (i.e. repo root).
const SYSIMAGE_PATH = joinpath(@__DIR__, "..", "SpaceAGORA.so")

# Optional: provide a precompile script so PackageCompiler can record which
# methods are actually called and include their native code.  Point this at a
# lightweight example that exercises the core code paths.
const PRECOMPILE_SCRIPT = joinpath(@__DIR__, "..", "src", "examples", "Earth.jl")

println("Building sysimage → $(abspath(SYSIMAGE_PATH))")
println("This will take several minutes the first time …")

create_sysimage(
    PACKAGES;
    sysimage_path = SYSIMAGE_PATH,
    # Use the Earth example as a precompile execution trace so that the
    # methods hot-called during a real run are compiled into the image.
    # Remove or replace this line if GRAM / SPICE data are not available.
    precompile_execution_file = isfile(PRECOMPILE_SCRIPT) ? PRECOMPILE_SCRIPT : nothing,
    project = joinpath(@__DIR__, "..", ".SpaceAGORA"),
)

println("Done.  Launch Julia with:")
println("  julia --project=.SpaceAGORA --sysimage SpaceAGORA.so <script.jl>")
