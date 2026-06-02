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
# full local sysimage. PackageCompiler needs them listed by symbol.
const FULL_PACKAGES = [
    :Arrow,
    :AssociatedLegendrePolynomials,
    :AstroTime,
    :CSV,
    :ComponentArrays,
    :ControlSystemsBase,
    :DataFrames,
    :DiffEqBase,
    :DiffEqCallbacks,
    :DifferentialEquations,
    :Interpolations,
    :LaTeXStrings,
    :LoopVectorization,
    :MatrixEquations,
    :NLsolve,
    :OrdinaryDiffEq,
    :PlotlyJS,
    :Plots,
    :Polyester,
    :PreallocationTools,
    :Quaternions,
    :Reexport,
    :Roots,
    :Rotations,
    :SPICE,
    :SatelliteToolbox,
    :SatelliteToolboxAtmosphericModels,
    :SatelliteToolboxGeomagneticField,
    :SatelliteToolboxGravityModels,
    :SatelliteToolboxTransformations,
    :SpecialFunctions,
    :StaticArrays,
]

# CI jobs benefit more from fast, reliable sysimage builds than from squeezing
# every optional plotting package into the image. Keep the runtime-heavy core
# stack and omit the largest optional frontends.
const CI_PACKAGES = [
    :Arrow,
    :AssociatedLegendrePolynomials,
    :AstroTime,
    :CSV,
    :ComponentArrays,
    :ControlSystemsBase,
    :DataFrames,
    :DiffEqBase,
    :DiffEqCallbacks,
    :DifferentialEquations,
    :Interpolations,
    :LaTeXStrings,
    :MatrixEquations,
    :NLsolve,
    :OrdinaryDiffEq,
    :Polyester,
    :PreallocationTools,
    :Quaternions,
    :Reexport,
    :Roots,
    :Rotations,
    :SPICE,
    :SatelliteToolbox,
    :SatelliteToolboxAtmosphericModels,
    :SatelliteToolboxGeomagneticField,
    :SatelliteToolboxGravityModels,
    :SatelliteToolboxTransformations,
    :SpecialFunctions,
    :StaticArrays,
]

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CI_MIN_SYSIMAGE = get(ENV, "SPACEAGORA_CI_MIN_SYSIMAGE", "0") == "1"
const PACKAGES = CI_MIN_SYSIMAGE ? CI_PACKAGES : FULL_PACKAGES

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
println("Sysimage package profile: $(CI_MIN_SYSIMAGE ? "ci-min" : "full")")
println("Packages baked in: $(length(PACKAGES))")
if PRECOMPILE_SCRIPT === nothing
    println("Precompile execution trace: skipped")
else
    println("Precompile execution trace: $(abspath(PRECOMPILE_SCRIPT))")
end

kwargs = Dict{Symbol,Any}(
    :sysimage_path => SYSIMAGE_PATH,
    :project => PROJECT_ROOT,
)

if PRECOMPILE_SCRIPT !== nothing
    kwargs[:precompile_execution_file] = PRECOMPILE_SCRIPT
end

PackageCompiler.create_sysimage(PACKAGES; kwargs...)

println("Done.  Launch Julia with:")
println("  julia --project=. --sysimage SpaceAGORA.so <script.jl>")
