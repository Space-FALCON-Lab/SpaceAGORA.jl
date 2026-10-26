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
# The default precompile trace runs examples/AGORA_Basic_GRAMEarth.jl, which
# continues in the examples environment: prepare it first with
# `julia --startup-file=no examples/setup.jl`, or name another trace script in
# SPACEAGORA_SYSIMAGE_PRECOMPILE_FILE.
#
# The resulting sysimage is written to SpaceAGORA.so in the repository root.
# That file is listed in .gitignore and should NOT be committed to version
# control.  Rebuild it whenever the package set in Project.toml
# changes.

using PackageCompiler

# All packages from the root Project.toml that should be baked into the
# full local sysimage. PackageCompiler needs them listed by symbol, and each
# name must be :SpaceAGORA or a root Project.toml [deps] entry
# (test/gates/ci_sysimage_package_lists_gate.jl checks both lists).
# :SpaceAGORA itself is included so the package's own ~40+ src/ files get
# compiled into the image too -- without it, PackageCompiler only bakes in
# the dependencies, and `using SpaceAGORA` still recompiles the package's own
# code from scratch on every fresh process (measured: ~45-50s) even with the
# sysimage loaded.
const FULL_PACKAGES = [
    :SpaceAGORA,
    :Arrow,
    :AstroTime,
    :CSV,
    :ComponentArrays,
    :DataFrames,
    :DiffEqBase,
    :DiffEqCallbacks,
    :Interpolations,
    :LoopVectorization,
    :OrdinaryDiffEq,
    :Polyester,
    :Reexport,
    :SPICE,
    :SatelliteToolboxGeomagneticField,
    :SpecialFunctions,
    :StaticArrays,
]

# CI jobs benefit more from fast, reliable sysimage builds than from squeezing
# every optional plotting package into the image. Keep the runtime-heavy core
# stack and omit the largest optional frontends.
const CI_PACKAGES = [
    :SpaceAGORA,
    :Arrow,
    :AstroTime,
    :CSV,
    :ComponentArrays,
    :DataFrames,
    :DiffEqBase,
    :DiffEqCallbacks,
    :Interpolations,
    :OrdinaryDiffEq,
    :Polyester,
    :Reexport,
    :SPICE,
    :SatelliteToolboxGeomagneticField,
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
    # src/examples/Earth.jl never existed in this layout; examples live under
    # examples/ at the repo root. AGORA_Basic_GRAMEarth.jl is the smallest
    # real, GRAM-backed path (SPICE + GRAM atmosphere + a short solve), so it
    # exercises the packages most worth specializing ahead of time.
    default_script = joinpath(PROJECT_ROOT, "examples", "AGORA_Basic_GRAMEarth.jl")

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

# The default trace example stops without the examples environment; say so before
# the build starts rather than after it.
PRECOMPILE_SCRIPT !== nothing && endswith(PRECOMPILE_SCRIPT, joinpath("examples", "AGORA_Basic_GRAMEarth.jl")) &&
    !isfile(joinpath(PROJECT_ROOT, "examples", "Manifest.toml")) &&
    error("The default precompile trace runs an example: prepare the examples environment first with " *
          "`julia --startup-file=no examples/setup.jl`, or name another trace script in SPACEAGORA_SYSIMAGE_PRECOMPILE_FILE.")

println("Building sysimage → $(abspath(SYSIMAGE_PATH))")
println("This will take several minutes the first time …")
println("Sysimage package profile: $(CI_MIN_SYSIMAGE ? "ci-min" : "full")")
println("Packages baked in: $(length(PACKAGES))")
if PRECOMPILE_SCRIPT === nothing
    println("Precompile execution trace: skipped")
else
    println("Precompile execution trace: $(abspath(PRECOMPILE_SCRIPT))")
end

# --cpu-target=native (PackageCompiler's default) crashes on this machine
# with "LLVM ERROR: Cannot select: ... i64 = vscale" from HostCPUFeatures.jl's
# aarch64 SVE-detection code (Polyester/VectorizationBase, both real
# dependencies -- SpaceAGORA's own setup.jl uses Polyester directly, so this
# isn't avoidable by trimming packages). Reproduced identically at both
# ~15 MiB and ~141 MiB free memory, so it isn't OOM -- it's an AOT-codegen
# incompatibility with scalable-vector intrinsics on this LLVM/Apple Silicon
# combination. `cpu_target="generic"` sidesteps the SVE-feature-detection
# codegen path entirely; the sysimage's baked-in code is less
# CPU-tuned as a result, but Julia still JIT-specializes hot loops normally
# at runtime, so this only costs a bit of the "already precompiled" benefit
# for SIMD-heavy kernels, not correctness.
const CPU_TARGET = get(ENV, "SPACEAGORA_SYSIMAGE_CPU_TARGET", "generic")

kwargs = Dict{Symbol,Any}(
    :sysimage_path => SYSIMAGE_PATH,
    :project => PROJECT_ROOT,
    :cpu_target => CPU_TARGET,
)

if PRECOMPILE_SCRIPT !== nothing
    kwargs[:precompile_execution_file] = PRECOMPILE_SCRIPT
end

PackageCompiler.create_sysimage(PACKAGES; kwargs...)

println("Done.  Launch Julia with:")
println("  julia --project=. --sysimage SpaceAGORA.so <script.jl>")
