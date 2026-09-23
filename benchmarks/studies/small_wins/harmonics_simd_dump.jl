# Reproduces the "after" (current tree) half of the WS11f item 1 dump/cmp
# proof for `@simd ivdep` on the flat harmonics batch kernel
# (`_harmonics_flat_batch_kernel!`, src/dynamics/coupled/perturbations.jl).
#
# `@simd ivdep` was already applied and shipped at this branch's base commit
# (cd212833aa, "The batched harmonics loops vectorise at every batch size" --
# see `git log -1 cd212833aa`), so there is no source change for this item to
# make; the work was re-verifying bit-identity holds at this tip and fixing a
# stale doc comment that still said the kernel carried no `@simd` at all.
#
# The dump/cmp tool is `benchmarks/studies/third_body_cost/variants.jl --dump`
# (per the WS11 contract, read-only). This script just wraps the two
# invocations for the contract's reference cases -- the 256-spacecraft L50
# vacuum constellation, and the 4096-spacecraft P2 case
# (`gravity_4096sat_l50_vacuum_5800s` in
# benchmarks/studies/parallelization_performance/cases.jl) -- so they can be
# re-run identically later.
#
# Usage:
#   julia --project=. --threads=1 benchmarks/studies/small_wins/harmonics_simd_dump.jl <outdir>
#
# To reproduce the "before" (no-`@simd`) half: replace every
# `@inbounds @simd ivdep for b = 1:B` in `_harmonics_flat_batch_kernel!` with
# `@inbounds for b = 1:B` (three occurrences), re-run this script into a
# second `<outdir>`, restore the file, and `cmp` the two `*_vacuum.bin` pairs.
# That is exactly how the numbers in
# benchmarks/studies/small_wins/results/harmonics_simd_ratio.csv were made;
# see docs/architecture/small_wins_20260923.md for the measured bytes and
# cmp results.

outdir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "results", "harmonics_simd_dump")
mkpath(outdir)

const STUDY_DIR = normpath(joinpath(@__DIR__, "..", "third_body_cost"))
const VARIANTS = joinpath(STUDY_DIR, "variants.jl")

function run_dump(n::Int, mission::Float64, tag::String)
    dump_prefix = joinpath(outdir, tag)
    cmd = `julia --project=. --threads=1 $VARIANTS --n=$n --mission=$mission --variants=vacuum --repeats=1 --dump=$dump_prefix`
    println("running: ", cmd)
    run(cmd)
    return dump_prefix * "_vacuum.bin"
end

run_dump(256, 5800.0, "n256")
run_dump(4096, 5800.0, "n4096_p2")
println("dumps written under ", outdir)
