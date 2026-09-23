# Resident memory of one process-pool worker running Earth native GRAM.
#
# Context: a 16-worker process pool running the P6p trace
# (aero_4096sat_l50_gram_process_100s) was OOM-killed on a 250 GB machine while
# its workers were starting, each having just printed "Reading MERRA2 data from
# .../MERRA2All_01.bin". This measures, in one process shaped like one such
# worker (Earth, native GRAM, one-spacecraft sample, 100 s, L50 + aero), what the
# resident set is at each stage, and in particular across the first GRAM
# native update, which is where the MERRA2 read happens.
#
# Every figure is read from /proc/self/status: VmRSS is the resident set at that
# moment, VmHWM the peak resident set of the process so far. Nothing is
# estimated.
#
# Usage (one process, capped so a runaway is killed alone):
#   systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q -- \
#     julia --project=. --threads=1 benchmarks/studies/gram_thread_scaling/worker_rss.jl

using Printf

function proc_status_mb(key::String)::Float64
    for line in eachline("/proc/self/status")
        startswith(line, key * ":") || continue
        return parse(Float64, split(line)[2]) / 1024.0   # kB -> MB
    end
    return NaN
end

const STAGES = Tuple{String, Float64, Float64}[]
function stage!(label::String)
    GC.gc(); GC.gc()
    rss = proc_status_mb("VmRSS"); hwm = proc_status_mb("VmHWM")
    push!(STAGES, (label, rss, hwm))
    @printf("%-58s VmRSS %8.1f MB   VmHWM %8.1f MB\n", label, rss, hwm)
    flush(stdout)
end

stage!("0 julia started")

const GTS_DIR = @__DIR__
const GTS_REPO_ROOT = normpath(joinpath(GTS_DIR, "..", "..", ".."))
const GTS_PPC_DIR = joinpath(GTS_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")
include(joinpath(GTS_PPC_DIR, "cli.jl"))
include(joinpath(GTS_PPC_DIR, "modes.jl"))
include(joinpath(GTS_PPC_DIR, "cases.jl"))
stage!("1 SpaceAGORA + harness loaded")

ppc_ensure_gramsuite_loaded!()
stage!("2 GRAMSuite loaded")

const EM = SpaceAGORA.SimulationModel.EnvironmentModels
const MODEL = ppc_gram_atmosphere_model("earth")
stage!("3 Earth GRAMAtmosphereModel constructed")

function one_sample_config()
    planet = Earth("", PPC_SPICE_PATH)
    return ppc_build_config(
        planet=planet,
        spacecraft=[ppc_spacecraft(planet; id=1)],
        mission_time_s=100.0,
        orientation_sim=false,
        dynamic_effectors=(ppc_harmonics_model(planet, 50), AerodynamicCoefficientfM()),
        density_model=MODEL,
        dt_max_orbit=5.0
    )
end
run_sample() = withenv("SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "1",
                       "SPACEAGORA_VACUUM_GRAM_CACHE" => "0") do
    SimulationEngine.run_simulation(one_sample_config(); isolate_state=false)
end

# The first sample is where a fresh worker reads MERRA2: the "Reading MERRA2
# data from ..." line is printed by the first native update, inside the solve
# (a bare GRAM call before any solve cannot run -- the engine is what loads the
# SPICE kernels GRAM's injected ephemeris needs). Every later sample prints the
# line again -- each solve reaches a freshly initialized native atmosphere -- so
# the difference between sample 1 and sample 2 is the first-read cost, not the
# whole read.
println("--- sample 1 (first native GRAM update, MERRA2 read) ---"); flush(stdout)
run_sample()
stage!("4 after sample 1 (1 sat, 100 s; includes MERRA2 read)")
run_sample()
stage!("5 after sample 2 (same model; MERRA2 read again)")

# A second native instance in the same process, the way the isolated pool and
# per-satellite instances create them: does each instance carry its own MERRA2?
second = deepcopy(MODEL)
stage!("6 second instance deepcopied")
println("--- first evaluation of the second instance follows ---"); flush(stdout)
EM._gram_core_density_state(second.core, 300.0e3, 0.3, 0.5, 0.0, true, ReentrantLock(), 200.0)
stage!("7 after second instance's first evaluation")

# A process-pool worker in a 16-wide pool over 4096 samples serves 256 of them
# in one process. Each sample above printed its own "Reading MERRA2" line, so
# whether resident memory accumulates per sample is the question that scales
# with the sample count rather than with the worker count.
const GTS_EXTRA = parse(Int, get(ENV, "GTS_RSS_EXTRA_SAMPLES", "30"))
for k in 1:GTS_EXTRA
    run_sample()
    (k % 10 == 0 || k == GTS_EXTRA) && stage!("8 after $(k) further samples")
end

println()
@printf("sample 1 incl. MERRA2 read (4 - 3):         %8.1f MB resident\n", STAGES[5][2] - STAGES[4][2])
@printf("sample 2, repeat read (5 - 4):              %8.1f MB resident\n", STAGES[6][2] - STAGES[5][2])
@printf("second instance incl. its read (7 - 5):     %8.1f MB resident\n", STAGES[8][2] - STAGES[6][2])
@printf("growth over the further samples (last - 7):  %8.1f MB resident\n", STAGES[end][2] - STAGES[8][2])
@printf("peak resident set of the process (VmHWM):   %8.1f MB\n", STAGES[end][3])
