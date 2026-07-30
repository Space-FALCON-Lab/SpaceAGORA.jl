# Navigation examples

This folder contains the runners, analyses, demos, and plotting tools for the
spaceborne angles-only multi-target navigation example.


## Executables

The files intended to be executed directly are:

| File | Purpose |
|---|---|
| `examples/run_navigation.jl` | Run one simulation, a method comparison, or a Monte Carlo campaign |
| `examples/analyze_navigation.jl` | Analyze saved campaigns or generate plots and animations |
| `examples/navigation/demos/run_nominal_proposed.jl` | Short self-contained demo of the proposed method |
| `examples/navigation/demos/run_nominal_no_da.jl` | Short self-contained demo without data association |

`examples/navigation.jl`, the files under `examples/navigation/run/` and
`examples/navigation/analysis/`, and the Python files under
`scripts/plotting/navigation/` are implementation files called by these entry
points. They normally do not need to be executed directly.

## Quick start

Run the proposed-method demo:

```bash
julia --project=. examples/navigation/demos/run_nominal_proposed.jl
```

Run the no-data-association demo:

```bash
julia --project=. examples/navigation/demos/run_nominal_no_da.jl
```

Each demo creates a new temporary directory, so it cannot overwrite a previous
demo or a Monte Carlo campaign. At the end it prints commands that can be
copied to:

- open the numerical results;
- generate and open the static plots;
- generate and open the 3D animation.


## Running simulations

The main runner is:

```bash
julia --project=. examples/run_navigation.jl
```

With no command-line arguments, it reads `EDITOR_RUN_CONFIG` at the top of
`examples/run_navigation.jl`. This makes it possible to edit the configuration
and start the file with the VS Code Play button:

```julia
const EDITOR_RUN_CONFIG = (
    mode="single",
    case="proposed",
    output="",
    mission_time=600,
    targets=300,
    runs=1,
    resume=true,
    sweeps="false_alarm",
)
```

The fields have the following meaning:

| Field | Meaning |
|---|---|
| `mode` | Type of run: single, comparison, nominal Monte Carlo, or sensitivity analysis |
| `case` | Navigation method; multiple methods are provided as a comma-separated list |
| `output` | Output root. Leave empty to use the default location |
| `mission_time` | Simulated duration in seconds |
| `targets` | Number of generated target objects; the 16 observers are unchanged |
| `runs` | Number of realizations in either Monte Carlo campaign |
| `resume` | If `true`, reuse completed compatible realizations |
| `sweeps` | Fault families included in `sensitivity-monte-carlo` |



### Navigation cases

| Case | Meaning |
|---|---|
| `proposed` | Complete proposed distributed data-association and filtering method |
| `no-da` | Direct truth-based measurement routing, with the M2T and T2T association logic bypassed |
| `centralized_oracle` | Centralized reference estimator using known target identities |
| `independent_local_da` | Local data association and filtering without distributed track fusion |
| `distributed_oracle_da` | Distributed estimator with known target identities |
| `baseline_da` | Optional legacy baseline implementation |




### Nominal Monte Carlo campaign

Use `monte-carlo` for the nominal campaign:

```julia
mode="monte-carlo"
case="proposed,centralized_oracle,independent_local_da,distributed_oracle_da"
output=""
mission_time=10000
targets=300
runs=100
resume=true
```

With an empty `output`, results are saved under:

```text
output/navigation/nominal/
```


### Sensitivity Monte Carlo campaign

Use `sensitivity-monte-carlo` for the persistent measurement-fault campaign:

```julia
mode="sensitivity-monte-carlo"
case="proposed"
output=""
mission_time=10000
targets=300
runs=100
resume=true
sweeps="bias,misdetection,false_alarm"
```

With an empty `output`, results are saved under:

```text
output/navigation/sensitivity/
```

The available fault families are:

| Sweep | Model |
|---|---|
| `bias` | One independent angular-bias vector is generated for each observer and held constant over the realization, representing a fixed residual sensor misalignment |
| `misdetection` | Each available target detection is independently removed with the configured probability |
| `false_alarm` | The number of false detections at each observer epoch follows a Poisson distribution |

The Monte Carlo levels are:

- bias: `1e-5`, `5e-5`, and `1e-4` rad;
- misdetection probability: `0.05`, `0.10`, and `0.20`;
- false-alarm Poisson mean: `0.20`, `0.40`, and `0.80` false detections
  per observer epoch.

Every selected sensitivity campaign also includes the nominal configuration.
To run only one fault family, use for example:

```julia
sweeps="false_alarm"
```

## Analyzing saved results

The analysis entry point is:

```bash
julia --project=. examples/analyze_navigation.jl
```

For the normal workflow, edit only `EDITOR_ANALYSIS_MODE` at the top of the
file and press the VS Code Play button.

For the nominal campaign:

```julia
const EDITOR_ANALYSIS_MODE = "monte-carlo"
```

This reads:

```text
output/navigation/nominal/
```

For the measurement-fault campaign:

```julia
const EDITOR_ANALYSIS_MODE = "sensitivity"
```

This normally reads:

```text
output/navigation/sensitivity/
```

For backward compatibility, if that directory does not contain
`stress_run_status.csv`, the analyzer also checks the former
`output_navigation_stress_monte_carlo` location in the repository root and
under `~/SpaceAGORA_outputs/`. An explicit input always takes precedence:

```bash
julia --project=. examples/analyze_navigation.jl \
  --mode=sensitivity \
  --input=/path/to/sensitivity-campaign
```

These campaign analyses print the complete numerical summaries. Reported
intervals are the median and the 5th--95th percentile across realizations.
Sensitivity configurations retained in an older dataset but no longer present
in the current campaign definition are reported as ignored and are not included
in the generated summaries.

The saved pairwise IOD campaign is analyzed with:

```bash
julia --project=. examples/analyze_navigation.jl --mode=iod-pairwise
```

This reads `output/navigation/iod_pairwise/` by default.


## Output locations

The default persistent structure is:

```text
output/navigation/
├── nominal/
│   ├── run_<id>/<case>/
│   ├── mc_*.csv
│   └── iod_geometry_analysis/
├── sensitivity/
│   ├── run_<id>/<configuration>/<case>/
│   └── stress_*.csv
├── iod_pairwise/
├── timing/
└── figures/
```

Temporary runs are stored under the operating-system temporary directory:

```text
spaceagora_navigation/
├── demos/
└── runs/
    ├── single/
    ├── comparison/
    └── sensitivity_sweep/
```

Set `output` to a relative repository path or to an absolute path on Google
Drive or an external disk when a different storage location is required.
Status tables store paths relative to the campaign root, so an entire campaign
directory can be moved without breaking the analysis.

When `resume=true`, compatible realizations are reused and their existing status
entries are preserved. Optional IOD geometry and pairwise diagnostics are kept
in their dedicated outputs and are not required by the nominal or sensitivity
campaigns.

## Plotting dependencies

Install the Python plotting dependencies once:

```bash
python3 -m venv .venv-navigation
.venv-navigation/bin/python -m pip install \
  -r scripts/plotting/navigation/requirements.txt
```

`analyze_navigation.jl` detects `.venv-navigation` automatically; the
environment does not need to be activated for every plot command.

MP4 export additionally requires `ffmpeg`. If it is unavailable, the animation
script saves a GIF with the same base filename.
