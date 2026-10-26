# Your First Controller Package

Use the **`package` starter** to keep a reusable feedback controller in its own
Julia package. The starter's RPO planner is an example to replace: a controller
implements `AbstractControlEffectorModel` directly and needs no planning family,
service facade or change to SpaceAGORA.

This recipe packages the `RateDamper` from
[Feedback Controllers and Estimators](gnc_feedback.md#Worked-example-A:-a-rate-damper) as
`RateDamperDemo`. It keeps the controller, the simulation setup and the tests in
separate files. The example uses an ideal body torque, with no propellant use or
mounted hardware; it needs no optional companion, native GRAM or SPICE data.

## Generate the package

From a clean SpaceAGORA checkout named `SpaceAGORA.jl`, at the commit you intend
to use:

```sh
julia --startup-file=no --compiled-modules=existing scripts/new_contributor_repository.jl package ../RateDamperDemo --name RateDamperDemo
cd ../RateDamperDemo
```

The generator creates the files and a new package UUID. The directory must be
outside SpaceAGORA and new or empty. Keep the generated `SpaceAGORASource.toml`
and `.github/workflows/ci.yml`: their core revision and workflow revision agree.
Keep the generated package name, UUID, version and SpaceAGORA compatibility
entry in `Project.toml`.

Replace the planner example in these files:

| File | What it holds for this package |
|---|---|
| `src/RateDamperDemo.jl` | the `RateDamper` type and its public SpaceAGORA methods |
| `examples/configuration.jl` (new) | one short simulation setup shared by the example and tests |
| `examples/run_in_simulation.jl` | a run that checks the damping result |
| `test/runtests.jl` | public controller conformance on that simulation setup |
| `test/method_tests.jl` | checks of the damping law |
| `README.md` | what the controller does, its units, scope and setup commands |

Update `CITATION.cff` and the copyright holder in `LICENSE` for your package.
Keep the generated `.gitignore`. A package does not commit the separate test
environment or its Manifest.

## Declare the dependencies

The controller uses `SpaceAGORA` and `StaticArrays`. Add these entries to the
existing tables in `Project.toml`:

```toml
# In [deps]
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

# In [compat]
StaticArrays = "1"
```

Keep the existing `SpaceAGORA` dependency and compatibility entry. Remove
`LinearAlgebra` and `Random` from `[deps]` and `[compat]` after replacing the
planner code, since this recipe does not import them. Keep `Test` in `[extras]`
and the test target: it belongs to tests and the checked example, not the
controller's imports.

Do not add a `[sources]` entry for SpaceAGORA. The separate environment below
loads the core checkout at the recorded revision by path.

## Move the controller into its module

Replace `src/RateDamperDemo.jl`. Begin it with:

```julia
module RateDamperDemo

using SpaceAGORA
using StaticArrays

export RateDamper
```

Immediately after that header, copy the **complete controller code block** from
“Worked example A: a rate damper” in
[Feedback Controllers and Estimators](gnc_feedback.md#Worked-example-A:-a-rate-damper). In the matching
source, `examples/gnc_algorithm_env/feedback_controller.jl`, this is the region
between `# snippet: controller` and its `# end snippet`; omit the marker lines.
Copy the type and every method in that block, including the run-start reset,
fresh binding, inactive-command release, `checkpoint_state` and
`restore_checkpoint_state!`. Finish the file with:

```julia
end # module RateDamperDemo
```

Keep the simulation and `using Test` out of this module. Extending methods as
`SpaceAGORA.update_control!`, `SpaceAGORA.control_wrench` and the other public
names requires no import from an internal module.

The gain, torque limit, period and spacecraft scope are configuration: the
initializer and binding preserve them. `torque_body` and `tick_s` are owned run
state, not settings for an initial command. A fresh run resets those two fields;
set the gain or limit to configure the law rather than editing its held output.

The tick computes and holds the torque. The output hook reads that held value
and checks `x.spacecraft_index`. The two checkpoint methods save `torque_body`
and `tick_s` and restore them after initialization, so a resumed run holds the
same command until the next tick. Preserve this state handling when changing
the law.

## Share one nonzero-spin simulation setup

Create `examples/configuration.jl`. Assemble it from this header, the existing
example's configuration block, and the footer below:

```julia
using SpaceAGORA, StaticArrays, RateDamperDemo

function rate_damper_case(; damper=RateDamper(), duration_s=12.0,
                         initial_rate_rad_s=0.05)
```

Inside the function, copy the `run` code block from the rate-damper example,
starting at `root = Link(...)` and ending at the closing parenthesis of
`config = SimulationConfiguration(...)`. Do **not** copy
`solution = run_simulation(...)`. Make these three changes to the copied block:

1. Replace `ω0 = 0.05` with `ω0 = Float64(initial_rate_rad_s)`.
2. Remove `damper = RateDamper()`: use the function's `damper` argument.
3. Replace `mission_s = 120.0` with `mission_s = Float64(duration_s)`.

Keep the `ControlModel` containing that same `damper` object, with
`control_rates=[damper.period_s]`. Then append:

```julia
    return (; damper, configuration=config, initial_rate_rad_s=ω0,
             inertia_zz=root.inertia[3, 3], duration_s=mission_s)
end

# For this cube, initial spin and unsaturated single-axis damping only.
# The first tick is at one period, so the initial interval has zero torque.
expected_final_rate(case) = case.initial_rate_rad_s *
    (1 - case.damper.gain_n_m_s * case.damper.period_s / case.inertia_zz) ^
    (floor(Int, case.duration_s / case.damper.period_s) - 1)
```

This is the same public setup as the worked example, shortened for tests. Its
nonzero initial spin matters: a controller checked only while its command stays
zero can conceal a missing checkpoint state. The analytic expression here is
for the default controller and a duration that is a whole number of control
periods; changing the initial command, saturating the torque or changing the
physical case needs a corresponding method check.

## Test the interface and the method

Replace `test/runtests.jl` with:

```julia
using Test, SpaceAGORA, RateDamperDemo
using SpaceAGORA.Conformance: conformance_fixture, test_control_model

include(joinpath(@__DIR__, "..", "examples", "configuration.jl"))

@testset "RateDamperDemo" begin
    case = rate_damper_case()
    report = test_control_model(case.damper;
        fixture=conformance_fixture(configuration=case.configuration))
    checkpoint = only(filter(check -> check.id === :checkpoint_round_trip,
                             report.checks))
    @test checkpoint.status === :pass

    include("method_tests.jl")
end
```

The fixture must contain the exact controller passed to `test_control_model`,
not a second equal-looking instance. Run with at least two threads so the
concurrent-run check applies. Read any skipped or advisory checks in the report;
the explicit assertion above requires a successful checkpoint round trip even
when the pinned core reports that check as advisory.

Replace `test/method_tests.jl` with:

```julia
using Test, SpaceAGORA, StaticArrays, RateDamperDemo

@testset "the damping law" begin
    limited = RateDamper(gain_n_m_s=2.0, max_torque_n_m=0.1)
    state = (sc=[(ω=SVector(1.0, 0.0, 0.0),)],)
    update_control!(limited, state, nothing, 1.0, 1)
    @test held_command(limited, 1) == SVector(-0.1, 0.0, 0.0)
    @test held_command_time(limited, 1) == 1.0

    case = rate_damper_case()
    solution = run_simulation(case.configuration; return_solution=true)
    final_rate = solution.u[end].sc[1].ω[3]
    @test final_rate ≈ expected_final_rate(case) rtol=1e-6
    @test 0 < final_rate < case.initial_rate_rad_s
end
```

Conformance checks scheduling, state ownership, reset, binding and continuation.
The separate method test checks the torque limit and the expected physical
result. Add the scientific checks your own law needs; interface conformance does
not establish its accuracy or suitability for a mission.

Replace `examples/run_in_simulation.jl` with:

```julia
using SpaceAGORA, RateDamperDemo, Test
include(joinpath(@__DIR__, "configuration.jl"))

case = rate_damper_case()
solution = run_simulation(case.configuration; return_solution=true)
final_rate = solution.u[end].sc[1].ω[3]
@test final_rate ≈ expected_final_rate(case) rtol=1e-6
@test 0 < final_rate < case.initial_rate_rad_s
println("Body-rate damping agrees with the sampled-control analytic result.")
```

This preserves the example path already named in the generated CI workflow.

## Build the environment and run the checks

After editing the dependency declarations, from `RateDamperDemo/`:

```sh
julia --startup-file=no --compiled-modules=existing --project=@stdlib ../SpaceAGORA.jl/scripts/contributor_environment.jl check package .
julia --startup-file=no --compiled-modules=existing --project=@stdlib ../SpaceAGORA.jl/scripts/contributor_environment.jl scan package .
julia --startup-file=no --compiled-modules=existing --project=@stdlib ../SpaceAGORA.jl/scripts/contributor_environment.jl package . ../RateDamperDemo-env --manifest core
julia --startup-file=no --compiled-modules=existing --threads=2 --depwarn=error --project=../RateDamperDemo-env -e 'using Pkg; Pkg.test("RateDamperDemo"; julia_args=["--compiled-modules=existing", "--threads=2", "--depwarn=error"])'
julia --startup-file=no --compiled-modules=existing --threads=2 --depwarn=error --project=../RateDamperDemo-env examples/run_in_simulation.jl
```

The fourth line, `Pkg.test`, is quoted for bash or zsh (on Windows, Git Bash). In
PowerShell or Command Prompt, start Julia with the same options,
`julia --startup-file=no --compiled-modules=existing --threads=2 --depwarn=error --project=../RateDamperDemo-env`,
and enter the code between the single quotes at the `julia>` prompt; the other
lines run there as written.

The environment directory must be new or empty and outside both repositories.
Setup may download dependencies missing from the Julia depot. Once prepared,
these examples need no external simulation data. The core checkout must still
be clean and match `SpaceAGORASource.toml`; tests run from the separate environment,
not by activating the package directory alone.

The generated workflow calls SpaceAGORA's reusable contributor CI. Keep
`kind: package` and `examples: examples/run_in_simulation.jl`. It checks the pin,
scans for internal names, builds the environment, precompiles, and runs the tests
and example with strict deprecation handling and two threads.

To adopt a newer core, use a clean sibling checkout of the chosen commit and,
from `RateDamperDemo/`, run:

```sh
julia --startup-file=no --compiled-modules=existing --project=@stdlib ../SpaceAGORA.jl/scripts/contributor_environment.jl repin package .
```

`repin` updates the recorded revision, workflow call and SpaceAGORA compatibility
entry together. Build a new separate environment, then rerun the checks before
committing the pin change. Keep the full commit pin; do not replace it with a
moving branch name.

## Before your first pull request

Open this pull request in your new package repository. Using the existing
feedback hooks does not require a pull request to SpaceAGORA.

- [ ] The package source imports public SpaceAGORA names only; the internal-name
      scan passes.
- [ ] The planner implementation, planner tests and planner claims have been
      replaced, including the README and configured example.
- [ ] Units, frames, torque limits, tick period and intended use are documented.
- [ ] Run state resets correctly, bound copies own their state, and inactive
      spacecraft receive the controller's off command.
- [ ] Conformance has run on a nonzero-spin case with two threads, and
      `checkpoint_round_trip` is `:pass`; skipped and advisory checks have been read.
- [ ] Method tests check the physical result as well as the interface, and both
      local test and example commands succeed.
- [ ] The core pin, workflow revision and compatibility entry agree; the tested
      environment is separate and no generated environment files are committed.
- [ ] The README has the setup commands, license and citation details are yours,
      and the PR states what the controller implements and which checks ran.
