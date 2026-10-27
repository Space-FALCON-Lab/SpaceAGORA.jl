# Your First Controller Package

The `package` starter gives a reusable feedback controller a Julia package of
its own. The starter's example is an RPO planner, which you replace: a
controller subtypes `AbstractControlEffectorModel` directly and needs no
planning family, service facade or change to SpaceAGORA.

This recipe packages the `RateDamper` from
[Feedback Controllers and Estimators](gnc_feedback.md#Worked-example-A:-a-rate-damper)
as `RateDamperDemo`, with the controller, the simulation setup and the tests in
separate files. The controller applies an ideal body torque, with no propellant
use or mounted hardware, and the package needs no optional companion, native
GRAM or SPICE data.

## Generate the package

Start in a clean SpaceAGORA checkout named `SpaceAGORA.jl`, at the commit you
intend to use, and run:

```sh
julia --startup-file=no --compiled-modules=existing scripts/new_contributor_repository.jl package ../RateDamperDemo --name RateDamperDemo
cd ../RateDamperDemo
```

The generator writes the starter's files with a new package UUID and pins the
checkout's commit. The target directory must be outside SpaceAGORA and either
new or empty, and the generator refuses a checkout with uncommitted changes to
tracked files. Keep the generated `SpaceAGORASource.toml` and
`.github/workflows/ci.yml`, which name the same core commit, and keep the
generated package name, UUID, version and SpaceAGORA compatibility entry in
`Project.toml`.

Replace the planner example in these files:

| File | What it holds for this package |
|---|---|
| `src/RateDamperDemo.jl` | the `RateDamper` type and its public SpaceAGORA methods |
| `examples/configuration.jl` (new) | one short simulation setup shared by the example and tests |
| `examples/run_in_simulation.jl` | a run that checks the damping result |
| `test/runtests.jl` | public controller conformance on that simulation setup |
| `test/method_tests.jl` | checks of the damping law |
| `README.md` | what the controller does, its units, scope and setup commands |

Put your own details in `CITATION.cff` and the copyright holder in `LICENSE`.
The generator writes the holder into `LICENSE` for you when you pass it
`--holder "<name>"`. Keep the generated `.gitignore`: a package commits no
Manifest, and the test environment you build below lives outside the
repository.

## Declare the dependencies

The controller uses `SpaceAGORA` and `StaticArrays`. Add these entries to the
existing tables in `Project.toml`:

```toml
# In [deps]
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

# In [compat]
StaticArrays = "1"
```

Keep the existing `SpaceAGORA` dependency and compatibility entry. Once the
planner code is gone, remove `LinearAlgebra` and `Random` from `[deps]` and
`[compat]`, since this recipe does not import them. Keep `Test` in `[extras]`
and in the test target: the tests and the example use it, the controller does
not.

Do not add a `[sources]` entry for SpaceAGORA; the environment script refuses
one. The separate environment built below loads the core checkout at the
recorded commit by path.

## Move the controller into its module

Replace `src/RateDamperDemo.jl`, starting with this header:

```julia
module RateDamperDemo

using SpaceAGORA
using StaticArrays

export RateDamper
```

Right after the header, paste the whole controller code block from
[Worked example A: a rate damper](gnc_feedback.md#Worked-example-A:-a-rate-damper).
In its source, `examples/gnc_algorithm_env/feedback_controller.jl`, that block
is the region between `# snippet: controller` and its `# end snippet`; leave out
the marker lines. Copy the type and every method in the block, including the
reset at the run's start (`initialize_control!`), `bind_spacecraft`,
`release_command!`, `checkpoint_state` and `restore_checkpoint_state!`. End the
file with:

```julia
end # module RateDamperDemo
```

Keep the simulation and `using Test` out of this module. The methods extend
public names with qualified definitions, such as `SpaceAGORA.update_control!`
and `SpaceAGORA.control_wrench`, so nothing has to be imported from an internal
module.

The gain, torque limit, period and spacecraft index are configuration. The
initializer leaves them alone, and `bind_spacecraft` copies the gain, limit and
period into the copy it makes for a spacecraft. `torque_body` and `tick_s` are
run state the model owns, not settings for an initial command. Every run starts by
resetting them, so configure the law through the gain or the limit rather than
by editing its held output.

The tick computes and holds the torque. The output hook returns that held value
and checks `x.spacecraft_index`. The two checkpoint methods save `torque_body`
and `tick_s` and restore them after initialization, so a resumed run holds the
same command until the next tick. Keep this state handling when you change the
law.

## Share one nonzero-spin simulation setup

Create `examples/configuration.jl` from three parts: the header below, the
configuration code of the worked example, and the footer that follows.

```julia
using SpaceAGORA, StaticArrays, RateDamperDemo

function rate_damper_case(; damper=RateDamper(), duration_s=12.0,
                         initial_rate_rad_s=0.05)
```

Inside the function, paste the worked example's code block that builds the run,
from `root = Link(...)` through the closing parenthesis of
`config = SimulationConfiguration(...)`. Leave out
`solution = run_simulation(...)`. Then make three changes to the pasted code:

1. Replace `ω0 = 0.05` with `ω0 = Float64(initial_rate_rad_s)`.
2. Delete `damper = RateDamper()`, so the code uses the function's `damper`
   argument.
3. Replace `mission_s = 120.0` with `mission_s = Float64(duration_s)`.

The `ControlModel` keeps that same `damper` object, with
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

This is the worked example's public setup, shortened for tests. The nonzero
initial spin matters: a controller checked only while its command stays zero
can hide a missing checkpoint state. The closed-form expression holds for the
default controller and a duration of a whole number of control periods. If you
change the initial command, let the torque saturate or change the physical
case, write a method check to match.

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

The fixture must hold the very controller object passed to
`test_control_model`, not a second instance that looks the same, because the
check finds the model by identity. Run with at least two threads so that the
concurrent-run check applies. Read any skipped or advisory checks in the
report. The explicit assertion above requires the checkpoint round trip to pass
even if the pinned core reports that check as advisory.

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

Conformance checks the interface: the hooks; repeated runs, runs after a
`deepcopy` and concurrent runs; the reset at a run's start; pure output hooks;
the copies `bind_spacecraft` makes; the checkpoint round trip; and the outputs
for spacecraft outside the model's scope. The method test checks the torque limit
and the expected physical result. Add the scientific checks your own law needs;
interface conformance says nothing about its accuracy or its fitness for a
mission.

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

The generated CI workflow already names this example path.

## Build the environment and run the checks

After editing the dependency declarations, run these from `RateDamperDemo/`:

```sh
julia --startup-file=no --compiled-modules=existing --project=@stdlib ../SpaceAGORA.jl/scripts/contributor_environment.jl check package .
julia --startup-file=no --compiled-modules=existing --project=@stdlib ../SpaceAGORA.jl/scripts/contributor_environment.jl scan package .
julia --startup-file=no --compiled-modules=existing --project=@stdlib ../SpaceAGORA.jl/scripts/contributor_environment.jl package . ../RateDamperDemo-env --manifest core
julia --startup-file=no --compiled-modules=existing --threads=2 --depwarn=error --project=../RateDamperDemo-env -e 'using Pkg; Pkg.test("RateDamperDemo"; julia_args=["--compiled-modules=existing", "--threads=2", "--depwarn=error"])'
julia --startup-file=no --compiled-modules=existing --threads=2 --depwarn=error --project=../RateDamperDemo-env examples/run_in_simulation.jl
```

The first line checks the repository and its core pin, the second scans the
package for internal names, and the third builds the separate environment. The
fourth line, `Pkg.test`, is quoted for bash or zsh (on Windows, Git Bash). In
PowerShell or Command Prompt, start Julia with the same options,
`julia --startup-file=no --compiled-modules=existing --threads=2 --depwarn=error --project=../RateDamperDemo-env`,
and enter the code between the single quotes at the `julia>` prompt; the other
four lines run there as written.

The environment directory must be new or empty and outside both repositories.
Building it may download dependencies missing from the Julia depot; after that,
these examples need no external simulation data. The core checkout must still
be clean and at the commit in `SpaceAGORASource.toml`. The tests run in the
separate environment, not by activating the package directory alone.

The generated workflow calls SpaceAGORA's reusable contributor CI. Keep
`kind: package` and `examples: examples/run_in_simulation.jl`. The workflow
checks the pin, scans for internal names, builds the environment, precompiles,
and runs the tests and the example with `--depwarn=error` and two threads.

To move to a newer core, use a clean sibling checkout of the chosen commit and
run, from `RateDamperDemo/`:

```sh
julia --startup-file=no --compiled-modules=existing --project=@stdlib ../SpaceAGORA.jl/scripts/contributor_environment.jl repin package .
```

`repin` updates the commit and version in `SpaceAGORASource.toml`, the commit in
the workflow call and, when it does not admit the new core, the SpaceAGORA
compatibility entry. Build a new separate environment and rerun the checks
before you commit the change. Keep the full commit in the pin; a branch name is
not a pin.

## Before your first pull request

Open this pull request in your new package's repository. Using the existing
feedback hooks needs no pull request to SpaceAGORA.

- [ ] The package source imports public SpaceAGORA names only; the internal-name
      scan passes.
- [ ] The planner code, planner tests and planner claims are gone, from the
      README and the example too.
- [ ] Units, frames, torque limits, tick period and intended use are documented.
- [ ] Run state resets at a run's start, bound copies own their state, and
      inactive spacecraft receive the controller's off command.
- [ ] Conformance has run on a nonzero-spin case with two threads, and
      `checkpoint_round_trip` is `:pass`; skipped and advisory checks have been
      read.
- [ ] Method tests check the physical result as well as the interface, and the
      local test and example commands both succeed.
- [ ] The core pin, the workflow's commit and the compatibility entry agree; the
      tested environment is separate, and no generated environment files are
      committed.
- [ ] The README has the setup commands, the license and citation details are
      yours, and the pull request states what the controller implements and
      which checks ran.
