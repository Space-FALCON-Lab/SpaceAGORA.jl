# StarterName

An RPO reference planner for [SpaceAGORA](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl).
Replace this paragraph with what the planner does and when to use it.

## The planner

`StraightLinePlanner` is an example; replace it with the package's own method.
It subtypes `AbstractRPOPlanner`, so it belongs to SpaceAGORA's RPO reference
family: the family builds each request, checks every proposed reference with
its own trusted code and installs only a copy that passed. The planner declares
what it supports (`planner_capabilities`) and whether it draws from its random
stream (`consumes_rng`), and proposes references (`plan_rpo!`). A package cannot
add a family of its own; a new family goes into SpaceAGORA through a pull
request.

## Set up and test

The package is developed against the SpaceAGORA commit in `SpaceAGORASource.toml`,
in a separate environment that develops that commit by path. Use Julia 1.12,
matching the pinned core. Run the setup commands from this generated repository's
root, beside `SpaceAGORASource.toml`. The sibling directories
`../StarterName-core` and `../StarterName-env` must be new; these commands do not
reuse or change another SpaceAGORA checkout:

```bash
git clone --filter=blob:none -c lfs.fetchexclude="*" https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git ../StarterName-core
git -C ../StarterName-core checkout 0000000000000000000000000000000000000000
julia --startup-file=no ../StarterName-core/scripts/contributor_environment.jl package . ../StarterName-env --manifest core
```

`--filter=blob:none` skips the contents of old file versions; Git downloads a
commit's files when you check it out. `-c lfs.fetchexclude="*"` keeps this
clone, its checkouts and later pulls from downloading SpaceAGORA's large data
files (Git LFS), which a contributor repository does not use. Without Git LFS
installed the setting does nothing, and nothing is downloaded either.

The checkout command holds the same commit as `SpaceAGORASource.toml`. Setup
writes the separate environment and can download packages.

Then run `julia --project=../StarterName-env -e 'using Pkg; Pkg.test("StarterName")'`.
That line is quoted for bash or zsh (on Windows, Git Bash). In PowerShell or
Command Prompt, start `julia --project=../StarterName-env` and enter
`using Pkg; Pkg.test("StarterName")` at the `julia>` prompt; the three commands
above run there as written.

`test/runtests.jl` runs SpaceAGORA's planner conformance checks on the family's
representative requests and a short simulation; `test/method_tests.jl` holds the
package's own tests of its method. `examples/run_in_simulation.jl` runs the
planner in the RPO pilot; CI runs it after the tests.

To move to another SpaceAGORA commit, check that commit out in
`../StarterName-core` and run its `repin` from this repository's root:

```bash
julia --startup-file=no ../StarterName-core/scripts/contributor_environment.jl repin package .
```

`repin` updates the pin file, both commit references in the workflow call and
`[compat]`, but not this README: put the new `rev` in the checkout command
yourself, then rebuild the environment in a new directory.

People who use the package install SpaceAGORA themselves first: the pin file
only serves this package's own development and CI. `[compat]` admits one
SpaceAGORA release series.

## License

`LICENSE` is the MIT license, as for SpaceAGORA, with the year and the copyright
holder the generator was given (`--holder`). Replace `Your name or organization`
if it is still there, or replace the file to choose another license.

## Claims

Passing the conformance checks shows that the planner meets the family's
interface; the method's correctness and every claim made from it are the
author's work.
