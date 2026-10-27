# StarterName

A companion package for [SpaceAGORA](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl)
that serves its external-propagation boundary: spacecraft that an outside integrator
owns for a whole run. Replace this paragraph with the engine it brings.

## The boundary

The companion reaches SpaceAGORA only through the versioned facade
`SpaceAGORA.ExternalPropagationServices`. It checks the facade's exact contract
version (`CORE_CONTRACT_VERSION`) before it defines any method, and again when it
loads, so an incompatible core is refused before any of the companion's methods
exist. `CircularOrbitPropagator` is an example owner; replace it with the
companion's own engine. It implements every function in the facade's
`IMPLEMENTED_FUNCTIONS`.

A companion for a different boundary needs that boundary's facade in SpaceAGORA
first, added and reviewed there; then point this starter at it.

## Set up and test

Use Julia 1.12, matching the pinned core. Run the setup commands from this
generated repository's root, beside `SpaceAGORASource.toml`. The sibling
directories `../StarterName-core` and `../StarterName-env` must be new; these
commands do not reuse or change another SpaceAGORA checkout:

```bash
git clone --filter=blob:none -c lfs.fetchexclude="*" https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git ../StarterName-core
git -C ../StarterName-core checkout 0000000000000000000000000000000000000000
julia --startup-file=no ../StarterName-core/scripts/contributor_environment.jl companion . ../StarterName-env --manifest core
```

`--filter=blob:none` skips the contents of old file versions; Git downloads a
commit's files when you check it out. `-c lfs.fetchexclude="*"` keeps this
clone, its checkouts and later pulls from downloading SpaceAGORA's large data
files (Git LFS), which a contributor repository does not use. Without Git LFS
installed the setting does nothing, and nothing is downloaded either.

Then run `julia --project=../StarterName-env -e 'using Pkg; Pkg.test("StarterName")'`.
That line is quoted for bash or zsh (on Windows, Git Bash). In PowerShell or
Command Prompt, start `julia --project=../StarterName-env` and enter
`using Pkg; Pkg.test("StarterName")` at the `julia>` prompt; the three commands
above run there as written. Setup writes the separate environment and can
download packages; it does not install the companion's native engine for you.

The checkout command above holds the same commit as `SpaceAGORASource.toml`. To
move to another SpaceAGORA commit, run `repin` from a checkout of that commit:

```bash
julia --startup-file=no ../StarterName-core/scripts/contributor_environment.jl repin companion .
```

`repin` updates the pin file, the workflow call and `[compat]`. It does not
change this README or `CORE_CONTRACT_VERSION`: put the new `rev` from
`SpaceAGORASource.toml` in the checkout command yourself, and update
`CORE_CONTRACT_VERSION` when the facade's contract has moved.

The tests check the facade's shape at the companion's contract version, a short
run with the example owner, and that a session that loads only SpaceAGORA loads
none of the companion. `test/incompatible_core.jl` (an incompatible contract is
refused) and `test/load_order.jl core-first` / `companion-first` run in fresh
processes; CI runs all of them.

## License

`LICENSE` is the MIT license, as for SpaceAGORA, with the year and the copyright
holder the generator was given (`--holder`). Replace `Your name or organization`
if it is still there, or replace the file to choose another license.

## Claims

Passing the conformance checks shows that the companion meets the boundary; the
engine's correctness and every claim made from it are the author's work.
