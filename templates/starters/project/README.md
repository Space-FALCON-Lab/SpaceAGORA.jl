# StarterName

A research study built on [SpaceAGORA](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl).
Replace this paragraph with the study's question.

## Layout

- `src/StarterName.jl`: the study's own models (here `LowThrustPush`, an example
  force model) and its case (`case_configuration`, `run_case`).
- `scripts/run_case.jl`: the worked example; CI runs it after the tests.
- `test/runtests.jl`: conformance checks of the study's models and a short run.
- `SpaceAGORASource.toml`: the SpaceAGORA commit the study is tested against.
- `Manifest.toml`: committed, so every result has its exact package versions.
  The generator starts it as a copy of the pinned commit's manifest.

## Set up and test

The study is tested in a separate environment that develops the pinned SpaceAGORA
commit by path. Use Julia 1.12, matching the pinned core. Run the setup commands
from this generated repository's root, beside `SpaceAGORASource.toml`.
The sibling directories `../StarterName-core` and `../StarterName-env` must be
new; these commands do not reuse or change another SpaceAGORA checkout:

```bash
git clone --filter=blob:none -c lfs.fetchexclude="*" https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git ../StarterName-core
git -C ../StarterName-core checkout 0000000000000000000000000000000000000000
julia --startup-file=no ../StarterName-core/scripts/contributor_environment.jl project . ../StarterName-env --manifest committed
```

`--filter=blob:none` skips the contents of old file versions; Git downloads a
commit's files when you check it out. `-c lfs.fetchexclude="*"` keeps this
clone, its checkouts and later pulls from downloading SpaceAGORA's large data
files (Git LFS), which a contributor repository does not use. Without Git LFS
installed the setting does nothing, and nothing is downloaded either. These
commands, and the test command below, run as written in bash or zsh, PowerShell
and Command Prompt.

The checkout command holds the same commit as `SpaceAGORASource.toml`. Setup
writes the separate environment and can download packages.

Then, from `test/`:

```bash
julia --startup-file=no --project=../../StarterName-env runtests.jl
```

CI (`.github/workflows/ci.yml`) does the same through SpaceAGORA's reusable
workflow.

To move to another SpaceAGORA commit, check that commit out in
`../StarterName-core` and run its `repin` from this repository's root:

```bash
julia --startup-file=no ../StarterName-core/scripts/contributor_environment.jl repin project .
```

`repin` updates the pin file, both commit references in the workflow call and
`[compat]`, but not this README: put the new `rev` in the checkout command
yourself. Then rebuild the environment in a new directory and commit the
updated `Manifest.toml`.

## Results and claims

Keep results out of this repository (Git ignores `output/` and `results/`) and
record the SpaceAGORA commit with each one. The conformance tests check that the
study's models meet SpaceAGORA's interfaces; validating the models, the case and
every claim made from them is the author's work.
