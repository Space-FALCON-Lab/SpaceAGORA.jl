# StarterName

A research study built on [SpaceAGORA](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl).
Replace this paragraph with the study's question.

## Layout

- `src/StarterName.jl`: the study's own models (here `LowThrustPush`, an example
  force model) and its case (`case_configuration`, `run_case`).
- `scripts/run_case.jl`: the working example; CI runs it after the tests.
- `test/runtests.jl`: interface conformance of the study's models and a short run.
- `SpaceAGORASource.toml`: the SpaceAGORA commit the study is tested against.
- `Manifest.toml`: committed, so every result has its exact package versions.

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

`--filter=blob:none` skips the contents of old file versions: Git downloads the
files of a commit when you check it out. `-c lfs.fetchexclude="*"` keeps this
clone, its checkouts and later pulls from downloading SpaceAGORA's large data
files (Git LFS), which a contributor repository does not use. Without Git LFS
installed the setting does nothing, and nothing is downloaded either. These
commands, and the test command below, run as written in bash or zsh, PowerShell
and Command Prompt.

The generated checkout command contains the same exact revision as
`SpaceAGORASource.toml`. SpaceAGORA's `repin` does not rewrite this README:
after a repin, use the new `rev` from `SpaceAGORASource.toml` in the checkout
command. Setup writes the separate environment and can access package servers.

Then, from `test/`:

```bash
julia --startup-file=no --project=../../StarterName-env runtests.jl
```

CI (`.github/workflows/ci.yml`) does the same with SpaceAGORA's reusable workflow.
To move to a newer SpaceAGORA commit, change the commit in `SpaceAGORASource.toml`
and in the workflow call together, rebuild the environment and commit its Manifest.

## Results and claims

Keep results out of this repository (`output/` and `results/` are ignored) and
record the SpaceAGORA commit with each one. The conformance tests check that the
study's models meet SpaceAGORA's interface contracts; validating the models, the
case and every claim made from them is the author's work.
