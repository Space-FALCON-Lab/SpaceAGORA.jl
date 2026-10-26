# Extension walkthrough scripts

The complete scripts of the two extension walkthroughs, run by CI from this
separate project with the public `SpaceAGORA` surface only:

| Script | Page |
|---|---|
| `my_effector.jl` | [Adding a Force or Torque of Your Own](../../docs/src/user/custom_effector.md) |
| `my_actuator.jl` | [Adding an Actuator of Your Own](../../docs/src/user/custom_actuator.md) |

From the repository root, prepare the environment once:

```sh
julia --project=examples/walkthrough_env examples/walkthrough_env/setup.jl
```

Setup creates a local `Project.toml` from the tracked `Project.template.toml`,
seeds `Manifest.toml` from the repository manifest and adds this checkout by
path, keeping the recorded dependency versions. Both generated files, and the
`output/` folder the scripts write their results to, are ignored. Setup
downloads any of those versions your Julia depot does not hold yet, so on a
fresh depot it needs network access. After setup the scripts run offline, as CI
runs them, and need no optional package, native GRAM or SPICE kernels. Then run
one, or both as CI does:

```sh
julia --project=examples/walkthrough_env --depwarn=error examples/walkthrough_env/my_effector.jl
julia --project=examples/walkthrough_env --depwarn=error examples/walkthrough_env/smoke.jl
```

Each page shows its script between the `# snippet: script` and `# end snippet`
lines, and its shorter excerpts from the nested regions, without the marker
lines; `test/gates/ci_hpc_extensibility_docs_gate.jl` compares them, so change
the page and the script together. The checks after the region test what the
page states about the result.
