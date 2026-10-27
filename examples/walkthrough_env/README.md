# Extension walkthrough scripts

The complete scripts of the two extension walkthroughs. CI runs them from this
separate project, which uses only the public `SpaceAGORA` API.

| Script | Page |
|---|---|
| `my_effector.jl` | [Adding a Force or Torque of Your Own](../../docs/src/user/custom_effector.md) |
| `my_actuator.jl` | [Adding an Actuator of Your Own](../../docs/src/user/custom_actuator.md) |

From the repository root, prepare the environment once:

```sh
julia --project=examples/walkthrough_env examples/walkthrough_env/setup.jl
```

Setup copies the tracked `Project.template.toml` to a local `Project.toml`,
starts `Manifest.toml` from the repository manifest so the recorded package
versions are kept, and adds this checkout by path. Git ignores both generated
files and the `output/` folder the scripts write their results to. On a fresh
Julia depot setup downloads the packages it lacks, so it needs network access
the first time. After that the scripts run offline, as in CI, with no optional
package, native GRAM or SPICE kernels. Run one script, or both as CI does:

```sh
julia --project=examples/walkthrough_env --depwarn=error examples/walkthrough_env/my_effector.jl
julia --project=examples/walkthrough_env --depwarn=error examples/walkthrough_env/smoke.jl
```

Each page shows its whole script, the region between `# snippet: script` and
its `# end snippet`, and quotes shorter excerpts from the regions nested inside
it, all without the marker lines. `test/gates/ci_hpc_extensibility_docs_gate.jl`
compares the pages with the scripts, so change both together. The checks after
the script region test what the page says about the result.
