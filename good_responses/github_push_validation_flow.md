# What GitHub does after you click Push

This repo’s validation is driven by:

- [.github/workflows/julia-ci.yml](.github/workflows/julia-ci.yml)
- [test/runtests.jl](test/runtests.jl)
- [test/README.md](test/README.md)

The flow is:

1. You push the commit.
   - GitHub accepts the git update first.
   - A successful push does not guarantee CI passes.

2. GitHub checks the workflow trigger.
   - It reads [.github/workflows/julia-ci.yml].


3. GitHub starts a fresh runner.
   - The job definition is in [.github/workflows/julia-ci.yml], under `jobs` and `runs-on: ubuntu-latest`.

4. GitHub checks out the repo.
   - This is the `actions/checkout` step in [.github/workflows/julia-ci.yml].

5. GitHub installs Julia and dependencies.
   - The workflow runs `setup-julia` and `Pkg.instantiate()` in [.github/workflows/julia-ci.yml].

6. GitHub runs the test suite.
   - The command is:
     `julia --startup-file=no --project=. test/runtests.jl`
   - That file is [test/runtests.jl].

7. GitHub runs smoke checks.
   - Extra jobs in [.github/workflows/julia-ci.yml]run smoke tests such as clean-depot and threaded smoke.
   - These are under [test/smoke](test/smoke).

8. If any step fails, the workflow fails.
   - GitHub marks the check as failed even though the git push already succeeded.

9. GitHub shows the result in the UI.
   - This appears in the Actions tab, commit status, and PR checks.

In short: GitHub reads the YAML workflow, then executes the repo’s Julia test entrypoint and smoke checks. The push is one thing; CI pass/fail is a separate thing.
