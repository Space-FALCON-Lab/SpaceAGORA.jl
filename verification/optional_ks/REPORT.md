# Current-upstream optional KS delivery report

The bounded candidate meets the stated fine trajectory and tested physical-time
lifecycle targets. It is reviewable development source, not an adopted release.
Performance and final delivery review remain pending. A supported verdict is
specific to the inputs, settings, source and sampled epochs below.

## Source boundary

SpaceAGORA upstream base is a664365499c08cd84b1421047cb2208bf2b46e5b. The
independent checkpoint-table repair is PR576, production615a447c; its review-only
follow-up a5682903 changes no production. This adapter is a separate dependent
draft. Its numerical source and current-upstream scientific replay are at
c091191a1ae0dc602ae50025889d34d4b2800a96 (guard source0037984e0).

The scientific replay installed KS8cfc936fe2ba403996ed7682aa700b56d6e04c5d.
Complete-candidate CI then found a Julia1.9/1.10 load failure: the error-display
annotation bound Core.IO before the IO facade. Separate KS draft25 qualifies
that annotation as Base.IO. Independent review proves intended dispatch/output
equivalence: all other48 tracked source files are unchanged. Numerical equations,
controls, tolerances and sampling are unchanged. Original scientific source
identities remain explicit; no trajectory rerun is inferred from a version bump.

The installation pin is the complete source05d969c24a6595c2752f7ff947af77c504053623
(version0.2.0, KS draft24). Its complete source/test/Project tree matches the
normal package run with4,170 passing assertions. A genuine installed GitHub
copy matches97 checked source/test/Project files and passes316 adapter checks.
The type-only delta requires no further trajectory rerun according to independent
review. New coupled physical-time sensitivities are not implemented.

## Physics and independent references

The four24h cases use one500kg rigid link, inertial Cartesian SI state, a7,000km
semimajor axis, eccentricity0.01, inclination0.7rad, node0.3rad and periapsis0.2rad.
The specified Earth model has mu3.986004418e14m3/s2 and a spherical6378136.3m
altitude datum; J2 is0.00108262668 about inertial+z, rotation7.292115e-5rad/s.
Drag uses deterministic exponential density, a co-rotating atmosphere, no wind,
constant CD2.2,10m2 area and the frozen attitude/Maxwellian thermal assumptions.
See baseline-inputs.toml for every parameter and heat-unit convention.

The Cartesian reference separately expresses central/J2 acceleration, rotating
exponential drag and the returned heat law; it calls no KS transformation or
SpaceAGORA force function. High-order reference refinement is constrained to1%
of the physical budgets. Two-body additionally agrees with analytical Kepler.
Every comparison uses the same1,441 physical epochs0:60:86400s. No result claims
a uniform continuous error bound, real force-model fidelity or near-collision
validity. Energy conservation is not an acceptance criterion with drag/thrust.

## Trajectory verdicts

| Case | Fine maximum position error(m) | Fine maximum velocity error(m/s) | Verdict |
|---|---:|---:|---|
| Two-body |0.536640386|0.000581055192|Pass|
| J2 only |0.564589426|0.000612789387|Pass|
| Drag only |0.426574111|0.000462111432|Pass|
| J2 plus drag |0.439535783|0.000477593594|Pass|

Fine means relative1e-12, absolute1e-13 and maximum physical step5s, applied to
the phase/auxiliary tolerance fields as the fixture specifies. Targets are1m
and0.001m/s plus mass/heat/clock gates. These settings are case-specific. All
eight coarse/medium rows remain FAIL in evidence/baseline-results.toml; their
maximum position errors range roughly9.2–263.1m. They miss the chosen accuracy
targets and are not automatically another implementation defect.

## Physical-time lifecycle verdicts

The unchanged6000s case uses13s prescribed control,17s output, area changes,
an ideal burn/current-mass propellant flow, a velocity impulse, state/time force,
a descending continuous observer event and997.3s checkpoints. The442s case
coincides impulse/control/output/checkpoint at221s. Separate processes execute
direct, segmented, deliberately interrupted and resumed runs.

All required histories, epochs, mass/heat and physical trajectory gates pass.
An independent SciPy DOP853 implementation recomputed all eight new raw tables
and both interrupted snapshots, with52 additional checks. Main KS direct
maximum errors are0.01093393m and1.19122e-5m/s. Across common physical samples,
segmented/resumed versus direct KS differs by at most5.719e-6m/6.193e-9m/s;
coincident resume differs by at most1.613e-7m/1.756e-10m/s. Independent root
epoch refinement is2.73e-12s against its1e-6s gate.

Current upstream restores valid saved output history. Main resume contains360
rows, coincident resume27, both beginning at0; every nominal epoch is present
and timestamps are unique. The restored60-row prefix through997.3s is bitwise
identical to the segmented prefix. This full-prefix expectation was declared
before the new replay; the older prototype's resume-only row counts were not
forced on current upstream. Missing/incompatible history follows core's separate
fallback contract. No arbitrary callback/event class is qualified by these cases.

## Guards, source review and installation

Independent review found two port guard gaps: newer mounted data owners bypassed
preflight, and a manually resolved inactive central scope was counted. Both were
reproduced before repair and now fail clearly. Normal public scoped-to-spacecraft1
gravity remains admitted. Tests pass4data-owner,6scope,68actual force/heat parity,
20other selected guard/cache and218 dense-history assertions. Independent source
diagnostics add620 assertions and4 root-reference checks. These are distinct
engineering checks, not replacements for the trajectory/lifecycle evidence.

A fresh normal GitHub checkout, including actual LFS assets, ran the committed
installer in separate with-KS and without-KS environments. Runtime316/9 checks
pass with JULIA_LOAD_PATH=@:@stdlib. The complete new KS pin also passes the316
installed checks. Direct Pkg.add of the SpaceAGORA Git URL failed in Julia's
LibGit2 submodule checkout before loading; that failed route is preserved, not
claimed to work. Use the verified checkout plus installer route in README.md.
The KS repository currently requires authorized access; no registry release is
created. OrdinaryDiffEqCore is restricted to the4.19 minor series.

## Ownership and remaining delivery

SpaceAGORA owns configuration, physical state/time, force interfaces/conventions,
scheduling/events/checkpoints/outputs and propagator selection. KSPropagator
remains separate and owns KS coordinates/equations/derivatives, numerical methods
and their verification. Integrated physics comes from the configured SpaceAGORA
RHS. Shared mathematical ownership may follow explicit comparison/review;
Mohammed's MPC prediction method remains unchanged.

The route initially covers one spacecraft and the verified force/lifecycle set.
Integrated sails, multiple spacecraft, interlinks, mounted systems and new
physical-time sensitivities remain deferred. The separate TrajectoryRecorder
cadence and split/multirate paths are outside this table-output qualification.

The bounded matched-accuracy study and final findings review are still pending.
No speed advantage is claimed. The companion upstream-delivery evidence retains
raw CSV/binary results, complete logs, environments, source hashes and failed
attempts; selected verdicts/reviews are included in evidence/. Reproduction
scripts and frozen inputs are in this directory.
