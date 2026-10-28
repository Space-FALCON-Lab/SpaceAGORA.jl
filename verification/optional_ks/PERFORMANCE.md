# Bounded matched-accuracy findings

The tested KS adapter is much slower in the qualified two-body timing case.
J2-plus-drag wall-time performance is inconclusive because its monitor contains
long observation gaps. These results support retaining an experimental optional
route, not choosing this implementation as the default. No general speed or
physical-validity claim follows.

## Source, boundaries and achieved accuracy

Numerical sources remain SpaceAGORA c091191a1ae0dc602ae50025889d34d4b2800a96 and
installed KS05d969c24a6595c2752f7ff947af77c504053623, Julia1.12.7, Apple M3Max,
48GiB RAM, one Julia/BLAS thread. Both methods use the same24h physical inputs
and1,441 epochs. The independent refined Cartesian reference and declared
1m/0.001m/s plus mass/heat/endpoint/clock gates are unchanged. Fine KS settings
are relative1e-12, absolute1e-13, maximum physical step5s. DP8 calibration selects
relative3e-9, absolute3e-10, maximum600s for both cases by error closeness, not
speed. All18 calibration candidates remain;10 miss common targets and6 are
more accurate than the matching band. Current KS fine errors are recomputed on
this source and reproduce their historical values exactly.

| Case | KS position(m)/velocity(m/s) | DP8 position(m)/velocity(m/s) | Actual worst error ratio(position/velocity) |
|---|---|---|---|
| Two-body |0.536640386/0.000581055192|0.644833842/0.000628369455|1.201613/1.081428|
| J2 plus drag |0.439535783/0.000477593594|0.375206414/0.000331547088|1.171451/1.440500|

Every actual stable cross-pair meets the factor2 position AND velocity match,
complete86400s endpoint and all common budgets. No accuracy tolerance is adjusted
to hide a disagreement. `stable_execution` in raw worker files certifies only
residual compilation<=1%wall, not quiet or statistically stable machine load.

Timing covers the normal run_simulation call with Cartesian dense solution
retained, results files disabled, including setup/storage. Configuration,
references, sampling and record writing are outside. Package load is separate;
first use uses timed eval to include compilation plus eval overhead. Warm calls
use the direct call. Compilation seconds are reported, never subtracted from
wall time. All3 accepted warm calls per method/case record zero compilation.

## Two-body observed-machine comparison

| Method | Warm median wall(s) | Warm range(s),3 calls | First-use wall(s)/recorded compile(s) | Package-load wall(s) |
|---|---:|---|---|---:|
| DP8 |0.004805334|0.004522375–0.005073750|12.153551375/12.135105916|10.294380750|
| KS |106.051936625|105.865250459–108.972006959|120.462750542/14.597614333|11.985661084|

The observed warm median wall ratio is22,069.63 KS/DP8. It applies only to this
full-route experiment at the specified achieved accuracy. Warm KS allocations
are about197.88GB versus2.22MB for DP8 (allocated bytes across a call, not peak
resident memory); GC is about20.45–21.07s versus zero reported. The monitor's
largest interval is2.211s, consistent with its nominal2s sampling. Three repeats,
millisecond DP8 calls, OSnoise and substantial recorded swap use prevent a
machine-independent or statistically precise ratio claim. The large observed
cost and allocation difference warrant investigation before expanding support.

## J2 plus drag: accuracy passes, wall-time verdict inconclusive

Retain ALL raw calls: KS warm196.230100167,197.079256958 and803.240453167s;
DP8 warm0.011795042–0.012042459s. KS first-use wall7721.427115334s includes
only21.160181329s recorded compilation. The monitor has10 gaps above10s,
maximum5034.119146s, despite a terminal clean process-activity status. That
status does not establish observation coverage across the gaps. No cause is
assigned, no slow call is dropped and no J2 speed ratio is reported. Even the
similar later warm calls cannot rescue a qualified comparison by selecting them.

## Exclusions, provenance and interpretation

Calibration attempt1 was interrupted after a premature launch; attempt2 was
invalidated by another task's active Julia process. Timing's first DP8 launch
failed before spawning Julia because the monitor parent directory was absent.
J2KS timing attempt1 was invalidated when the documented dormant process's CPU
counter advanced0:29.93 to0:29.94. Each failure and retry is preserved. Earlier
valid workers are retained by a prospective selection plan, never by speed.
The per-worker dormant counter remains strictly unchanged in every selected
worker; no other Julia process is accepted. These records do not establish an
otherwise empty machine or detect all shorter-than2s background tasks.

The original seven-script experiment identity remains unchanged. An additive,
independently reviewed retry-aware postprocessor separately hashes its source,
prospective plan and explicit selection manifest; checks source/calibration,
solver settings, exact raw/monitor hashes, launch-directory path association,
all raw physical budgets, ordered repetitions and actual paired errors; and
withholds the J2 ratio because of observation gaps. The10s diagnostic is an
explicit post-measurement evidence qualification, not a retroactive numerical
acceptance change. Original hardcoded-layout summarizer is preserved unchanged.
Raw results, calibration, scripts, monitors, environments, exclusions and reviews
are in the companion archive. See performance-reproduction.md.

Future work should profile the implicit stage/derivative allocation cost and
rerun J2 timing in a continuously observed environment under a prospectively
reviewed policy. That is follow-up research, not a reason to erase this phase's
inconclusive verdict or claim a KS advantage. Integrated sails, multiple craft,
interlinks and physical-time sensitivities stay deferred; MPC is unchanged.
