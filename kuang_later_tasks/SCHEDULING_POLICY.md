# How the scheduling policy works

This guide explains the existing SpaceAGORA scheduler used by these experiments,
especially `DataBacklogPolicy` and the `production_competition` case. The fixture
code configures and independently checks the scheduler; it does not introduce a
new scheduling algorithm.

The numerical examples below come from the saved runs described in
[REPORT.md](REPORT.md). Writing this guide did not rerun or change the simulations.
Data units are decimal: $1\,\mathrm{Mbit} = 10^6\,\mathrm{bit}$.

## 1. The main idea: availability, score, then matching

The scheduler answers three different questions:

1. **Allowed:** Was this link registered as a candidate connection?
2. **Available:** Can this registered link operate at the current state?
3. **Selected:** Does activating it belong to the best combination of links?

A link can be allowed and available but still not selected.

For the data fixtures, availability includes range and planet occultation.
Ground-station links instead use their station geometry, including the elevation
mask. Hardware and spacecraft activity can also constrain availability in the
engine, but these simplified fixtures do not establish realistic radio or power
performance.

The scheduling sequence is:

```text
Current spacecraft states and previous active links
    |
    v
Update which registered links are available
    |
    v
Calculate one policy score for each available link
    |
    v
Find the maximum-total-score set with no terminal conflicts
    |
    v
Replace the active-link set and record the selection
    |
    v
Continue simulation, update stores, and repeat
```

Unavailable links score zero. Only available links with strictly positive scores
can be selected.

## 2. A policy and a matching algorithm are different things

The **policy** defines what makes an individual link valuable.

The **matching algorithm** chooses the best compatible combination of those
scored links.

Formally, let $w_e$ be an edge's score, $E_{\mathrm{avail}}$ the set of available
edges, and $M$ the selected set. Selection solves:

$$
\begin{aligned}
\underset{M \subseteq E_{\mathrm{avail}}}{\operatorname{maximize}}
&\quad \sum_{e \in M} w_e \\
\text{subject to}
&\quad w_e > 0 && \text{for every } e \in M, \\
&\quad \sum_{\substack{e \in M \\ v \in e}} 1 \leq 1
&& \text{for every terminal } v.
\end{aligned}
$$

Graph vertices represent **terminals**, not whole spacecraft. A spacecraft with
two distinct terminals can participate in two simultaneous links, provided the
terminal assignments do not conflict. The busy relay uses this arrangement.

### Why this is not just "pick the biggest link first"

Consider this illustrative terminal graph:

$$
P \xleftrightarrow{\;6\;} Q
\xleftrightarrow{\;10\;} R
\xleftrightarrow{\;6\;} S
$$

Choosing the middle edge first gives score 10 and blocks both outer edges.
Choosing both outer edges gives score $6 + 6 = 12$. The scheduler seeks the latter,
maximum-total-score matching.

In the two-sender experiment, however, both links share one destination
terminal. There are only three possible selections: neither link, A's link, or
B's link. With positive scores, choosing the larger score is exactly optimal.

### How the implementation finds the matching

- Tree/forest components use dynamic programming in Julia.
- Cyclic candidate components require the optional PythonCall/NetworkX weighted
  Blossom backend.
- The implementation converts the represented finite `Float64` scores to exact
  integer weights before comparing matching totals. This avoids extra rounding
  in those comparisons; it does not make the underlying physical scores exact.
- Equal-optimal matchings are possible. Tie resolution is deterministic for a
  given backend version, but is not a fairness or round-robin rule.

The backlog fixtures have forest terminal graphs and do not test the cyclic
backend. The forty-spacecraft momentum run verifies disjoint feasible edges;
the ground-station run verifies star-shaped feasible components.

## 3. The full DataBacklogPolicy formula

For an available directed link from sender $i$ to receiver $j$:

$$
S_{i \to j}
= \min\left(
R,\;
\frac{Q_i}{H} + P_i,\;
\frac{F_j}{H} + O_j
\right).
$$

| Symbol | Meaning | Unit |
|---|---|---|
| $R$ | Candidate link's data rate at the current state | $\mathrm{bit/s}$ |
| $Q_i$ | Data currently held by the sender | $\mathrm{bit}$ |
| $H$ | Positive, finite backlog horizon | $\mathrm{s}$ |
| $P_i$ | Sender's current local production rate | $\mathrm{bit/s}$ |
| $F_j$ | Free storage at the receiver | $\mathrm{bit}$ |
| $O_j$ | Receiver's carried outflow from previously active, still-available links, evaluated at the current state | $\mathrm{bit/s}$ |

The minimum represents three possible bottlenecks:

1. **Link capacity:** the connection cannot carry more than $R$.
2. **Sender supply estimate:** stored data spread over $H$, plus ongoing production.
3. **Receiver acceptance estimate:** free space spread over $H$, plus current drainage.

For example, a sender holding 2 Mbit and producing 0.1 Mbit/s, with $H = 5\,\mathrm{s}$,
has a sender supply term of:

$$
\frac{2\,\mathrm{Mbit}}{5\,\mathrm{s}} + 0.1\,\mathrm{Mbit/s}
= 0.5\,\mathrm{Mbit/s}.
$$

If its link supports 1 Mbit/s but the receiver's acceptance term is only
0.3 Mbit/s, the score is $\min(1, 0.5, 0.3)\,\mathrm{Mbit/s} = 0.3\,\mathrm{Mbit/s}$.

### What the horizon does, and does not, mean

$H$ converts a stored amount into a rate-like priority.

- A smaller $H$ makes backlog and free-space terms larger; link-rate caps may
  then dominate and create ties.
- A larger $H$ makes those terms smaller relative to production and outflow.
- $H$ is **not** the solver time step.
- $H$ is **not** how long a selected link is reserved.
- The scheduler does **not** simulate the next $H$ seconds to find the best plan.
- It samples current production and estimates that rate continuing. It does not
  predict future production changes or add incoming relay traffic to $P_i$.

The original backlog cases use $H = 5\,\mathrm{s}$; the ground-station case uses $H = 60\,\mathrm{s}$.

### Empty and full stores

The scoring code respects the engine's empty/full thresholds: $10^{-9}$ times
store capacity. A sender counted as empty contributes zero stored bits; a
receiver counted as full contributes zero free space.

Consequently:

- Empty sender, no production: score zero.
- Empty sender, positive production: potentially positive score.
- Full receiver, no previous feasible outflow: score zero.
- Full receiver that is draining: potentially positive score.

Ground stations are unbounded endpoints in this model, so receiver storage
does not limit a downlink score. A spacecraft without a store is likewise
treated as unbounded for the relevant storage term. All spacecraft in the
original backlog experiments have finite stores.

## 4. Worked example: production_competition

The fixture in [backlog.jl](backlog.jl) creates:

$$
\begin{aligned}
A\;(\text{spacecraft }1)
&\xrightarrow{\;1\,\mathrm{Mbit/s}\;}
\text{destination terminal }(3,1), \\
B\;(\text{spacecraft }2)
&\xrightarrow{\;1\,\mathrm{Mbit/s}\;}
\text{destination terminal }(3,1).
\end{aligned}
$$

| Setting | A | B | Destination |
|---|---:|---:|---:|
| Initial stored data | 0 Mbit | 8 Mbit | 0 Mbit |
| Local production | 0.2 Mbit/s | 0 | 0 |
| Store capacity | 100 Mbit | 100 Mbit | 1,000 Mbit |
| Terminals | 1 | 1 | 1 |

Both links remain available throughout the saved 20-second run. Destination
space is ample, so receiver capacity does not limit either score.

With $Q_A$ and $Q_B$ expressed in Mbit, the following numerical formulas give
scores in Mbit/s:

$$
\begin{aligned}
S_A &= \min\left(1,\frac{Q_A}{5}+0.2\right), \\
S_B &= \min\left(1,\frac{Q_B}{5}\right).
\end{aligned}
$$

### Step 1: B wins initially

$$
\begin{aligned}
S_A &= \min\left(1,\frac{0}{5}+0.2\right) = 0.2, \\
S_B &= \min\left(1,\frac{8}{5}\right) = 1.0.
\end{aligned}
$$

B has the larger score, so B's link is selected. A's positive production does
not override that comparison.

### Step 2: A accumulates data while B drains

Before the first switch, B sends at 1 Mbit/s and A only produces. With $t$
expressed in seconds, backlogs in Mbit, and scores in Mbit/s:

$$
\begin{aligned}
Q_A(t) &= 0.2t, \\
Q_B(t) &= 8-t, \\
S_A(t) &= \min(1,0.04t+0.2), \\
S_B(t) &= \min(1,1.6-0.2t).
\end{aligned}
$$

B's score stays capped at 1 until $t = 3\,\mathrm{s}$, then decreases. A's score increases.
The continuous score crossing on this initial B-only trajectory is:

$$
\begin{aligned}
0.04t+0.2 &= 1.6-0.2t, \\
0.24t &= 1.4, \\
t &= \frac{1.4}{0.24} = 5.833333\ldots
\end{aligned}
$$

The resulting crossing time is $5.833333\ldots\,\mathrm{s}$.

These formulas describe the initial B-only phase, not the entire run.

### Step 3: the recorded first switch is at 6.082066 s

The engine reschedules at initialization and after each accepted solver step.
It does not continuously locate score-equality events.

| Accepted time (s) | A score (Mbit/s) | B score (Mbit/s) | Selection |
|---:|---:|---:|---|
| 5.832066 | 0.433283 | 0.433587 | B |
| 6.082066 | 0.443283 | 0.383587 | A |

The first sample is just before the crossing. At the next scheduling sample,
A has overtaken B. Thus the recorded switch occurs at 6.082066 s, not exactly
5.833333 s.

### Step 4: sending changes the scores, so the winner changes again

While A has backlog and is selected, it sends at 1 Mbit/s but produces only
0.2 Mbit/s:

$$
\begin{aligned}
\frac{dQ_A}{dt}
&= (0.2-1)\,\mathrm{Mbit/s}
= -0.8\,\mathrm{Mbit/s}, \\
\frac{dS_A}{dt}
&= \frac{-0.8\,\mathrm{Mbit/s}}{5\,\mathrm{s}}
= -0.16\,\mathrm{Mbit/s^2}.
\end{aligned}
$$

The score derivative applies while the sender-supply term is uncapped and
the sender still has backlog.

B is waiting and has no production, so its score stays constant.
Eventually B wins again.

While B is selected:

$$
\begin{aligned}
\frac{dS_B}{dt}
&= \frac{-1\,\mathrm{Mbit/s}}{5\,\mathrm{s}}
= -0.2\,\mathrm{Mbit/s^2}, \\
\frac{dS_A}{dt}
&= \frac{0.2\,\mathrm{Mbit/s}}{5\,\mathrm{s}}
= 0.04\,\mathrm{Mbit/s^2}.
\end{aligned}
$$

Again, these are the uncapped score derivatives.

That reverses their relative priorities again. This alternation comes from
rescoring changing backlogs, not a programmed round-robin rotation.

### Step 5: actual selection intervals

These intervals are from the saved
[link schedule](results/backlog/production_competition/link_schedule.csv).
Times below are rounded to six decimal places; intervals run from their start
up to the next boundary.

| Start (s) | End (s) | Selected sender |
|---:|---:|---|
| 0 | 6.082066 | B |
| 6.082066 | 6.582066 | A |
| 6.582066 | 6.832066 | B |
| 6.832066 | 7.082066 | A |
| 7.082066 | 7.332066 | B |
| 7.332066 | 7.832066 | A |
| 7.832066 | 8.082066 | B |
| 8.082066 | 8.332066 | A |
| 8.332066 | 8.582066 | B |
| 8.582066 | 20 | A |

### Step 6: why B's last data does not get sent

After B's last selected interval, it still holds approximately 0.917934 Mbit:

$$
S_B = \frac{0.917934\,\mathrm{Mbit}}{5\,\mathrm{s}}
\approx 0.183587\,\mathrm{Mbit/s}.
$$

A's score cannot fall below 0.2 Mbit/s in this scenario, even when A is empty:

$$
S_A = \frac{0\,\mathrm{Mbit}}{5\,\mathrm{s}}
+ 0.2\,\mathrm{Mbit/s}
= 0.2\,\mathrm{Mbit/s}.
$$

A therefore stays selected for the rest of the run. At approximately
8.852582 s its backlog empties, and it then sends its continuing production.

More generally, under these ample-receiver-space and equal-link assumptions,
an idle nonproducing B loses to an empty producing A whenever:

$$
\frac{Q_B}{H} < P_A
\quad\Longleftrightarrow\quad
Q_B < P_A H.
$$

Here that threshold is $(0.2\,\mathrm{Mbit/s})(5\,\mathrm{s}) = 1\,\mathrm{Mbit}$. Equality would be a tie, not a
guaranteed switch. Continued starvation beyond this run would depend on the
same geometry, production, and capacity conditions continuing.

The final recorded totals are:

- A sends approximately 4 Mbit, its total production.
- B sends approximately 7.082066 Mbit.
- B retains approximately 0.917934 Mbit.
- The destination receives approximately 11.082066 Mbit.
- No production loss occurs to the fixture's tolerance.

## 5. A scheduling score is not a commanded transfer rate

This is the most important distinction when reading the results.

| Quantity | Purpose |
|---|---|
| Link rate | Capacity supplied by the link model |
| Policy score | Weight used to choose the active matching |
| Actual carried rate | Flow the active link can carry given store and link constraints |

At A's first selection, its score is about 0.443283 Mbit/s, but it has enough
backlog to transmit at the full 1 Mbit/s link rate.

Once A is empty, the same selected link carries only 0.2 Mbit/s because that is
all A produces. The scheduler does not automatically use the spare capacity
for B: the destination's single terminal is still assigned to A.

The runtime flow calculation separately enforces empty/full-store constraints,
including relay inflows and outflows. The score is a selection heuristic, not
that runtime flow solution.

## 6. Previous feasible outflow: why the busy relay is different

For an incoming link to a full relay, $F_j = 0$. Its receiver term becomes
the relay's previous feasible outflow $O_j$.

The engine evaluates this outflow:

1. After refreshing current availability.
2. From links that were active before the new selection and remain available.
3. Using current-state carried-flow caps, not just their nominal capacities.
4. Once for the scoring pass, before choosing a new matching.

It does **not** assume an outgoing link will be selected merely because that
would make an incoming link score better. It also gives no credit for a
previous outgoing link that has gone out of range.

In the busy-relay fixture:

| Scheduling situation | Relay input score | Direct score | Relay output score | Winning combination |
|---|---:|---:|---:|---|
| Initial state, no previous output | 0 | 0.5 | 1 | Direct + relay output: 1.5 |
| Next step, previous feasible output | 1 | 0.5 | 1 | Relay input + relay output: 2 |

Scores are in Mbit/s. Source-to-relay input and the direct route share the
source's terminal, so they cannot both be selected. Relay input and output use
different relay terminals, so they can coexist.

The second combination has a larger sum of link scores, but only relay output
delivers to the destination. Relay input is an intermediate transfer.
Moreover, at the full relay, incoming data displaces its own production, which
overflows.

That is why the direct-route control can deliver more and lose less despite
the busy-relay run maximizing its instantaneous matching score.

**Maximum-weight matching is exact for the supplied scores. Those scores are
not an end-to-end routing objective, a no-loss guarantee, or a packet-origin
tracking system.**

## 7. Scheduling time, solver accuracy, and switching

The scheduling callback runs at initialization and after accepted integration
steps. A selection recorded at $t_k$ applies to the following interval until
the next scheduling update. A final-time selection has no subsequent interval
inside the run.

The original data cases use a maximum step of 0.25 s. Adaptive integration can
take smaller steps, including near store events. This is not a fixed 0.25-second
scheduling timer.

The saved [step-refinement comparison](results/backlog/step_refinement.csv)
shows:

| Competition run | Maximum step | Selection changes | Final delivery |
|---|---:|---:|---:|
| Original | 0.25 s | 9 | 11.082065734343070 Mbit |
| Refined | 0.125 s | 19 | 11.082065734343072 Mbit |

Delivery is essentially identical here, while switching differs significantly.
That is evidence about these two runs, not a general convergence guarantee.

Tight solver error tolerances do not make score crossings or visibility
boundaries continuously located scheduling events. The backlog policy also
contains no switching penalty, minimum dwell time, or hysteresis.

## 8. What this policy does not guarantee

- **Fairness:** no waiting-time bonus, round-robin quota, or equal airtime.
- **Complete drainage:** a smaller nonproducing backlog can remain unserved.
- **Maximum destination delivery:** intermediate transfers also earn scores.
- **No loss:** production can overflow a full store.
- **Future optimality:** no look-ahead over future contacts, source changes, or routes.
- **Stable switching frequency:** rescoring and step resolution affect switches.
- **Hardware realism:** the fixtures use fixed-rate links and simplified terminals,
  not validated RF budgets, acquisition delays, or battery constraints.

Adding fairness, deadlines, switching costs, or destination-aware routing would
require an explicitly different objective or scheduling design. Those features
are not implemented by these fixtures.

## 9. How the other experiments fit

### Empty producer control

With no competitor, an empty source producing 0.2 Mbit/s has a positive score,
remains selected, and delivers 4 Mbit over 20 seconds.

### Ground-station downlinks

The ground receiver is unbounded, so a feasible downlink's score reduces to:

$$
S_{i \to \mathrm{ground}}
= \min\left(R,\frac{Q_i}{60\,\mathrm{s}}+P_i\right).
$$

Each station has one antenna. In the verified geometry, a spacecraft sees at
most one station at once, so each feasible component is a star. Selection
therefore picks one maximum-positive-score link per station.
Visibility alone does not guarantee selection.

See [GROUND_REPORT.md](GROUND_REPORT.md) for the separate geometry and accounting.

### Momentum links

The forty-spacecraft momentum experiment uses a different policy:
`SemimajorAxisRatePolicy(targets=31:40)`.

Its score measures the target semimajor-axis rate:

$$
\frac{da}{dt}
= \frac{2a^2\,\mathbf{v}\cdot\mathbf{F}}{\mu m}.
$$

Here $a$ is semimajor axis, $\mathbf{v}$ is target velocity, $\mathbf{F}$ is
the force on the target, $\mu$ is the planet's gravitational parameter, and
$m$ is target mass.

The units are m/s, not bit/s. The same matching machinery chooses compatible
positive-score links, but backlog, production, and storage are irrelevant to
that objective. Raw scores from different physical units should not simply be
added without a deliberately scaled combined objective.

## 10. Code and evidence map

| Item | Where to read |
|---|---|
| Policy definition, backlog formula, scoring, and matching | [interlink_scheduling.jl](../src/models/gnc/control/interlinks/interlink_scheduling.jl) |
| Initialization and accepted-step scheduling | [interlink_scheduler.jl](../src/orchestration/simulation/callbacks/interlink_scheduler.jl) |
| Actual carried data rates and store limits | [terminals.jl](../src/models/vehicle/interlinks/terminals.jl) |
| Experiment setup, independent scores, exhaustive matching checks | [backlog.jl](backlog.jl) |
| Shared configuration and matching assertions | [run.jl](run.jl) |
| Scores for both selected and rejected competition links | [candidate_scores.csv](results/backlog/production_competition/candidate_scores.csv) |
| Every competition scheduling decision | [accepted_steps.csv](results/backlog/production_competition/accepted_steps.csv) |
| Competition contact and active intervals | [link_schedule.csv](results/backlog/production_competition/link_schedule.csv) |
| Competition store contents and counters | [node_accounting.csv](results/backlog/production_competition/node_accounting.csv) |
| Previous-outflow probes | [outflow_probes.csv](results/backlog/outflow_probes.csv) |
| Results, assumptions, and tolerances | [REPORT.md](REPORT.md) |

For CSV interpretation, `1.1->3.1` means spacecraft 1, terminal 1, transmitting
to spacecraft 3, terminal 1. Thus A is `1.1->3.1` and B is `2.1->3.1`.
An `available`/contact record describes a candidate opportunity; an `active`
record describes the selected connection. `changed=true` marks a selected-set
change, not necessarily a change in visibility.

**Summary:** the backlog policy estimates each link's usefulness from current
capacity, stored data, production, receiver room, and previous feasible outflow.
The matching algorithm maximizes the sum of those scores without terminal
conflicts. The simulation then determines actual data flow, and the next
accepted step starts the process again.
