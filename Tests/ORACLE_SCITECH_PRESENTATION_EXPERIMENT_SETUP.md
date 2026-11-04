# Experiment Setup in the ORACLE SciTech Presentation

## Source and scope

Source: the supplied 35-page **ORACLE Scitech Presentation.pdf**, titled
*Assessing Small Satellite Maneuverability in Proliferated LEO Constellations
with Open Cavity Laser (OCL) Interlinks*, by Nihal Simha, Kuang Sun, and
Giusy Falcone, for AIAA SciTech Forum 2027.

This note summarizes the **presentation**, not an assumed configuration of
any particular repository snapshot. Page references below are PDF page
numbers, which match the printed slide numbers on the cited slides.
Both extracted text and rendered plots were inspected, including plot titles,
axes, and legends that were absent from the extracted text.

Evidence is distinguished as follows:

- **Stated:** explicitly written in slide text or a parameter table.
- **Visible in plots:** recoverable from titles, legends, or axes.
- **Derived:** calculated or interpreted from those settings.
- **Unspecified:** not established by the presentation; requires source
  configuration or original data to reproduce exactly.

## 1. Main study in one paragraph

The main experiment is a **single-target maneuverability parameter sweep**
using range-limited OCL links to an evenly distributed equatorial helper
population at **1000 km altitude**. It varies helper count, target altitude,
and target inclination, and compares two separately applied instantaneous
objectives: **increasing semimajor axis** and **increasing inclination**.
Each case runs for **80 orbits, approximately six days**. Outputs are signed
cumulative laser-induced velocity increments in the target's radial,
transverse, and normal frame, and changes in orbital elements. A separate
two-spacecraft case checks the mechanical balance laws.

The presentation is not a three-code comparison and does not identify
which of the local numbered code snapshots produced the plots.

## 2. Common spacecraft and laser parameters

The parameter table on **slide 14** states:

| Parameter | Value | Interpretation |
|---|---:|---|
| Maximum laser-link range | 200 km | Distance between linked spacecraft, not their altitude difference |
| Laser input power | 10 kW | The force model uses pair power |
| OCL amplification factor | 100 | Denoted by $B$ |
| Geometric transmission/loss factor | 1.0 | Denoted by $\beta$; no attenuation from this factor |
| Spacecraft mass | 227 kg | Presented as a common spacecraft parameter |

Slides 6-9 model translational motion using position and velocity states,
central gravity, a $J_2$ acceleration term, and switched OCL force:

$$
\dot{\mathbf{r}}_i=\mathbf{v}_i,\qquad
\dot{\mathbf{v}}_i=
\mathbf{a}_{i,\mathrm{Kep}}+\mathbf{a}_{i,J_2}
+\sum_{j\ne i}K_{ij}\frac{\mathbf{F}_{i\leftarrow j}}{m_i}.
$$

The slide-8 force law is:

$$
\mathbf{F}_{i\leftarrow j}
=\frac{\beta_{ij}B P_{ij}}{c}\hat{\boldsymbol{\rho}}_{i\leftarrow j},
\qquad
\hat{\boldsymbol{\rho}}_{i\leftarrow j}
=\frac{\mathbf{r}_i-\mathbf{r}_j}{\lVert\mathbf{r}_i-\mathbf{r}_j\rVert}.
$$

For an active pair, the forces are equal and opposite and directed along
the separation line. Both spacecraft respond; helpers should not be treated
as unaffected prescribed trajectories when reproducing these dynamics.

**Derived scale:** taking $c=299{,}792{,}458\ \mathrm{m/s}$, the listed
parameters give about **0.003336 N (3.336 mN)** per endpoint and an acceleration
magnitude of about **$1.46945\times10^{-5}\ \mathrm{m/s^2}$** at 227 kg.
These are derived values, not additional settings printed in the table.
Do not insert another factor of two merely because the force is optical:
the displayed law is $\beta BP/c$.

### Environment qualification

The presentation's governing model includes $J_2$, and the verification
plot on slide 15 explicitly says **“J2 active.”** This is strong evidence
for the intended environment model. However, individual sweep figures do
not expose their runtime $J_2$ flag. Exact reproduction should confirm that
flag rather than substitute the no-$J_2$ setting of our recent comparisons.

No drag or solar-radiation-pressure term appears in the displayed equations.
Their absence from those equations is not a complete record of every
runtime option. Numerical values for Earth radius, $\mu$, and $J_2$ are
not supplied.

## 3. Main geometry and sweep

### Stated setup on slide 17

| Quantity | Stated setting |
|---|---|
| Number of targets | One |
| Helper orbit altitude | Fixed at 1000 km |
| Helper inclination | Fixed at 0 degrees |
| Helper population | Varied from 1 to 300, evenly distributed |
| Target altitude | Varied from 850 to 1150 km |
| Target inclination | Varied from 0 to 1.5 degrees |
| Duration | 80 orbits, approximately six days |

“Fixed” here describes the helper population's configured shell/plane, not
evidence that helper states are held fixed during integration.

The motivation mentions multi-shell Walker-Delta constellations, but the
explicit experiment uses an equatorial helper population and a single
target whose shell and plane are varied. No complete Walker specification
(total satellites, plane count, and phasing parameter) is supplied.

### Discrete values actually displayed

The plots on **slides 18-21**, reinforced by slides 22-26, show:

| Axis or grouping | Visible values |
|---|---|
| Target altitude columns | 1150, 1050, 1000, 950, 850 km |
| Initial inclination-offset curves | 0, 0.5, 1.0 degrees |
| Helper counts in RTN plots | 1, 50, 100, 150, 200, 250, 300 |
| Helper counts in orbital-element plots | 50, 100, 150, 200, 250, 300 |
| Policies | SMA-rate maximizing and inclination-rate maximizing |

The legend uses $\Delta i$ for the **initial target/helper inclination
offset**. Because helpers are equatorial, those offsets correspond to
target inclinations of 0, 0.5, and 1 degree. This is distinct from the
output $\Delta i=i_f-i_0$ plotted on the inclination-change axes.

**Important discrepancy:** slide 17 states an inclination range through
1.5 degrees, but the displayed legends show only 0, 0.5, and 1 degree.
The PDF does not establish whether 1.5-degree runs were performed but omitted,
or whether the setup text was not updated.

**Conditional case count:** interpreting the visible RTN plots as a full
cross-product gives $5\times3\times7\times2=210$ parameter/policy combinations.
The orbital-element panels display the six-count subset, or 180 combinations.
These counts are inferred from the plots, not an explicit run manifest;
the two groups can represent outputs from the same underlying runs.
Do not interpret “1 to 300” as proof that every integer helper count was run.

### Initial eccentricity and phase

Plot titles on **slides 22-24, 26, and 27** include:

$$
e=0.0,\qquad \nu=1.0^\circ.
$$

This supports a circular initial orbit and an initial anomaly parameter of
1 degree for those examples. In context these appear to describe the target,
but the headings do not give separate complete initial-element records for
target and helpers. Applying them to every sweep cell is an inference.

For a circular orbit, true anomaly relative to periapsis is conventional:
the actual position also depends on the chosen argument of periapsis and
node. The slides do not specify those angles or the helper-ring phase origin.
In particular, **do not assume the target began at 0 degrees**, as in our
recent runs.

## 4. Link constraints and scheduling rules

Slide 8 defines symmetric binary activation, $K_{ij}=K_{ji}\in\{0,1\}$,
limited by pair range. Slide 9 states a one-to-one link assumption:

$$
\sum_{j\ne i}K_{ij}\leq1.
$$

For the single-target case, slide 16 selects the feasible helper having
the largest **positive instantaneous** objective rate. If no feasible helper
has a positive rate, the target does not activate a link.

### SMA-rate objective

$$
S_j^{(a)}=\dot a_j=
\frac{2}{n\sqrt{1-e^2}}
\left(e\sin\nu\,A_{R,j}+\frac{p}{r}A_{T,j}\right).
$$

### Inclination-rate objective

$$
S_j^{(i)}=\dot i_j=
\frac{r\cos u}{na^2\sqrt{1-e^2}}A_{N,j}.
$$

Here $A_R,A_T,A_N$ are components of the candidate perturbing acceleration;
$p=a(1-e^2)$, $n=\sqrt{\mu/a^3}$, and $u=\omega+\nu$.
The slides use orbital-element **GVE form**, rather than the Cartesian
vector-rate form now used in our v9 implementation.

The slide equations describe physical rates. They do not prove whether
the producing code evaluated actual acceleration or ranked unit directions.
With equal candidate force magnitudes and one target mass, either approach
has the same ranking when the positive scaling is common.

These are instantaneous selection objectives, not a demonstrated global
optimum over the full 80-orbit trajectory. The slides do not specify callback
cadence, tie-breaking, terminal-level graph matching, or extra visibility,
pointing, activity, power-budget, or thermal constraints. A 200 km altitude
difference alone is insufficient to establish availability: actual
three-dimensional spacecraft separation matters.

## 5. Measured quantities

Slide 17 defines two output families.

### Signed RTN velocity increments

$$
\Delta v_d(t)=\int_0^t A_d(\tau)\,d\tau,\qquad d\in\{R,T,N\}.
$$

Slides 18 and 20 label these as OCL contributions and plot them in **m/s**.
They are signed component integrals, so opposing contributions can cancel.
Although slide 17 calls them “cumulative impulse,” they are more precisely
**accumulated velocity increments**, not force impulse in N s and not
propulsion $I_{\mathrm{sp}}$ in seconds.
Nor are they necessarily the components of the endpoint difference in
inertial velocity: the RTN frame rotates during the orbit.

Inclination raising does not require a positive unweighted
$\int A_N\,dt$: the physical inclination rate weights $A_N$ by $r\cos u/h$.
This explains why a signed normal-velocity integral can be negative in an
inclination-maximizing case without itself proving a scheduler error.

### Orbital-element differences

$$
\Delta q=q_f-q_0,\qquad q\in\{a,e,i,\Omega\}.
$$

Slides 19 and 21 display SMA change in metres, eccentricity change
(dimensionless), inclination change in degrees, and RAAN change in degrees.

The presentation does not say these are differences **against a no-laser
baseline**. Under its written definition, they are end-minus-start changes
and may include gravitational perturbation effects. The distinction matters
especially for RAAN when $J_2$ is active.

The RAAN panels show nonzero-inclination curves near the upper 300-degree
range. Without the angle-unwrapping convention, one should not interpret
those numbers as literal hundreds of degrees of laser-induced rotation:
modulo-$360^\circ$ subtraction can turn a small negative change into a large
positive number. RAAN is also undefined at exactly zero inclination.
These are interpretation cautions, not a diagnosis of the plotting code.

## 6. Separate mechanical-verification experiment

**Slide 15 is not the 80-orbit parameter sweep.** It states:

| Setting | Evidence |
|---|---|
| Helper | One spacecraft at 1000 km |
| Target | One spacecraft at 1100 km |
| $J_2$ | Plot explicitly labeled active |
| Displayed duration | Plot x-axis spans approximately 0-50 orbits |
| Diagnostic shown | Per-spacecraft change in z angular momentum |

Its duration is inferred from the plotted axis; the slide does not give a
separate exact time specification or the initial inclinations, phases,
and scheduler choice for this case.

It checks net internal OCL force, the summed OCL angular-momentum change,
and energy change minus optical work. The slide reports:

- Net OCL force: zero to the displayed precision.
- Angular-momentum residual: at most $1.30\times10^{-12}\ \mathrm{kg\,m^2/s}$.
- Energy residual: at most $6.08\times10^{-1}\ \mathrm{J}$.

These are the presentation's reported values, not results independently
recomputed here. Zero internal OCL force/torque does not mean the total
constellation momentum vector is constant under external gravity.

## 7. What the highlighted figures compare

| Slides | Comparison shown |
|---|---|
| 18-19 | RTN and orbital-element outcomes under SMA-rate scheduling |
| 20-21 | The same output families under inclination-rate scheduling |
| 22 | SMA scheduling and transverse authority; same-shell geometry highlighted |
| 23 | SMA scheduling and radial authority; shell separation highlighted |
| 24 | Inclination scheduling and normal authority; plane separation highlighted |
| 26, left | SMA policy, target at 1000 km, helper-count effect on transverse increment |
| 26, right | Inclination policy, target at 850 km, helper-count effect on normal increment |
| 27 | Policy comparison across target altitudes |

The two panels on slide 26 use **different target altitudes**, so they are
not a controlled same-geometry comparison of scheduler performance.
Its normal-component red curve is also not monotonic with helper count:
the general saturation conclusion should not be read as a guarantee that
every signed output increases with each additional helper.

Slide 27 is headed “Orbital element change,” but the embedded axes are
$\Delta v_R$ and $\Delta v_N$ in m/s. Those panels should be described as
RTN velocity-increment results. The actual orbital-element panels are
slides 19 and 21.

## 8. Differences from our recent 50-helper experiments

| Setting | Presentation | Recent local comparison runs |
|---|---|---|
| Duration | 80 orbits | 50 initial target periods |
| Helpers | Sweep of populations up to 300 | 50 |
| Target altitude | Five displayed values from 850 to 1150 km | 1050 km |
| Target inclination | Several offsets; displayed 0, 0.5, 1 degree | 0 degrees for SMA experiment; 1 degree for inclination experiment |
| Example anomaly parameter | 1 degree in plot titles | 0 degrees |
| Gravity | $J_2$ in model; explicitly active in verification | Central gravity only; $J_2$ disabled |
| Main purpose | Geometry and policy effects on maneuverability | Agreement of prototype, v1, and v9 |
| Outputs | Full RTN increments and changes in $a,e,i,\Omega$ | Final orbital metrics and helper-contact histories |

The common 1000 km helper altitude, 200 km range, 10 kW power, factor 100,
and 227 kg mass do not make these identical experiments.
Our existing values should therefore not be expected to reproduce the
slide curves directly or by simply multiplying by $80/50$; phasing,
perturbations, and changing encounter geometry affect the integrated result.

Related local reports:
[SMA experiment](SMA_SCHEDULER_COMPARISON_50_HELPERS_50_ORBITS.md) and
[inclination experiment](INCLINATION_COMPARISON_50_HELPERS_50_ORBITS.md).

## 9. Missing information needed for an exact reproduction

1. **Run manifest:** exact discrete sweep, whether 1.5-degree cases exist,
   and which panels omit the single-helper case.
2. **Initial states:** helper eccentricities, all phase offsets, target/helper
   RAAN and argument-of-periapsis conventions, epoch, and reference frame.
3. **Duration convention:** which orbit period defines 80 orbits, whether
   duration changes with target altitude, and whether counting uses an
   initial period or completed revolutions.
4. **Dynamics flags/constants:** confirmation of $J_2$ for every sweep run,
   Earth constants, other enabled forces, and mass treatment.
5. **Numerical settings:** solver, tolerances, maximum step, event detection,
   scheduler cadence, and continuous-force versus discrete-kick application.
6. **Link feasibility:** line of sight, pointing/acquisition, minimum
   separation, and treatment of coincident or co-orbital spacecraft.
7. **Policy edge cases:** equal-score tie-breaking and equatorial
   inclination-rate treatment.
8. **Post-processing:** RTN integration method, angle wrapping/unwrapping,
   osculating versus mean elements, and whether a control trajectory is used.
9. **Provenance:** code revision, input files, raw results, and environment.

**Bottom line:** the slides establish an 80-orbit, single-target,
range-limited OCL maneuverability sweep and a separate two-spacecraft
balance-law verification. They give the principal physical parameters and
selection objectives, but not a complete executable specification.
