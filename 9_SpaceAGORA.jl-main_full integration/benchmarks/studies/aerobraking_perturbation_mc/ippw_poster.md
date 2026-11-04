# IPPW Poster Presentation Plan
## *Perturbation-Assisted Aerobraking: Dimensionless Regimes and Equilibrium Periapsis Corridors Across Celestial Bodies*

Evan D. Yu and Giusy Falcone — University of Michigan

---

## 1. The Central Argument (What the Poster Must Prove)

The poster answers a single question:

> **When do gravitational perturbations help—rather than hurt—an aerobraking campaign, and how can mission designers predict this in advance across any target body?**

The poster's claim is: a unified, dimensionless perturbation framework predicts whether a given orbit geometry places the spacecraft in a drag-assisting, drag-fighting, or near-equilibrium regime, validated by high-fidelity numerical survey across Mars, Venus, Earth, and Titan. This claim is new because prior work treats gravitational perturbations as error sources to be bounded, not as designable resources to be exploited.

---

## 2. Exact Contributions

These are stated precisely, not as activities.

**Contribution 1 — Dimensionless perturbation parameter set** $(\Pi_{3b}, \Pi_{J_2}, \Pi_H, \Pi_D)$
A compact, cross-body nondimensionalization of the four dominant perturbation families (third-body tides, $J_2$ oblateness, higher-order harmonics, and drag) that normalizes each by local primary gravity. The contribution is not the individual parameters—some are known—but the unified set that enables direct cross-body comparison and regime classification on a single dimensionless axis.

**Contribution 2 — Periapsis evolution scaling law with equilibrium condition**
An analytical expression:
$$\Delta r_p \sim -a\,\Pi_D\,\mathcal{F}_D(e) - a\bigl[\Pi_{J_2}\,\Xi_{J_2} + \Pi_H\,\Xi_H + \Pi_{3b}\,\Xi_{3b}\bigr]$$
that links each perturbation family to mean per-orbit periapsis change through dimensionless parameters and geometry factors. The equilibrium condition $\Delta r_p \approx 0$ defines corridors where gravitational perturbations balance drag-driven descent. The contribution is the structural insight that equilibrium corridors exist and are predictable from the $\Pi$-parameters alone, without running a full campaign simulation.

**Contribution 3 — Multi-body numerical regime survey and regime maps**
A systematic sampling of $(a, e, i, \Omega, \omega)$ across representative aerobraking corridors for Mars, Venus, Earth, and Titan using SpaceAGORA.jl, producing regime maps that classify each orbit as drag-dominant (net lowering), gravitational-assist (net raising or equilibrium), or hazardous (accelerated lowering). The contribution is the maps themselves and the validation that analytical trends from Contributions 1–2 correctly predict the numerically observed regime boundaries.

**Contribution 4 — Design implications and cross-body regime comparison**
A demonstration that the regime maps unify qualitatively familiar but previously disconnected phenomena: solar-tide corridor effects at Venus, gravity-anomaly-driven excursions at Mars, and strong third-body influence at Titan and high-apoapsis Earth regimes. The contribution is the interpretive framework that lets a mission designer read regime maps for a new body and extract actionable guidance on periapsis targeting conservatism and maneuver frequency.

---

## 3. What Results Must Appear on the Poster

Each result below is tied to a specific contribution and a specific poster element.

### R1. The $\Pi$-parameter table (Contribution 1)

A four-row comparison table giving the closed-form expression, evaluation point, and physical interpretation for each parameter. This should also include a numerical column showing representative values for Mars, Venus, Earth, and Titan at a reference aerobraking orbit. This makes the abstraction concrete immediately.

| Parameter | Expression | Evaluated at | Representative range |
|---|---|---|---|
| $\Pi_{3b}$ | $3(\mu_k/\mu_P)(r_a/r_k)^3$ | Apoapsis | Venus: ~1×10⁻⁶; Titan: ~1×10⁻³ |
| $\Pi_{J_2}$ | $(3/2)J_2(R_P/a)^2$ | Semimajor axis | Mars: ~2×10⁻³ |
| $\Pi_H$ | $\max_{\ell\geq3}|\bar{C}_\ell|(R_P/a)^\ell$ | Semimajor axis | Body-specific |
| $\Pi_D$ | $\rho_p r_p / \beta$ | Periapsis | Campaign-dependent |

The numerical values must be computed and filled in for each body at a representative orbit from the survey.

### R2. The periapsis evolution equation and equilibrium diagram (Contribution 2)

Two parts:
- The scaling equation $\Delta r_p \sim -a\Pi_D\mathcal{F}_D(e) - a[\Pi_{J_2}\Xi_{J_2} + \Pi_H\Xi_H + \Pi_{3b}\Xi_{3b}]$ displayed prominently as the central analytical result.
- A schematic diagram showing the equilibrium condition graphically: a 2D plot of the drag term vs. the total gravitational term, with the equilibrium line $\Pi_D\mathcal{F}_D = -[\ldots]$ drawn. Points above the line are drag-dominant (periapsis drops), points below are gravitational-assist (periapsis stable or rises). This is the single most communicable result on the poster.

### R3. Regime maps for each target body (Contribution 3) — the core numerical result

For each of Mars, Venus, Earth, and Titan: a 2D colored scatter plot or heatmap in a representative two-parameter slice (e.g., $\omega$ vs. $e$, or $\omega$ vs. inclination, holding other elements near campaign-typical values), where each point is colored by regime:
- **Blue**: net raising / gravitational assist ($\Delta r_p > 0$)
- **Red**: net lowering / drag-dominant ($\Delta r_p < 0$)  
- **Green/White**: near-equilibrium ($|\Delta r_p| < \epsilon$)

These four panels are the centerpiece of the poster. They show the regime structure visually and should be arranged in a 2×2 grid. Overlaid contours from the analytical prediction (Contribution 2) should match the numerically computed boundaries, validating the theory.

**Critical**: the regime maps must include the equilibrium corridors as well-defined loci, not just scattered points. The visual take-away is that equilibrium corridors exist and are orbit-geometry-dependent.

### R4. Per-orbit $\Delta h_p$ as a function of $\omega$ for a reference orbit (Contributions 2+3 combined)

A line plot for one reference orbit at each body showing $\Delta h_p(\omega)$ as AOP cycles through $0°$–$360°$. This directly illustrates:
- The sinusoidal periapsis evolution under $J_2$-only ($A\sin\omega$ shape at Mars)
- The $\sin(2\omega)$ modulated ZLK direct channel at Venus/Titan
- The crossing from positive to negative at specific AOP values — precisely the equilibrium corridor boundaries

This plot is easy to read and directly connects the equation to the operational picture of a campaign evolving through an AOP cycle.

### R5. Cross-body regime comparison summary (Contribution 4)

A simple visual element — either a 2×4 comparison table or a labeled axis — placing all four bodies on a spectrum from "drag-dominant" to "gravitational-assist-dominant" at their characteristic aerobraking orbits. This should directly reference the three previously qualitative examples now unified by the framework:
- Venus: solar-tide ZLK dominant → explained by large $\Pi_{3b}$, near-zero $f_\text{geoid}$
- Mars: $J_2$-dominant indirect channel → explained by large $\Pi_{J_2}$, oblate geoid
- Titan: Saturn-tidal dominant, campaign-timescale ZLK → largest $\Pi_{3b}$ of the four
- Earth: intermediate; both channels relevant

This panel closes the loop between the abstract framework and the familiar operational history, giving the audience an "aha" moment.

---

## 4. Poster Layout Plan

A standard IPPW poster is typically A0 or 36×48 inches, landscape or portrait. The following layout assumes portrait A0 (841 × 1189 mm), divided into a header and 3 rows of panels.

```
┌───────────────────────────────────────────────────────────┐
│  HEADER: Title / Authors / University / Abstract sentence │
│  [optional: visual hook — AOP sinusoid sketch or globe]   │
├──────────────┬──────────────┬──────────────────────────────┤
│  ROW 1       │              │                              │
│  Motivation  │  Framework   │  Central Equation +          │
│  + Gap       │  Schematic   │  Equilibrium Diagram (R2)    │
│  (3 bullets) │  (force      │                              │
│              │  diagram)    │                              │
├──────────────┴──────────────┴──────────────────────────────┤
│  ROW 2 — REGIME MAPS (centerpiece, full width, 2×2 grid)  │
│                                                            │
│  [Mars]  [Venus]       [Earth]  [Titan]                   │
│   (Δhp color map, ω vs e or ω vs i, analytical overlay)  │
│                                                            │
├──────────────┬──────────────┬──────────────────────────────┤
│  ROW 3       │              │                              │
│  Δhp vs ω   │  Π-parameter │  Design Implications +       │
│  line plot   │  table (R1)  │  Cross-body comparison (R5)  │
│  (R4)        │              │                              │
├──────────────┴──────────────┴──────────────────────────────┤
│  FOOTER: Key takeaway sentence + QR code + Acknowledgments │
└───────────────────────────────────────────────────────────┘
```

**Design notes:**
- The regime maps (Row 2) should occupy ~40% of the poster height. They are the primary result and the most visually striking element.
- The equilibrium diagram (Row 1, right panel) and the $\Delta h_p(\omega)$ plots (Row 3, left panel) are the two analytical "teaching" visuals that support the maps.
- The $\Pi$-parameter table (Row 3, center) should fit in one column and be dense but readable — it anchors the framework numerically.
- Color scheme: use a single consistent diverging colormap (e.g., blue–white–red) for the regime maps, referenced in a shared legend. Avoid per-panel color schemes.

---

## 5. How Each Result Should Be Presented

### The $\Pi$-parameter table
Present as a styled table, not as equations in running text. Include a "Dominant at" column naming the body/regime where each parameter is largest. Add a small bar chart or sparkline next to each row showing relative magnitude across the four bodies. This makes the cross-body comparison instantaneous.

### The central equation + equilibrium diagram
Display the equation in a light-colored box (highlighted, not in a wall of text). Immediately below it, place the equilibrium schematic: a 2D plane with axes "drag term $\Pi_D\mathcal{F}_D$" (horizontal) and "gravitational term $\sum\Pi_i\Xi_i$" (vertical), with the equilibrium line as a diagonal and the three regions (lowering, raising, equilibrium) shaded. Label each region with its operational meaning ("periapsis drops — need raise maneuver", "stable corridor", "periapsis rises — can delay maneuver"). This is the conceptual core of the poster and should be readable in 10 seconds.

### Regime maps
Each panel should be a scatter plot or 2D heatmap in a well-chosen projection. Recommended axes:
- **Mars**: $\omega$ (x) vs. $e$ (y), fixed inclination near 93° (Odyssey-like), fixed $a$ near campaign entry.
- **Venus**: $\omega$ (x) vs. $i$ (y), fixed $a$ near Venus Express campaign entry.
- **Earth**: $\omega$ (x) vs. $\Pi_{3b}$ (proxy for apoapsis altitude) (y), representative LEO-to-GEO transfer.
- **Titan**: $\omega$ (x) vs. $e$ (y), fixed $a$ near representative capture orbit.

Overlay the analytically predicted equilibrium locus as a solid line. Agreement between the analytical prediction and numerical boundary validates Contribution 2. If disagreement exists in certain regions, note it as a known limitation (e.g., short-period terms not captured by the secular scaling).

Include a shared colorbar spanning the full $\Delta h_p$ range, with explicit contour values in km/orbit labeled on the maps.

### $\Delta h_p$ vs. $\omega$ line plot
A single figure with four curves (one per body), each normalized by its own campaign-typical $|a\Pi_D|$ so that all four bodies appear on a common vertical axis of "normalized periapsis change per orbit." This comparison immediately shows:
- Mars: smooth $\cos\omega$ sinusoidal shape, zero crossings at 90° and 270°
- Venus: $\sin(2\omega)$ shape, zero crossings at 0°, 90°, 180°, 270°
- Titan: large amplitude, potentially nonlinear shape if ZLK self-amplification is included
- Earth: mixed shape, intermediate between Mars and Venus

Mark the campaign-relevant $\omega$ range for each body as a shaded region on the x-axis.

### Cross-body comparison summary
Present as a horizontal "regime spectrum" axis or a 4×4 comparison table with qualitative regime labels (Drag-dominant / $J_2$-indirect dominant / ZLK-direct dominant / Mixed). Include one-sentence operational implications per body:
- *Mars*: "Phase I/III passively exploits $J_2$ sinusoid; equilibrium corridor lasts half the AOP cycle."
- *Venus*: "ZLK direct eccentricity change dominates; $J_2$ indirect term negligible ($f_\text{geoid} \approx 0$)."
- *Titan*: "Campaign-timescale ZLK oscillation; both channels require combined treatment."
- *Earth*: "Intermediate regime; third-body tides shift the AOP equilibrium without eliminating the $J_2$ sinusoid."

---

## 6. What Is NOT a Result (Do Not Show)

The following are framework elements or derivation steps, not results. They should not occupy poster real estate:

- Full derivations of the GVE-based $\Delta r_p$ expression (point to the reference document or technical paper in preparation)
- The Kozai Hamiltonian or Lidov-Kozai formulation in full — the poster only needs the scalar result $\langle\dot e\rangle_{3b}$ and its effect on $\Delta h_p$
- Detailed discussion of short-period oscillation handling (a one-line note in the Methods column is sufficient: "numerical propagation retains short-period terms; secular model used for regime boundaries only")
- The three-phase campaign architecture in full detail (it appears only as context in the motivation panel)
- Any RL/autonomy system architecture — this poster is about the perturbation physics, not the meta-RL framework

---

## 7. Key Takeaway (For Footer and QR Callout)

> *Gravitational perturbations are not just noise — they define natural corridors where periapsis is self-stabilizing. Four dimensionless parameters predict which regime any orbit occupies, across any body, without a full campaign simulation.*

The QR code should link to either the arXiv preprint (if available) or the SpaceAGORA.jl repository / technical reference document.

---

## 8. Questions to Resolve Before Producing the Poster

The following items must be computed or confirmed before the final figures can be generated:

1. **Numerical values of $\Pi_{J_2}$, $\Pi_{3b}$, $\Pi_H$, $\Pi_D$** for reference orbits at all four bodies — these are the numbers that fill the parameter table and calibrate the regime maps.
2. **SpaceAGORA.jl survey output** — the actual $\Delta h_p$ scatter data for the four bodies. Without this, the regime maps are schematic only.
3. **Equilibrium locus coordinates** — the analytical prediction of the $\Delta r_p = 0$ locus for each body must be computed and extracted in the same $(e, \omega)$ or $(i, \omega)$ plane used for the regime maps, to enable the overlay validation comparison.
4. **Whether short-period oscillations produce systematic offsets** between the secular scaling and the one-orbit numerical $\Delta h_p$ — this determines whether an "analytical overlay" or "analytical boundary only" is the correct claim.
5. **Titan $T_\text{ZLK}$ estimate** — confirm whether the ZLK period is campaign-timescale (months to years) or much longer (decades+). This determines whether the Titan regime map uses the secular ZLK model or requires a more careful treatment.

---

## 9. Contribution Framing for the Poster Headline

The poster should be legible with one reading of the header. The recommended framing:

**Headline (title as-is):**
*Perturbation-Assisted Aerobraking: Dimensionless Regimes and Equilibrium Periapsis Corridors Across Celestial Bodies*

**Subtitle or hook sentence (below authors, before body panels):**
*We show that four dimensionless parameters classify any aerobraking orbit as drag-dominant, gravitational-assist, or equilibrium — enabling perturbation-aware corridor design without full campaign simulation.*

This tells an expert exactly what was done, why it matters, and what was shown, before they read a single panel.