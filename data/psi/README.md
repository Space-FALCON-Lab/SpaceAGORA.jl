# Plume-surface interaction field tables

Gridded plume fields for the plume-surface interaction effector
(`PlumeSurfaceInteractionModel`). Each file is a `PlumeFieldTable` written by

```
julia --project=. scripts/dev/psi/build_plume_field.jl [--engine NAME] [--out PATH]
```

and read back with `load_plume_field("data/psi/<name>.json")`. They are
regenerated from the engine's own parameters, so they are reproducible rather
than measured; the point of versioning them is that a run does not have to
redo the source-flow integration.

Using one:

```julia
cfg = PlumeSurfaceConfig(field=load_plume_field("data/psi/apollo_lmde.json"))
plume = PlumeSurfaceInteractionModel(control, terrain; config=cfg)
```

Without a `field` the effector keeps its analytic Gaussian footprint, so an
existing run is unchanged unless it asks for a table.

## What is in the tables

Axes are the nozzle exit plane above the ground over the exit diameter
(`h/D_e`, 24 logarithmic nodes from 0.5 to 200) and the ground radius over that
height (`r/h`, 64 nodes uniform in incidence angle out to 76 degrees). At each
node the table holds the surface pressure, the wall shear stress, and the
density, speed, temperature and Mach number of the wall jet, all nondimensional
and carrying the thrust scaling explicitly (see `PlumeFieldTable`), so one table
serves the whole throttle band. A query outside either axis clamps.

## The model

The plume is the Simons source-flow model of a rocket exhausting into vacuum
with Boynton's nozzle-boundary-layer correction:

- G. A. Simons, "Effect of nozzle boundary layers on rocket exhaust plumes",
  *AIAA Journal* 10(11), 1972, pp. 1534-1535.
- F. P. Boynton, "Exhaust plumes from nozzles with wall boundary layers",
  *J. Spacecraft and Rockets* 5(10), 1968, pp. 1143-1147.
- The explicit equations implemented are section 2.1, equations (3) to (14), of
  P. J. Herraiz, J. M. Fernandez and J. R. Villa, "Development of a MATLAB plume
  impingement tool for fast system analysis", 8th European Conference for
  Aeronautics and Aerospace Sciences (EUCASS), 2019, DOI 10.13009/EUCASS2019-661.

The plume meets the ground through classical Newtonian impact theory for the
wall pressure and the oblique normal-shock relations for the state behind the
surface shock; the wall shear stress is the wall jet's dynamic pressure times
the rough-surface drag coefficient of 0.2 that L. Roberts used for the
meteorite-gardened lunar surface (Roberts 1966, quoted by Morris 2012 section
4.4.2, below).

Classical Newtonian is used rather than modified Newtonian because it satisfies
the axial-momentum balance on the ground exactly. The built tables recover it:
the integral of the surface pressure over the ground comes to 1.008 times the
engine thrust at every height, the residual being the difference between the
source flow's limiting velocity and the engine's vacuum thrust.

## Engine parameters and where they come from

`PlumeNozzle` carries them all, with the source or the derivation on every
field. In short, for the Apollo lunar module descent engine (LMDE):

| quantity | value | status |
|---|---|---|
| specific impulse, vacuum | 305 s | sourced (Morris 2012 Table 2.2.1; TRW data via "Rocket Propulsion Evolution" 9.42) |
| area ratio | 47.5 | sourced (same two) |
| exit Mach number | 5.03 | sourced (Morris 2012 Table 3.1.1) |
| chamber temperature | 2730 K | sourced (Morris 2012 Table 3.1.1) |
| exit diameter | 1.50 m | sourced (59 in nozzle of the area-ratio 47.5 engine Apollo 11 flew) |
| surface drag coefficient | 0.2 | sourced (Roberts 1966, via Morris 2012 section 4.4.2) |
| ratio of specific heats | 1.30 | derived from Morris 2012 Table 3.1.1 |
| gas constant | 461 J/(kg K), 18.0 kg/kmol | derived from Morris 2012 Table 3.1.1 |
| limiting exhaust speed | 3303 m/s | derived |
| divergent length | 1.91 m | derived (80 percent bell) |
| exit boundary-layer thickness | 0.0595 exit radii | derived from the exit Reynolds number of Morris 2012 Table 3.1.1 |
| nozzle exit lip half-angle | 10 degrees | **assumption** |
| mean limit-speed ratio | 0.75 | **assumption** (the source bounds it only between 0.5 and 1) |

A consistency check that is not an input: the chamber pressure the model
implies from the choked throat is 6.81 bar, against the 103.4 psia (7.13 bar)
the TRW data sheet gives for full thrust, 4.5 percent low.

Morris 2012 is A. B. Morris, *Simulation of rocket plume impingement and dust
dispersal on the lunar surface*, PhD dissertation, University of Texas at
Austin, 2012.

## The two shipped tables

Both are the same engine and differ only in where the point source sits.

| file | source position | peak wall shear at 5 m, 13.34 kN | near the ground |
|---|---|---|---|
| `apollo_lmde.json` | nozzle exit plane | 41.4 Pa | well behaved to contact |
| `apollo_lmde_virtual_source.json` | 1.8 m downstream of the exit plane | 100.0 Pa | freezes below about 2.6 m |

The reference for that column is the only published number found for this engine
at a stated condition: Morris 2012 section 4.4.2 reports a peak laminar
smooth-wall shear stress of 92 Pa for the LMDE hovering 5 m above the surface at
13.34 kN, from a DSMC simulation. The same work (section 4.6) measures the
plume's effective point source about 1.8 m below the exit plane, where the
nozzle's internal compression wave meets the axis, and identifies a source at
the exit plane as the reason Roberts' theory under-predicts the surface stress
below about ten nozzle diameters.

So the default table is the published source-flow model taken literally and is
2.2 times below the one DSMC number available; moving the source to the offset
that work measured brings it to within 9 percent of it, at the cost of a field
that has to be frozen within a few meters of contact. The default ships the
conservative one. Neither was tuned to the 92 Pa: the offset is that paper's own
measurement, used as published.

## Against the analytic field

With the Apollo 11 approach thrust (11.5 kN) and `apollo_lmde.json`:

| height | stagnation pressure, analytic / table | footprint radius | peak wall shear |
|---|---|---|---|
| 5 m | 673 / 821 Pa | 2.33 / 2.20 m | 5.78 / 35.6 Pa |
| 10 m | 168 / 205 Pa | 4.66 / 4.40 m | 1.44 / 9.18 Pa |
| 31 m | 17.5 / 21.4 Pa | 14.5 / 13.6 m | 0.150 / 0.984 Pa |

The two agree on the pressure to 22 percent and on the footprint radius to
6 percent over the whole descent, which is the real check that the Gaussian
footprint was a reasonable stand-in. They do not agree on the shear stress: the
table reports six to seven times more, because Roberts' rough-surface drag
coefficient acting on the wall jet's dynamic pressure is a much stronger
coupling than the 0.01 skin-friction coefficient acting on the static pressure.
That moves the erosion onset height for this thrust from 31 m to 80 m at the
unchanged 0.15 Pa threshold. The threshold is a calibrated parameter of the
erosion closure, not of the plume, and it is not retuned here. Re-evaluating the
Apollo 11 descent on the tabulated field with the erosion constants left where
they are moves 29 tonnes of regolith against the analytic field's 429 kg, far
above the tonne-scale estimates from the landings: the calibration does not
transfer between the two fields, and closing that is the erosion model's job,
not the plume field's. See `docs/src/user/lunar_landing.md`.

## Against the Apollo measurements

Lane and Metzger, "Estimation of Apollo lunar dust transport using optical
extinction measurements", *Acta Geophysica* 63(2), 2015, 568-599, inverted the
Apollo 12 descent films for the plume's action on the ground. Two of their
results bear directly on this field: the radius of the eroding region runs
about 1.3 to 2.2 times the nozzle height (their Table 2, equation 11), an
effective half-angle of 52 to 66 degrees; and their fitted surface shear stress
is `tau(h) = 6.21 exp(-0.123 h)` Pa (their equation 14 and Figure 11). Both
results are taken here from the validation study's reading of that paper
(`benchmarks/studies/psi_validation/`); they have not been checked against the
primary source from this side, and the comparisons below inherit that.

At the 11.5 kN approach throttle, with `apollo_lmde.json`:

| height | eroding radius / height | pressure footprint / height | Lane and Metzger |
|---|---|---|---|
| 5 m | 1.92 (62 deg) | 0.44 (24 deg) | 1.3 to 2.2 |
| 10 m | 1.48 (56 deg) | 0.44 (24 deg) | 1.3 to 2.2 |
| 15 m | 1.27 (52 deg) | 0.44 (24 deg) | 1.3 to 2.2 |
| 31.5 m | 0.95 (44 deg) | 0.44 (24 deg) | 1.3 to 2.2 |

The computed field reproduces the measured width of the eroding region at the
heights the films cover, and it does so without widening the plume: the
momentum-carrying pressure core stays at about 24 degrees, the same width the
analytic field assumes, while the region that can move soil is two to four
times wider. The two decouple because the wall shear stress is the skin
friction of the radial wall jet, `C_D` times half `rho u^2`, and the wall jet
keeps accelerating outward (165 m/s at 1 m from the stagnation point at a 10 m
height, 2063 m/s at 8 m, 3133 m/s at 30 m) while its density thins, so the
shear peaks well outside the pressure footprint and falls off far more slowly
than the pressure does. A shear stress proportional to the local static
pressure, which is what the analytic field uses, cannot do that: it forces the
eroding region and the pressure footprint to have the same width.

The shear stress itself:

| height | Lane and Metzger fit | model peak | model mean inside 1.3 h | model mean inside 2.2 h |
|---|---|---|---|---|
| 5 m | 3.36 Pa | 35.7 Pa | 12.3 Pa | 4.51 Pa |
| 10 m | 1.82 Pa | 9.19 Pa | 3.16 Pa | 1.16 Pa |
| 15 m | 0.98 Pa | 4.14 Pa | 1.42 Pa | 0.52 Pa |
| 31.5 m | 0.13 Pa | 0.95 Pa | 0.33 Pa | 0.12 Pa |

Read that with care, because the comparison is only meaningful if the two sides
are the same quantity. Against the model's *peak* the fit is 5 to 11 times low.
Against the model averaged over the eroding region it is far closer: the average
inside 1.3 h is 1.5 to 3.7 times the fit and the average inside 2.2 h is 0.5 to
1.3 times it, so the fit lies between the two averages at 10, 15 and 31.5 m and
1.3 times below even the wider average at 5 m -- with no constant fitted
anywhere in the field. Which comparison is like-for-like depends on whether
equation 14 of that paper is a peak or an area average, which has not been
checked here. The analytic field averaged
the same way over 1.3 h gives 0.38, 0.17 and 0.039 Pa at 10, 15 and 31.5 m,
3 to 6 times *below* the fit, so on an area-averaged basis the computed field is
much the closer of the two and the analytic field's apparent agreement comes
from comparing its peak against an average.

The model's shear decays as `h^-2`, the fit as `exp(-0.123 h)`. Neither shape
matches the other, and the DSMC work of Morris 2012 section 4.6 measures a
third, `h^-2.825`. That disagreement is unresolved here.

The two published references are also not consistent with each other: at 5 m
Morris 2012's DSMC gives 92 Pa (at 13.34 kN) and the Lane and Metzger fit gives
3.4 Pa (at 11.5 kN), a factor of 25 apart for nearly the same condition. The
model's peak, 41 Pa scaled to 13.34 kN, falls between them. They are different
quantities -- a computed smooth-wall boundary-layer stress against a stress
inverted from a dust cloud through an erosion model -- so the gap is not
necessarily an error in either, but it does mean no single number here can be
treated as ground truth.

### Erosion onset

The height at which the peak shear first clears the soil's threshold, at
11.5 kN, against the roughly 31.5 m at which the Apollo crews reported blowing
dust. The 0.057 Pa is the Shields-Bagnold threshold this round's regolith work
derived from the Lunar Sourcebook soil through Shao and Lu, *J. Geophys. Res.*
105(D17), 2000, equation 22; it is quoted here from that work, not rederived.

| threshold | analytic field | computed field |
|---|---|---|
| 0.15 Pa (the fitted value) | 31.0 m | 80.0 m |
| 0.057 Pa (from the soil, not fitted) | 50.3 m | 130.0 m |

The analytic field reproduces the 31.5 m only because 0.15 Pa was fitted to make
it do so. The computed field with a soil-derived threshold puts the onset 2.5 to
4 times higher, and that residual is not resolved here: the reference is the
height at which dust became *visible*, which is an obscuration threshold rather
than a measurement of when grains first move, and the Lane and Metzger fit
itself crosses 0.15 Pa at 30.3 m, so their inferred shear and this model's
differ by the same 5 to 11 times seen above. Closing it needs either the primary
source for equation 14 or a wall-jet boundary-layer treatment finer than a
constant drag coefficient. No constant was retuned to hide it.

### For the erosion models built on this field

An energy-flux erosion law needs a density and a temperature that are not
restatements of the pressure. The tabulated field's are not: the density comes
from the plume's own mass flux per steradian divided by the limiting speed and
the square of the slant range, and the temperature from the surface shock, so
between 1 m and 20 m from the stagnation point at a 10 m height the pressure
falls by a factor 18,000 while the density falls by 3,500 and the temperature by
4.3. `PlumeAnalyticField` cannot offer that -- with no thermodynamics of its own
it has to reconstruct density and temperature from its own pressure -- so an
energy-flux law should be evaluated on a table.

One caution that is physics rather than a modeling artifact: for a point source
over an infinite plane the integral of the kinetic energy flux over the ground
is independent of height, exactly as the integral of the pressure is. A law
that integrates the flux with no threshold will therefore be altitude-blind
whatever field it is given; the altitude dependence lives in the area over which
the flux clears the threshold.
