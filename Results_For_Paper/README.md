# Paper figures

For the **closed-form overlays**, run `julia Results_For_Paper/regenerate/regenerate.jl`
from the project root. To redraw from the existing prediction tables, run
`julia Results_For_Paper/regenerate/plot_overlays.jl`. Both use the original Julia
paper renderer for fonts, axes, subplot spacing and legends. The
[regenerate directory](regenerate/) contains the figures, [captions](regenerate/captions.md),
the raw prediction grid and the predictions used in the plots. The overlay legends
separate simulation markers (circles for a-max, triangles for i-max) from prediction
lines (solid for a-max, dotted for i-max). Each section puts Δa-max on the first row
and Δi-max on the second, with Δi = 0.0°, Δi = 0.5°, and Δi = 1.0° across each row.
A thin black frame and dividers separate the legend sections and columns.
The overlay legend is centered on the full image by default. To nudge both
legends left, change `legend_x_offset=0.0` to `-0.01` in
[plot_overlays.jl](regenerate/plot_overlays.jl) and redraw. This moves the whole
box by 1% of the footer width without moving the data panels; positive moves right.
Section and rule labels are centered in their cells. Markers are centered over the
line samples, and each Δi label aligns down its column. The footer adds space above
and below each row without reducing the font size. The OE legend labels the
dashed continuous-access lines in a merged left cell as
**Threshold,** $N_{\mathrm{cont}}$.
At 1000 km, a-max below 232 helpers uses `leapfrog`: the target starts at +1 degree
and the trailing helper at 0 degrees throughout this sparse grid. The `leapfrog_rtn`
rows use a numerical extension of the scalar leapfrog formula. At each step it resolves
the trailing helper's along-track gap, checks the 3D range, and projects the push onto
the local radial, tangential and normal axes. The two normal branches are integrated
separately. Tangential impulse and semi-major-axis gain are recomputed consistently;
inclination and node changes are not predicted by this fixed-plane model.
The new lines use no simulated impulses or fitted values. The node comparison
subtracts the matching Output 4 links-off endpoint from each simulated endpoint.

The original leapfrog [Δv figure](../Output_6_Close_Form_Solution/dv_RTN_original_leapfrog_100orbits.png)
and [orbital-element figure](../Output_6_Close_Form_Solution/OE_original_leapfrog_100orbits.png)
use the same 100-orbit J2-on simulation markers and the same prediction grid. For
1000 km Δa-max cases below 232 helpers, the original formula predicts only tangential
Δv and semi-major-axis change (Δa). Gaps in the radial, signed-normal, inclination,
and node lines mean **not predicted**, not zero. Eccentricity shows simulation markers only.
Run `julia Output_6_Close_Form_Solution/plot_original_leapfrog.jl` from the project
root to recreate both original-model figures as PNG and PDF without changing the extended-model figures.

The [combined 100-orbit CSV](combined_paper_data_100orbits.csv) contains one row
per case for both figures. Regenerate it from this directory with
`python3 combine_paper_data.py`. Inclination is the initial target inclination;
the negative normal impulse retains its sign, and blank Delta Omega means
RAAN is undefined for the equatorial orbit.

The renderer also accepts optional `policies`, `policy_labels`, `specification`,
`orbit_window`, `endpoint_row` and `source_limits` keywords for the Output 3/4/5
campaigns. Defaults retain this folder's original data, styles and rendering.
Alternate campaigns supply their own complete grids, correctly named legends
and fitted limits; they do not overwrite the figures or CSVs here.

[plot_paper_results.jl](plot_paper_results.jl) contains a local adaptation of the
original grid renderer. It recreates only the selected 100-orbit OE and split-N
delta-v figures from the three copied CSVs here, without importing the original
plotting code or rereading the simulation trajectories.

Run from the project root:

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Results_For_Paper/plot_paper_results.jl
```

Edit the settings at the top:

- `LINE_COLORS`: six original colors, ordered by inclination 0, 0.5, 1 degree.
- `LINE_STYLES`: original solid SMA / dotted INC; accepts `:solid`, `:dash`,
  `:dot`, `:dashdot`, etc. The legend follows these settings.
- `OE_SIZE` / `DV_SIZE`: subplot region width/height and legend footer height.
  With `figure_size=nothing`, total width is three subplot widths and total
  height is four subplot heights plus the footer. An explicit
  `figure_size=(width,height)` overrides the canvas and rescales the regions;
  it does not preserve their previous absolute sizes. Regions include axes and
  margins, not just the inner data rectangles. PNG resolution also scales with DPI.
- `FONT_FAMILY`, `FONT_SIZE`, `TICK_FONT_SIZE`, `GUIDE_FONT_SIZE`: Computer
  Modern, 7-point ticks, an 8-point legend and 9-point axis labels by default.
- `LINE_WIDTH`, `MARKER_SIZE`, `ACCESS_LINE_WIDTH`: thin curves, small markers
  and thinner dashed access thresholds.
- `X_TICK_ROTATION`: 0 degrees for horizontal helper-count tick labels.
- `LABELED_HELPERS`: 1, 200, 400, 600, 800, to avoid crowded tick labels.
  All nine helper populations are still plotted, with ticks/grid lines at each.
- `OE_SPACING` / `DV_SPACING`: independent settings for the two figures:
  - `panel_margins_mm`: left/right/top/bottom spacing around each panel,
    without later figure-level overrides.
  - `y_label_x_mm`: horizontal position of all four y-axis titles, measured
    from the left edge of that figure. Increase it to move the titles right.
    Titles remain aligned and vertically centered on their data panels.
  - `y_label_margin_mm`: extra space reserved beside the first column for titles.
  - `bottom_row_margin_mm` / `legend_top_margin_mm`: extra spacing below the
    bottom panels and above the footer. Defaults of -1 mm partly offset GR's
    built-in 2 mm padding while retaining tick/axis-label clearance.
  Editing one block does not affect the other.
  Programmatic calls can override a figure's block with `make_plot(:oe; spacing=...)`.
- `OE_LEGEND` / `DV_LEGEND`: independent legend placement and row spacing:
  - `x_offset`: horizontal shift as a fraction of the footer plotting width.
    Use `-0.03` to move all legend samples and labels left by 3% of that width.
  - `y_offset`: vertical shift as a fraction of the footer plotting height;
    positive moves up and negative moves down.
  - `row_spacing`: multiplier on the existing row gaps, anchored at the first
    row. `0.8` tightens them by 20%; `1.2` increases them by 20%.
  Defaults `(x_offset=0.0, y_offset=0.0, row_spacing=1.0)` preserve placement.
  These controls also move the OE access-threshold row. They move legend
  contents, unlike `legend_top_margin_mm`, which changes layout padding.
  Keep offsets modest to avoid clipping or overlap; increasing `footer_height`
  allocates more room for widely spaced rows. Changes to footer size also
  change the physical distance represented by a fractional offset.
  Calls may supply `make_plot(:oe; legend_settings=...)`.
- `DPI`: 300 for the PNG exports; PDF curves and text remain vector.
- `SHOW_OE_ACCESS=true` / `SHOW_DV_ACCESS=false`: independently control the
  threshold lines and their legend. The DV figure has neither by default.
- `Y_LIMITS`: original full-campaign shared ranges in source units (m, degrees,
  etc.), except RAAN. RAAN limits are automatically fitted to all non-equatorial
  100-orbit points in the local paper CSV (1000, 1050, 1100 km; both schedulers),
  then shared by its three panels. No 80-orbit data enter those limits.
  Rendering fails rather than silently clipping out-of-range values.
- `SPECS.oe.scales` / `labels`: display conversions and their matching labels.
  Semimajor-axis change is in km (source m / 1000), eccentricity change in units
  of 10^-4 (source value x 10000), inclination change in units of 10^-3 deg
  (source degrees x 1000), and RAAN change in degrees. CSVs retain their original
  units. Inclination tick values are unchanged from the former mdeg display.
  Label strings can be edited independently to adjust LaTeX spacing.

There is no overall title or descriptive footer caption. Each column has a
single header above its top panel: `Target Alt.: 1000 km`, etc., using
`FONT_SIZE`. The bottom x-axis identifies only the helper count. Columns remain
ordered 1000, 1050, 1100 km. All figures represent 100 initial target periods.

The compact reference-style layout retains four quantity rows and three altitude
columns, with a two-row, three-column legend (`Δa-max, 0.5°`, `Δi-max, 0.5°`, etc.).
Axes use typeset mathematical symbols/subscripts and plain numeric y-tick
labels to avoid repeated long scientific-notation labels narrowing the panels.
OE has an additional three-entry access-threshold legend row. Default canvases are
510 x 586 (OE) and 510 x 572 (DV), approximately 180 x 207 mm and 180 x 202 mm
in the exported PDFs. Inspect the PDF at its intended placement size before
submission; resizing it in a paper also resizes its fonts.

OE changes include J2 and laser effects; RAAN is unwrapped and equatorial cases
are omitted. DV rows show signed R, signed T, positive N and negative N,
in the instantaneous RTN frame. Negative contributions retain their signs.
The OE dashed thresholds remain a fixed-circular, range-only estimate.
Their colors identify inclination using the SMA palette; the estimate is the
same for both schedulers at a given initial geometry.
Each inclination (0, 0.5, 1 degree) has its own matching dashed legend sample.
`N_cont = ceil(pi / alpha_min)` is the estimated continuous-access helper count;
its values are stored in the threshold CSV and vary with altitude.

Data copies:
- [orbital_element_deltas_100orbits.csv](orbital_element_deltas_100orbits.csv):
  from OE folder `4_combined_1000_1050_1100km`.
- [cumulative_dv_RT_N_branches_100orbits.csv](cumulative_dv_RT_N_branches_100orbits.csv):
  from DV folder `6_signed_RT_positive_negative_N_1000_1050_1100km`.
- [continuous_access_thresholds.csv](continuous_access_thresholds.csv): from OE folder 4.

The original figures, plotting scripts and data are not modified.
