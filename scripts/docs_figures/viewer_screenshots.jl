# Writes the interactive-viewer screenshots into docs/src/assets/figures/.
# Needs the viewer environment (scripts/setup_visualizer.jl), Google Chrome and Node.js 22 or
# newer. Terminal, from the repository root:
#   julia --project=.visualizer-env scripts/docs_figures/viewer_screenshots.jl
using SpaceAGORA
using SpaceAGORAVisualizer
using CSV, DataFrames, JSON
import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "docs", "src", "assets", "figures")
const WORK = mktempdir()
mkpath(OUT)

# 1. The page from "Make your first page" on Interactive Visualization.
earth = joinpath(WORK, "viewer-earth")
SpaceAGORA.run_cli(["run", "--example=AGORA_Earth_NoGRAM.jl", "--visualize", "--output-dir=$earth"])

# 2. One Mars aerobraking pass with drag, drawn with the NASA Mars Odyssey model.
mars = joinpath(WORK, "viewer-mars-pass")
planet = make_no_gram_planet(:mars)
spacecraft = make_three_body_spacecraft(
    bus_dims=(2.2, 1.7, 2.6), panel_dims=(0.05, 2.0, 1.8),
    bus_mass=600.0, panel_mass_each=15.0, panel_offset_y=2.1,
    ic=InitialCondition(ra=planet.Rp_e + 400e3, rp=planet.Rp_e + 100e3,
                        i=93.0, ω=40.0, Ω=10.0, ν=200.0),
    prop_mass=60.0, id=1)
args = make_example_config(
    planet=planet, spacecraft=spacecraft, mission_time=5400.0,
    initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
    density_model=ExponentialAtmosphereModel(planet),
    ephemerides_model=SimpleEphemeridesModel(),
    orientation_sim=false, keplerian=true, EI_km=140.0, verbose=false,
    results_directory=mars)
run_simulation(with_visualization_scene(args, true))
export_visualization(joinpath(mars, "simulation_results");
    models=Dict(1 => joinpath(ROOT, "data", "scene", "models", "mars_odyssey_nasa_3d_resources.glb")),
    model_scale=1.0)
table = CSV.read(joinpath(mars, "simulation_results.csv"), DataFrame)
periapsis_fraction = table.time[argmin(table.sc1_altitude)] / table.time[end]

# 3. Screenshots. Controls are found by their data-role attributes in the viewer page.
role(js) = "(() => { const r = (n) => document.querySelector(`[data-role=\"\${n}\"]`); $js })()"
scrub(f) = role("const p = r('play'); if (p.textContent.trim() === 'Pause') p.click(); " *
                "const s = r('scrub'); s.value = $f; s.dispatchEvent(new Event('input', {bubbles: true}));")
shots = [
    Dict("url" => "file://" * joinpath(earth, "simulation_results_viewer.html"), "width" => 1400, "height" => 1000,
         "steps" => [Dict("wait" => 8000), Dict("eval" => scrub(0.85)),
                     Dict("eval" => role("const c = r('trailcolor'); c.value = 'altitude'; " *
                                         "c.dispatchEvent(new Event('change', {bubbles: true}));")),
                     Dict("wait" => 1000), Dict("drag" => [700, 300, 700, 520]), Dict("wait" => 1000),
                     Dict("wheel" => [700, 480, -120, 3]), Dict("wait" => 1500),
                     Dict("shot" => joinpath(OUT, "viewer_overview.jpg"))]),
    Dict("url" => "file://" * joinpath(mars, "simulation_results_viewer.html"), "width" => 1400, "height" => 1000,
         "steps" => [Dict("wait" => 9000), Dict("eval" => scrub(round(periapsis_fraction; digits=4))),
                     Dict("eval" => role("const l = r('labels'); if (l.checked) l.click(); r('follow').click();")),
                     Dict("wait" => 4500), Dict("wheel" => [700, 450, -120, 5]), Dict("wait" => 1000),
                     Dict("drag" => [900, 350, 820, 300]), Dict("wait" => 1500),
                     Dict("shot" => joinpath(OUT, "viewer_heating_closeup.jpg")),
                     Dict("click" => [775, 470]), Dict("wait" => 1500),
                     Dict("shot" => joinpath(OUT, "viewer_face_inspection.jpg"))]),
]
for (k, shot) in enumerate(shots)
    steps = joinpath(WORK, "steps_$k.json")
    write(steps, JSON.json(shot))
    run(`node $(joinpath(@__DIR__, "capture_viewer.mjs")) $steps`)
end
