# Writes the hand-drawn documentation diagrams into docs/src/assets/figures/.
# Standard library only. Terminal, from the repository root:
#   julia --startup-file=no scripts/docs_figures/diagrams.jl
include(joinpath(@__DIR__, "svg.jl"))
using .DocsSvg

const OUT = normpath(joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "figures"))

function configuration_anatomy(path)
    s = Svg(960, 660, "What a simulation is made of",
        "A SimulationConfiguration, built by make_example_config or by a MissionPlan builder, holds the epoch, " *
        "mission settings, result settings, asset paths, the environment, the guidance, navigation and control " *
        "models, the dynamics with its force and torque models and spacecraft, and the solver settings. Each " *
        "spacecraft holds its initial condition, links joined by joints, hardware mounted on the links and " *
        "optionally its own GNC. run_simulation or run_mission turns the configuration into results.")
    text!(s, 20, 30, "What a simulation is made of"; size=15, color=INK, bold=true)

    rect!(s, 20, 62, 200, 66)
    text!(s, 120, 84, "make_example_config"; size=11.5, color=INK, bold=true, mono=true, anchor="middle")
    text!(s, 120, 101, "a complete configuration from a"; size=10.5, anchor="middle")
    text!(s, 120, 115, "planet, spacecraft, duration, epoch"; size=10.5, anchor="middle")
    rect!(s, 20, 146, 200, 80)
    text!(s, 120, 168, "MissionPlan"; size=11.5, color=INK, bold=true, mono=true, anchor="middle")
    text!(s, 120, 185, "a named builder function that"; size=10.5, anchor="middle")
    text!(s, 120, 199, "returns a fresh configuration,"; size=10.5, anchor="middle")
    text!(s, 120, 213, "plus segments, phases, events"; size=10.5, anchor="middle")
    text!(s, 120, 250, "or build the parts yourself"; size=10.5, color=MUTED, anchor="middle", italic=true)
    line!(s, 220, 95, 246, 95; arrow=true)
    line!(s, 220, 186, 246, 186; arrow=true)

    rect!(s, 250, 46, 690, 524; fill="#ffffff", stroke=ACCENT, sw=1.5, rx=4)
    text!(s, 266, 68, "SimulationConfiguration"; size=13, color=ACCENT, bold=true, mono=true)
    row1 = [("InitialTime", "the epoch", 266, 124),
            ("MissionConfiguration", "when the run ends, output cadence", 400, 206),
            ("SimulationSettings", "result files, folder, verbosity", 616, 164),
            ("FilePaths", "licensed asset locations", 790, 134)]
    for (name, what, x, w) in row1
        rect!(s, x, 80, w, 46)
        text!(s, x + w / 2, 99, name; size=11, color=INK, bold=true, mono=true, anchor="middle")
        text!(s, x + w / 2, 116, what; size=10, anchor="middle")
    end

    rect!(s, 266, 140, 380, 92)
    text!(s, 278, 160, "EnvironmentModel"; size=11.5, color=INK, bold=true, mono=true)
    for (i, l) in enumerate(["planet: gravity parameter, radii, rotation",
                             "atmosphere: density model, winds, entry interface",
                             "ephemerides: analytic, or SPICE kernels",
                             "thermal model for heat rate and heat load"])
        text!(s, 278, 178 + 14 * (i - 1), l; size=10.5)
    end
    rect!(s, 660, 140, 264, 92)
    text!(s, 672, 160, "Guidance, navigation, control"; size=11.5, color=INK, bold=true)
    rich!(s, 672, 178, [("GuidanceModel", :mono), (", ", :plain), ("NavigationModel", :mono), (",", :plain)]; size=10.5)
    rich!(s, 672, 192, [("ControlModel", :mono), (": effectors, each with", :plain)]; size=10.5)
    text!(s, 672, 206, "its call rate in seconds; control"; size=10.5)
    text!(s, 672, 220, "effectors command the actuators"; size=10.5)

    rect!(s, 266, 246, 658, 254)
    text!(s, 278, 266, "DynamicsModel"; size=11.5, color=INK, bold=true, mono=true)
    rect!(s, 278, 278, 196, 210; fill="#ffffff")
    text!(s, 290, 298, "force and torque models"; size=11, color=INK, bold=true)
    for (i, l) in enumerate(["gravity: point mass, J2,", "  spherical harmonics, n-body,", "  gravity-gradient torque",
                             "aerodynamics: drag and lift", "  coefficients, mesh panels",
                             "solar radiation pressure", "your own effector"])
        text!(s, 290, 316 + 14 * (i - 1), l; size=10.5)
    end
    rect!(s, 500, 284, 410, 204; fill="#ffffff")
    rect!(s, 494, 278, 410, 204; fill="#ffffff")
    text!(s, 506, 298, "SpacecraftModel"; size=11.5, color=INK, bold=true, mono=true)
    text!(s, 626, 298, "one per vehicle"; size=10.5, color=MUTED)
    text!(s, 506, 316, "InitialCondition: the starting orbit and attitude"; size=10.5)
    rect!(s, 506, 326, 120, 34)
    text!(s, 566, 341, "root link"; size=10.5, color=INK, bold=true, anchor="middle")
    text!(s, 566, 354, "the bus"; size=10, anchor="middle")
    for y in (384, 426)
        rect!(s, 536, y - 12, 90, 30)
        text!(s, 581, y + 6, "panel link"; size=10.5, color=INK, anchor="middle")
        path!(s, "M 518 360 L 518 $(y + 3) L 536 $(y + 3)"; color=MUTED, sw=1)
        circle!(s, 518, y + 3, 3; fill="#ffffff", stroke=MUTED, sw=1)
    end
    text!(s, 506, 470, "links joined by joints"; size=10, color=MUTED, italic=true)
    text!(s, 646, 341, "mounted on each link:"; size=10.5, color=INK, bold=true)
    for (i, l) in enumerate(["thrusters, reaction wheels,", "torque rods, laser terminals,",
                             "a hinge, SRP facets,", "your own actuator kinds"])
        text!(s, 646, 357 + 14 * (i - 1), l; size=10.5)
    end
    text!(s, 646, 425, "optional own GNC, bound to"; size=10.5)
    text!(s, 646, 439, "this vehicle when the run starts"; size=10.5)
    text!(s, 646, 455, "optional compliant attachments"; size=10.5)

    rect!(s, 266, 512, 380, 44)
    rich!(s, 278, 530, [("SolverConfig", :monobold), (" and ", :plain), ("IntegrationTolerances", :monobold)];
          size=10.5, color=INK)
    text!(s, 278, 546, "the integrator, step limits, tolerances, parallelism"; size=10.5)
    rect!(s, 660, 512, 264, 44; fill=PALE, stroke=PALE_EDGE, dash="5 4")
    text!(s, 672, 530, "optional"; size=10.5, color=INK, bold=true)
    text!(s, 672, 546, "laser interlinks, external propagators"; size=10.5)

    line!(s, 595, 570, 595, 596; arrow=true)
    rect!(s, 380, 598, 430, 46; fill=ACCENT, stroke=ACCENT, rx=3)
    rich!(s, 595, 617, [("run_simulation", :mono), ("(config)  or  ", :plain), ("run_mission", :mono), ("(plan)", :plain)];
          size=11.5, color="#ffffff", anchor="middle")
    text!(s, 595, 634, "results table, the configuration that ran, written files"; size=10.5, color="#e5eef6", anchor="middle")
    save_svg(s, path)
end

function environments_map(path)
    rows = [
        ("Install SpaceAGORA or run the", "command-line tool", "Root project", "julia --project=.",
         ("Pkg.instantiate(julia_version_strict=true)", "from the repository root")),
        ("Run the examples and the", "plotting scripts", "Examples environment", "--project=examples",
         ("examples/setup.jl", "adds Plots and the other example tools")),
        ("Run the Odyssey tutorial with a", "frozen Mars atmosphere", "Odyssey example environment",
         "--project=examples/odyssey_surrogate_env",
         ("examples/odyssey_surrogate_env/setup.jl", "adds the GRAMSuite wrapper and EDG")),
        ("Open results in the 3D viewer", "", "Viewer environment", "--project=.visualizer-env",
         ("scripts/setup_visualizer.jl .visualizer-env", "adds the SpaceAGORAVisualizer package")),
        ("Plan with HYPR, run MuJoCo", "scenes or use EDG guidance", "One environment per companion",
         "--project=/path/to/hypr-env", ("scripts/optional_packages.jl", "setup hypr|mujoco|edg /path/to/env")),
        ("Write my own scenario outside", "the checkout", "Your own project", "--project=/absolute/my-project",
         ("setup_user_project.jl, listed on", "Installation & Environment")),
    ]
    n = length(rows)
    s = Svg(960, 92 + 74 * n + 40, "Which Julia environment to start",
        "Six tasks and the Julia environment each one uses: the root project for SpaceAGORA itself and the " *
        "command-line tool; the examples environment for the examples and plotting scripts; the Odyssey example environment for the frozen-atmosphere tutorial; the viewer " *
        "environment for the 3D viewer; one environment per optional companion for HYPR, MuJoCo and EDG; and " *
        "your own project for scenarios kept outside the checkout. Every environment loads SpaceAGORA from the " *
        "same checkout.")
    text!(s, 20, 30, "Which Julia environment to start"; size=15, color=INK, bold=true)
    for (x, label) in ((20, "I want to"), (262, "Start Julia with"), (568, "Set it up once with"))
        text!(s, x, 64, label; size=11, color=MUTED, bold=true)
    end
    y0 = 76
    for (k, (want1, want2, env, flag, setup)) in enumerate(rows)
        y = y0 + 74 * (k - 1)
        rect!(s, 20, y, 228, 62; fill=PALE, stroke=PALE_EDGE)
        text!(s, 32, isempty(want2) ? y + 35 : y + 27, want1; size=11.5, color=INK)
        isempty(want2) || text!(s, 32, y + 43, want2; size=11.5, color=INK)
        line!(s, 248, y + 31, 258, y + 31; arrow=true)
        rect!(s, 262, y, 292, 62)
        text!(s, 274, y + 25, env; size=11.5, color=INK, bold=true)
        text!(s, 274, y + 44, flag; size=10.5, color=ACCENT, mono=true)
        rect!(s, 568, y, 308, 62; fill="#ffffff")
        text!(s, 580, y + 25, setup[1]; size=10.5, mono=startswith(setup[1], r"Pkg|examples|scripts"))
        text!(s, 580, y + 43, setup[2]; size=10.5, mono=startswith(setup[2], "setup "))
        line!(s, 876, y + 31, 900, y + 31; color=MUTED, sw=1)
    end
    rect!(s, 900, y0, 40, 74 * n - 12; fill=ACCENT, stroke=ACCENT, rx=3)
    cx, cy = 920, y0 + (74 * n - 12) / 2
    push!(s.parts, "<text x=\"$cx\" y=\"$cy\" font-size=\"11.5\" fill=\"#ffffff\" text-anchor=\"middle\" " *
                   "transform=\"rotate(-90 $cx $cy)\">your SpaceAGORA checkout</text>")
    fy = y0 + 74 * n + 14
    text!(s, 20, fy, "Every environment uses Julia 1.12 and loads SpaceAGORA from the same checkout. " *
                     "Only the root project ships a committed manifest;"; size=10.5)
    text!(s, 20, fy + 16, "the others resolve their own during setup. Companion environments go outside " *
                          "the source checkouts; the viewer's default folder is ignored by Git."; size=10.5)
    save_svg(s, path)
end

function mission_workflow_timeline(path)
    s = Svg(960, 500, "One mission run: segments, phases and events on one timeline",
        "A single run_mission call integrates one continuous state. Operational segments initial_coast, " *
        "observation and final_coast change commands at transitions at 120 and 180 seconds, where the exit " *
        "action runs before the entry action. Descriptive phases initial_coast and later_coast only name time " *
        "windows. An event callback is built fresh for the run and, in this sketch, stops the run before the " *
        "300 second horizon. Afterwards run.history, the saved mission JSON file and the viewer show the " *
        "segments that actually ran.")
    x(t) = 200 + 2.4 * t                      # 0 to 300 s across 720 px
    rect!(s, x(260), 44, x(300) - x(260), 290; fill="#f1f3f5", stroke="none", sw=0)
    text!(s, (x(260) + x(300)) / 2, 62, "not run"; color=MUTED, anchor="middle")
    rich!(s, 20, 28, [("One mission run: ", :bold), ("MissionPlan", :monobold), (" and ", :bold),
                      ("run_mission", :monobold)]; size=15, color=INK)
    rect!(s, x(0), 48, x(260) - x(0), 22; fill=ACCENT, stroke=ACCENT, rx=3)
    text!(s, (x(0) + x(260)) / 2, 63, "one integration, one continuous state: position, velocity, attitude, mass";
          size=12, color="#ffffff", anchor="middle")

    lane_label(y, title, typename, l1, l2) = begin
        text!(s, 20, y, title; size=13, color=INK, bold=true)
        text!(s, 20, y + 17, typename; color=ACCENT, mono=true)
        text!(s, 20, y + 33, l1)
        text!(s, 20, y + 47, l2)
    end
    lane_label(104, "Operational segments", "MissionSegment", "commands change at", "each transition")
    for (a, b, name, l1, l2, mono1) in ((0, 120, "initial_coast", "MissionDeadline(120.0)", "120 s after the start of the run", true),
                                        (120, 180, "observation", "MissionDuration(60.0)", "60 s after entry", true),
                                        (180, 260, "final_coast", "no exit condition", "closed by the stop", false))
        rect!(s, x(a), 86, x(b) - x(a), 70)
        c = (x(a) + x(b)) / 2
        text!(s, c, 111, name; size=13, color=INK, bold=true, anchor="middle")
        text!(s, c, 129, l1; size=10.5, color=mono1 ? ACCENT : SOFT, mono=mono1, anchor="middle")
        text!(s, c, 146, l2; size=10.5, anchor="middle")
    end
    for t in (120, 180)
        line!(s, x(t), 80, x(t), 162; sw=2)
        text!(s, x(t), 175, "on_exit,"; size=10.5, color=ACCENT, mono=true, anchor="middle")
        text!(s, x(t), 189, "then on_enter"; size=10.5, color=ACCENT, mono=true, anchor="middle")
    end

    lane_label(212, "Descriptive phases", "MissionPhase", "names only; the run", "is unchanged")
    for (a, b, name, win) in ((0, 120, "initial_coast", "window_s=(0.0, 120.0)"),
                              (120, 300, "later_coast", "window_s=(120.0, 300.0)"))
        rect!(s, x(a), 198, x(b) - x(a), 50; fill=PALE, stroke=PALE_EDGE, dash="5 4")
        c = (x(a) + x(b)) / 2
        text!(s, c, 218, name; size=12, color="#374151", bold=true, anchor="middle")
        text!(s, c, 236, win; size=10.5, mono=true, anchor="middle")
    end

    lane_label(292, "Events", "MissionEvent", "callback built fresh", "for each run")
    circle!(s, x(0), 302, 5; fill="#ffffff", stroke=WARM, sw=2)
    text!(s, x(0) + 12, 298, "the factory builds the callback"; size=10.5)
    text!(s, x(0) + 12, 313, "before the integration starts"; size=10.5)
    line!(s, x(260), 44, x(260), 334; color=WARM, dash="4 3")
    path!(s, "M $(x(260)) 293 L $(x(260) + 9) 302 L $(x(260)) 311 L $(x(260) - 9) 302 z"; color=WARM, fill=WARM, sw=1)
    text!(s, x(260) - 16, 298, "a stop event ends the run at the first"; size=10.5, color=WARM_INK, anchor="end")
    text!(s, x(260) - 16, 313, "accepted step after its limit"; size=10.5, color=WARM_INK, anchor="end")

    line!(s, x(0), 348, x(300), 348; color=MUTED, sw=1)
    for t in 0:60:300
        line!(s, x(t), 348, x(t), 353; color=MUTED, sw=1)
        text!(s, x(t), 366, string(t); anchor="middle")
    end
    text!(s, 20, 366, "elapsed time, s")
    text!(s, x(300), 381, "mission horizon"; size=10.5, color=MUTED, anchor="end")

    text!(s, 20, 418, "What you inspect"; size=13, color=INK, bold=true)
    text!(s, 20, 434, "after the run")
    boxes = (("run.history", true, 12, "actual segment intervals and", "transition times, each with a reason", 200, 232),
             ("simulation_results_mission.json", true, 10.5, "saved with the run and linked", "to it in the result manifest", 464, 232),
             ("Viewer", false, 12, "jump to a segment or", "transition during playback", 728, 192))
    for (k, (name, mono, sz, l1, l2, bx, bw)) in enumerate(boxes)
        rect!(s, bx, 398, bw, 66)
        c = bx + bw / 2
        text!(s, c, 419, name; size=sz, color=INK, bold=true, mono=mono, anchor="middle")
        text!(s, c, 437, l1; size=10.5, anchor="middle")
        text!(s, c, 452, l2; size=10.5, anchor="middle")
        k < 3 && line!(s, bx + bw, 431, bx + bw + 30, 431; arrow=true)
    end
    rich!(s, 200, 488, [("run.phases", :mono), (" lists the phase names; phases and events are kept separate from the segment history.", :plain)])
    save_svg(s, path)
end

mkpath(OUT)
for (f, name) in ((configuration_anatomy, "configuration_anatomy.svg"),
                  (environments_map, "environments_map.svg"),
                  (mission_workflow_timeline, "mission_workflow_timeline.svg"))
    println("wrote ", relpath(f(joinpath(OUT, name)), normpath(joinpath(@__DIR__, "..", ".."))))
end
