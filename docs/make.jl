# Build from the docs environment, which takes SpaceAGORA from this checkout:
# `julia --startup-file=no docs/setup.jl`, then `julia --project=docs docs/make.jl`.
# The setup writes the environment's ignored manifest; a fresh checkout has none.
const DOCS_SETUP = "run `julia --startup-file=no docs/setup.jl`, then `julia --project=docs docs/make.jl`"
samefile(something(Base.active_project(), ""), joinpath(@__DIR__, "Project.toml")) ||
    error("Build the documentation from the docs environment: $(DOCS_SETUP).")
isfile(joinpath(@__DIR__, "Manifest.toml")) ||
    error("The docs environment is not prepared: $(DOCS_SETUP).")
using Documenter
using SpaceAGORA
using Dates
samefile(pkgdir(SpaceAGORA), dirname(@__DIR__)) ||
    error("The docs environment loaded SpaceAGORA from $(pkgdir(SpaceAGORA)), not from this checkout; rerun `julia --startup-file=no docs/setup.jl`.")

include(joinpath(@__DIR__, "public_api_symbols.jl"))
using .PublicAPISymbols
include(joinpath(@__DIR__, "extension_hooks.jl"))
using .ExtensionHooks

const REPO_URL = "https://github.com/Space-FALCON-Lab/SpaceAGORA.jl"
const DOCS_SRC = joinpath(@__DIR__, "src")
const GENERATED_SRC = joinpath(DOCS_SRC, "generated")
const GENERATED_API_PAGE = joinpath(GENERATED_SRC, "public_api.md")
const GENERATED_HOOK_REFERENCE_PAGE = joinpath(GENERATED_SRC, "extension_hook_reference.md")
const API_REFERENCE_PAGES = (
    (section="Simulation Setup", title="Simulation Setup API", file="simulation_setup_api.md"),
    (section="RPO planner pilot", title="RPO Planner API", file="rpo_planner_api.md"),
    (section="Burn-plan family", title="Burn-Plan Family API", file="burn_plan_family_api.md"),
    (section="Joint-reference family", title="Joint-Reference Family API", file="joint_reference_family_api.md"),
    (section="EDG guidance family", title="EDG Guidance Family API", file="edg_guidance_family_api.md"),
    (section="GNC algorithm contract", title="GNC Algorithm Contract API", file="gnc_algorithm_contract_api.md"),
    (section="GNC runtime", title="GNC Runtime API", file="gnc_runtime_api.md"),
    (section="Robotics and Planning", title="Robotics API", file="robotics_api.md"),
    (section="Compliant Multibody Dynamics", title="Compliant Multibody API", file="compliant_multibody_api.md"),
    (section="Terrain Queries", title="Terrain API", file="terrain_api.md"),
    (section="Mesh Aerodynamics", title="Mesh Aerodynamics API", file="mesh_aerodynamics_api.md"),
    (section="Visualization Scene", title="Visualization API", file="visualization_api.md"),
    (section="Lunar Landing", title="Lunar Landing API", file="lunar_landing_api.md"),
    (section="Actuators", title="Actuators API", file="actuators_api.md"),
    (section="Sensors", title="Sensors API", file="sensors_api.md"),
    (section="Model Extension Hooks", title="Extension Hooks API", file="extension_hooks_api.md"),
    (section="Contributor conformance", title="Contributor Conformance API", file="contributor_conformance_api.md"),
    (section="Interlinks", title="Interlinks API", file="interlinks_api.md"),
    (section="Constellations", title="Constellations API", file="constellations_api.md"),
)
const BUILD_DIR = joinpath(@__DIR__, "build")
const PUBLIC_API_REPORT = joinpath(BUILD_DIR, "undocumented_public_exports.txt")

function _write_generated_public_api_page()::String
    mkpath(GENERATED_SRC)
    separated_titles = [page.section for page in API_REFERENCE_PAGES]
    main_sections = filter(section -> section.title ∉ separated_titles, PUBLIC_API_SECTIONS)
    open(GENERATED_API_PAGE, "w") do io
        write(io, render_public_api_markdown(SpaceAGORA; sections=main_sections))
        for page in API_REFERENCE_PAGES
            println(io, "See the [$(page.title)]($(page.file)) for $(lowercase(page.section)).")
            println(io)
        end
    end
    for page in API_REFERENCE_PAGES
        sections = filter(section -> section.title == page.section, PUBLIC_API_SECTIONS)
        length(sections) == 1 || error("Expected one $(page.section) API section")
        open(joinpath(GENERATED_SRC, page.file), "w") do io
            write(io, render_public_api_markdown(SpaceAGORA; sections, title=page.title))
            println(io, "See the [Public API](public_api.md) for the rest of the exported interface.")
        end
    end
    return GENERATED_API_PAGE
end

_table_cell(text::AbstractString) = replace(text, "|" => "\\|", "\n" => " ")

function _page_title(page::AbstractString)::String
    for line in eachline(joinpath(DOCS_SRC, page))
        startswith(line, "# ") && return strip(line[3:end])
    end
    error("docs/src/$(page) has no title")
end

_page_link(page) = page === nothing ? "not yet" : "[$(_page_title(page))](../$(page))"

# The Extension Hook Reference page, rendered from docs/extension_hooks.jl.
function _write_extension_hook_reference_page()::String
    mkpath(GENERATED_SRC)
    open(GENERATED_HOOK_REFERENCE_PAGE, "w") do io
        print(io, """
            # Extension Hook Reference

            Every extension hook SpaceAGORA exports, by family: the function a model,
            effector, actuator, sensor or algorithm implements for its own type, when the
            engine calls it, and what happens for a type with no method of its own. Each
            hook links to its docstring in the API reference, and the last column to the
            page that explains it. [Extensibility](../extensibility.md) shows how the
            hooks fit together.

            This page is generated from `docs/extension_hooks.jl`, which
            `test/gates/ci_extension_hook_reference_gate.jl` checks against the package:
            every exported hook has a row, each signature is the first line of the hook's
            docstring, and each fallback agrees with the hook's methods.

            The **shape** is the kind of hook; the Called column says when each one runs.

            | Shape | Kind |
            |---|---|
            | T, tick | updates a GNC effector's held state, periodically, for one spacecraft |
            | R, right-hand side | a term of the equations of motion, for one spacecraft |
            | Q, query | answers an environment query |
            | D, declaration | declares a property of a model |
            | A, accessor | reads a held value for one spacecraft |
            | H, hardware | one actuator's or sensor's physics |
            | P, planner or algorithm | proposes, checks or installs a planning product |
            | L, configuration | prepares a model or an algorithm for a run |

            The **fallback** is what a type with no method of its own gets: `none` means
            every type that uses the hook needs its own method; otherwise the fallback
            returns the value given or throws. The **status** is `current`, or deprecated
            with its replacement and the release that removes it.

            """)
        for group in EXTENSION_HOOK_GROUPS
            println(io, "## ", group.title, group.core_owned ? " (core-owned)" : "")
            println(io)
            if group.note !== nothing
                println(io, group.note)
                println(io)
            end
            println(io, "| Hook | Shape | Signature | Called | Arguments and units | Fallback | Status | Explained on |")
            println(io, "|---|---|---|---|---|---|---|---|")
            for row in group.rows
                status = row.status == "current" ? "current" :
                    "deprecated → [`SpaceAGORA.$(row.replacement)`](@ref), removed in $(row.removed_in)"
                cells = ("[`SpaceAGORA.$(row.symbol)`](@ref)", row.shape, "`$(row.signature)`", row.called,
                         row.arguments, row.fallback, status, _page_link(row.page))
                println(io, "| ", join(_table_cell.(cells), " | "), " |")
            end
            println(io)
        end
        println(io, "## Versioned facades")
        println(io)
        println(io, "An optional companion package reaches SpaceAGORA through one versioned facade. ",
                    "Its page lists the functions the facade provides and consumes.")
        println(io)
        println(io, "| Facade | Contract | Explained on |")
        println(io, "|---|---|---|")
        for facade in VERSIONED_FACADES
            version = getproperty(SpaceAGORA, facade.symbol).CONTRACT_VERSION
            cells = ("[`SpaceAGORA.$(facade.symbol)`](@ref)", string(version), _page_link(facade.page))
            println(io, "| ", join(_table_cell.(cells), " | "), " |")
        end
    end
    return GENERATED_HOOK_REFERENCE_PAGE
end

function _remove_retired_contract_mirrors()
    # These files were generated by older builds, never contributor-owned sources.
    retired_contracts = joinpath(GENERATED_SRC, "contracts")
    isdir(retired_contracts) && rm(retired_contracts; recursive=true)
    return nothing
end

function _write_undocumented_public_api_report()::Vector
    mkpath(BUILD_DIR)
    missing = undocumented_public_api_specs(SpaceAGORA)
    open(PUBLIC_API_REPORT, "w") do io
        println(io, "SpaceAGORA undocumented public exports report")
        println(io, "generated_utc=$(Dates.now(Dates.UTC))")
        println(io)
        println(io, "Known placeholder docs (listed in docs/public_api_symbols.jl): ",
            join(sort!([string(owner, ".", symbol) for (owner, symbol) in KNOWN_PLACEHOLDER_DOCS]), ", "))
        println(io)
        if isempty(missing)
            println(io, "All documented.")
        else
            for spec in missing
                println(io, "$(spec.section): $(nameof(spec.owner_module)).$(spec.symbol)")
            end
        end
    end
    if isempty(missing)
        @info "All generated public API symbols have docstrings."
    else
        for spec in missing
            @warn "Undocumented public export" section=spec.section owner_module_name=String(nameof(spec.owner_module)) symbol=String(spec.symbol)
        end
    end
    return missing
end

_write_generated_public_api_page()
_write_extension_hook_reference_page()
_remove_retired_contract_mirrors()
missing_public_api = undocumented_public_api_specs(SpaceAGORA)

DocMeta.setdocmeta!(SpaceAGORA, :DocTestSetup, quote
    using SpaceAGORA
    const SimulationModel = SpaceAGORA.SimulationModel
    import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

    function spaceagora_no_gram_example_args()
        planet = SimulationModel.make_no_gram_planet(:earth)
        spacecraft = make_three_body_spacecraft(
            bus_dims=(1.6, 1.8, 1.2),
            panel_dims=(0.01, 0.9, 0.5),
            bus_mass=220.0,
            panel_mass_each=3.0,
            panel_offset_y=0.8,
            ic=SimulationModel.InitialCondition(
                ra=planet.Rp_e + 1_200e3,
                rp=planet.Rp_e + 400e3,
                i=28.5,
                ω=10.0,
                Ω=20.0,
                ν=175.0
            ),
            prop_mass=30.0,
            id=106
        )
        return make_example_config(
            planet=planet,
            spacecraft=spacecraft,
            mission_time=120.0,
            initial_time=SimulationModel.InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
            dynamic_effectors=(SimulationModel.InverseSquaredGravityModel(),),
            density_model=SimulationModel.NoAtmosphereModel(),
            ephemerides_model=SimulationModel.SimpleEphemeridesModel(),
            orientation_sim=false,
            keplerian=true,
            EI_km=120.0,
            verbose=false,
            results=false,
            results_directory=mktempdir()
        )
    end
end; recursive=true)

makedocs(
    sitename = "SpaceAGORA.jl",
    modules = [SpaceAGORA],
    source = DOCS_SRC,
    build = BUILD_DIR,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://space-falcon-lab.github.io/SpaceAGORA.jl",
        assets = ["assets/favicon.ico"],
        size_threshold_warn = 150 * 2^10,
    ),
    pages = [
        "Home" => "index.md",
        "Publications and Citing" => "publications.md",
        "Get Started" => Any[
            "Installation & Environment" => "user/installation_environment.md",
            "Quickstart" => "user/quickstart.md",
        ],
        "Build a Simulation" => Any[
            "Your First Script" => "user/first_script.md",
            "Configuring a Simulation" => "user/simulation_configuration.md",
            "Gravity and Forces" => "user/gravity_forces.md",
            "Adding a Force or Torque" => "user/custom_effector.md",
            "Adding an Actuator" => "user/custom_actuator.md",
            "Powering a Spacecraft" => "user/power_budget.md",
            "Stopping on a Condition" => "user/stop_conditions.md",
            "Mission Workflows" => "user/mission_workflows.md",
            "Solver Configuration" => "user/solver_configuration.md",
            "Parallel Execution" => "user/parallel_execution.md",
        ],
        "Vehicles and Models" => Any[
            "Articulated and Compliant Vehicles" => "user/articulated_compliant.md",
            "Mesh Aerodynamics" => "user/mesh_aerodynamics.md",
            "Constellations" => "user/constellations.md",
            "Plume Interaction" => "user/plume_interaction.md",
            "Interlinks" => "user/interlinks.md",
        ],
        "Atmosphere and Data" => Any[
            "Atmosphere Models" => "user/atmosphere_models.md",
            "GRAM Grid Surrogates" => "user/gram_grid_surrogates.md",
            "Assets & Modes" => "assets.md",
            "GRAMSuite Setup (native GRAM)" => "user/gramsuite_setup.md",
            "Lunar Imagery Sources" => "user/lunar_imagery_sources.md",
        ],
        "Guidance and Planning" => Any[
            "Writing a GNC Algorithm" => "user/gnc_algorithm.md",
            "Feedback Controllers and Estimators" => "user/gnc_feedback.md",
            "Planning Algorithms" => "user/gnc_planning.md",
            "Candidate Control Forecasts" => "user/control_forecasts.md",
            "Your First Controller Package" => "user/feedback_package.md",
            "Energy Depletion Guidance" => "user/edg_guidance.md",
            "Motion Planning with HYPR" => "user/hypr_planning.md",
            "RPO with Interchangeable Planners" => "user/rpo_planner_pilot.md",
        ],
        "Tutorials" => Any[
            "Odyssey Guidance with a Frozen Atmosphere" => "tutorials/odyssey_surrogate.md",
            "HyPR Paper Reconstruction" => "user/hypr_paper.md",
            "Lunar Landing" => "user/lunar_landing.md",
        ],
        "Results" => Any[
            "Simulation Outputs" => "user/outputs.md",
            "Output Columns" => "user/output_columns.md",
            "Interactive Visualization" => "user/visualization.md",
            "Orbital Scenes with MuJoCo" => "user/mujoco_scenes.md",
        ],
        "How It Works" => Any[
            "Concepts" => "user/concepts.md",
            "How a Simulation Fits Together" => "user/architecture.md",
            "The Integrated State" => "user/integrated_state.md",
        ],
        "Studies" => Any[
            "Verification Study" => "user/verification_study.md",
            "Studies and Benchmarks" => "user/studies_benchmarks.md",
        ],
        "Help" => Any[
            "Troubleshooting" => "user/troubleshooting.md",
            "Recipes" => "user/recipes.md",
            "Examples Catalog" => "user/examples_catalog.md",
            "Upgrading" => "user/upgrading.md",
        ],
        "Extend SpaceAGORA" => Any[
            "Extensibility" => "extensibility.md",
            "Extension Conventions" => "reference/extension_conventions.md",
            "Extension Hook Reference" => "generated/extension_hook_reference.md",
            "Migrating Extension Hooks" => "reference/migrating_extension_hooks.md",
            "Building Your Own Research Repository" => "user/research_repository.md",
            "Public API and Compatibility" => "public_api_policy.md",
            "Optional Adapter Contract" => "reference/adapter_contract.md",
            "EDG Service Contract" => "reference/edg_adapter_boundary.md",
            "External Propagation Services" => "reference/external_propagation_services.md",
            "RPO Planner Contract" => "reference/rpo_planner_contract.md",
            "HYPR Service Contract" => "reference/hypr_services.md",
            "Distributed and HPC" => "distributed_hpc.md",
        ],
        "API Reference" => Any[
            "CLI" => "cli.md",
            "Public API" => "generated/public_api.md",
            "Simulation Setup API" => "generated/simulation_setup_api.md",
            "RPO Planner API" => "generated/rpo_planner_api.md",
            "Burn-Plan Family API" => "generated/burn_plan_family_api.md",
            "Joint-Reference Family API" => "generated/joint_reference_family_api.md",
            "EDG Guidance Family API" => "generated/edg_guidance_family_api.md",
            "GNC Algorithm Contract API" => "generated/gnc_algorithm_contract_api.md",
            "GNC Runtime API" => "generated/gnc_runtime_api.md",
            "Robotics API" => "generated/robotics_api.md",
            "Compliant Multibody API" => "generated/compliant_multibody_api.md",
            "Terrain API" => "generated/terrain_api.md",
            "Mesh Aerodynamics API" => "generated/mesh_aerodynamics_api.md",
            "Visualization API" => "generated/visualization_api.md",
            "Lunar Landing API" => "generated/lunar_landing_api.md",
            "Actuators API" => "generated/actuators_api.md",
            "Sensors API" => "generated/sensors_api.md",
            "Extension Hooks API" => "generated/extension_hooks_api.md",
            "Contributor Conformance API" => "generated/contributor_conformance_api.md",
            "Interlinks API" => "generated/interlinks_api.md",
            "Constellations API" => "generated/constellations_api.md",
        ],
    ],
    doctest = true,
    linkcheck = true,
    # The PDS Geosciences Node often times out from CI runners (curl exit 28), which
    # failed the required docs checks of unrelated pull requests. Its two links on
    # docs/src/user/lunar_imagery_sources.md are stable archive addresses: check them
    # by hand when that page changes.
    linkcheck_ignore = [r"^https://pds-geosciences\.wustl\.edu/"],
    checkdocs = :exports,
    checkdocs_ignored_modules = Module[
        SpaceAGORA.RuntimeServices,
        SpaceAGORA.SimulationModel,
        SpaceAGORA.SimulationEngine,
        SpaceAGORA.SimulationCampaigns,
        SpaceAGORA.ParallelProfiles,
        SpaceAGORA.ParallelProcess,
        SpaceAGORA.TelemetryVerification,
        SpaceAGORA.SpaceAGORACLI,
    ],
    warnonly = false,
)

_write_undocumented_public_api_report()

if !isempty(missing_public_api)
    error("Public API documentation is incomplete. See $(PUBLIC_API_REPORT) for the missing exports report.")
end
stale_placeholders = stale_placeholder_entries(SpaceAGORA)
if !isempty(stale_placeholders)
    error("Remove these documented or no longer public names from KNOWN_PLACEHOLDER_DOCS in docs/public_api_symbols.jl: " *
          join([string(owner, ".", symbol) for (owner, symbol) in stale_placeholders], ", "))
end

# CI's strict build needs a GitHub token for link checks, but must not publish.
if !("--no-deploy" in ARGS)
    deploydocs(
        repo = "github.com/Space-FALCON-Lab/SpaceAGORA.jl.git",
        devbranch = "main",
        versions = nothing,
        push_preview = false,
    )
end
