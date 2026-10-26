module InterlinkContractTests

# The link contract (AbstractInterlink, link_range_m, link_available, link_force),
# the laser momentum link, a custom link kind in a run, and the deprecated names.

using Test
using Logging
using LinearAlgebra
using ComponentArrays
using StaticArrays
include(joinpath(@__DIR__, "..", "..", "..", "examples", "gve_sma_interlinks.jl"))

const SE = SpaceAGORA.SimulationEngine

# A kind with only a range: the default availability and no force.
struct RangeOnlyLink <: AbstractInterlink
    range_m::Float64
end
SpaceAGORA.link_range_m(link::RangeOnlyLink) = link.range_m

# Pulls its endpoints together with a constant force, and closes only while the
# first endpoint is on the +y side of the second.
struct PullLink <: AbstractInterlink
    force_n::Float64
    range_m::Float64
end
SpaceAGORA.link_range_m(link::PullLink) = link.range_m
SpaceAGORA.link_available(link::PullLink, a, b, ctx) =
    a.pos_ii[2] > b.pos_ii[2] && norm(a.pos_ii - b.pos_ii) <= link.range_m
function SpaceAGORA.link_force(link::PullLink, a, b, ctx)
    direction = normalize(b.pos_ii - a.pos_ii)
    return (link.force_n * direction, -link.force_n * direction)
end

# Pushes its first endpoint along its velocity and the second the other way: a
# test kind whose effect on the run is easy to predict.
struct AlongTrackLink <: AbstractInterlink
    force_n::Float64
end
SpaceAGORA.link_range_m(::AlongTrackLink) = 1e9
function SpaceAGORA.link_force(link::AlongTrackLink, a, b, ctx)
    push = link.force_n * normalize(a.vel_ii)
    return (push, -push)
end

# No link_range_m: refused at setup.
struct RangelessLink <: AbstractInterlink end

function pair_state(; first_y=1000.0)
    return ComponentVector(sc=[
        (pos=[7e6, first_y, 0.0], vel=[0.0, 7500.0, 0.0], mass=100.0, interlink_dv=zeros(3), interlink_delta_sma=0.0),
        (pos=[7e6, 0.0, 0.0], vel=[0.0, 7500.0, 0.0], mass=200.0, interlink_dv=zeros(3), interlink_delta_sma=0.0)])
end

# The warnings `f` gives; under --depwarn=error the deprecation throws instead.
function deprecation_messages(f)
    if Base.JLOptions().depwarn == 2
        try
            f()
            return String[]
        catch err
            return [sprint(showerror, err)]
        end
    end
    logger = Test.TestLogger(min_level=Logging.Warn)
    Logging.with_logger(f, logger)
    return [string(entry.message) for entry in logger.logs]
end

@testset "Interlink contract" begin
    @testset "defaults" begin
        context = LinkContext()
        @test isnan(context.time_s) && context.planet === nothing
        @test LinkContext(time_s=2, planet=:earth).time_s === 2.0
        state = pair_state()
        a = SM.InterLinkModels.link_endpoint(state.sc[1], nothing, 1)
        b = SM.InterLinkModels.link_endpoint(state.sc[2], nothing, 2)
        @test a.pos_ii == SVector(7e6, 1000.0, 0.0) && a.mass_kg == 100.0 && a.spacecraft_index == 1
        @test link_available(RangeOnlyLink(1000.0), a, b, context)
        @test !link_available(RangeOnlyLink(999.0), a, b, context)
        @test !link_available(RangeOnlyLink(1000.0), a, a, context)            # coincident endpoints
        @test link_force(RangeOnlyLink(1000.0), a, b, context) == (zeros(3), zeros(3))
        # A partial state without velocity or mass reads NaN there.
        partial = SM.InterLinkModels.link_endpoint((pos=[1.0, 2.0, 3.0],), nothing, 1)
        @test all(isnan, partial.vel_ii) && isnan(partial.mass_kg)
    end

    @testset "the laser momentum link" begin
        link = LaserMomentumLink()
        @test (link.power_w, link.magnification, link.range_m) == (10_000.0, 100.0, 200e3)
        @test link_range_m(LaserMomentumLink(range_m=5.0)) == 5.0
        @test propertynames(link) == (:power_w, :magnification, :range_m)
        terminal = LaserTerminal(boresight=[0.0, 0.0, 1.0], max_power_w=5.0, magnification=2.0, range_m=10.0)
        @test LaserMomentumLink(terminal) == LaserMomentumLink(power_w=5.0, magnification=2.0, range_m=10.0)
        for bad in ((power_w=-1.0,), (magnification=0.0,), (range_m=Inf,))
            @test_throws ArgumentError LaserMomentumLink(; bad...)
        end
        @test_throws ArgumentError LaserMomentumLink(power_w=1.0, P=1.0)
        state = pair_state()
        a = SM.InterLinkModels.link_endpoint(state.sc[1], nothing, 1)
        b = SM.InterLinkModels.link_endpoint(state.sc[2], nothing, 2)
        force_a, force_b = link_force(link, a, b, LinkContext())
        @test force_a ≈ (100.0 * 10_000.0 / 299_792_458.0) * SVector(0.0, 1.0, 0.0)
        @test force_b == -force_a                                                   # equal and opposite
        @test link_force(link, a, a, LinkContext()) == (zeros(3), zeros(3))
    end

    @testset "a custom link kind" begin
        spacecraft = [SpacecraftModel(id=index) for index in 1:2]
        model = InterlinkModel()
        key = register_link!(model, spacecraft, (2, 1), (1, 1), PullLink(3.0, 2000.0))
        @test key == ((1, 1), (2, 1)) && model.graph.edge_properties == [PullLink(3.0, 2000.0)]
        state = pair_state()
        @test update_availability!(model, spacecraft, state) == [true]
        @test update_availability!(model, spacecraft, pair_state(first_y=-1000.0)) == [false]   # its own rule
        @test update_availability!(model, spacecraft, state; is_active=[true, false]) == [false]  # the engine's
        update_availability!(model, spacecraft, state)
        model.graph.active[1] = true
        @test interlink_force_on_spacecraft(model, state, 1) ≈ SVector(0.0, -3.0, 0.0)
        @test interlink_force_on_spacecraft(model, state, 2) ≈ SVector(0.0, 3.0, 0.0)
        @test interlink_force_on_spacecraft(model, state, 1; is_active=[true, false]) == zeros(3)
        # The policies score by the link's own force: pulling the first
        # spacecraft against its velocity lowers its semimajor axis, while the
        # second is pulled along its velocity.
        mu = make_no_gram_planet(:earth).μ
        first_target = SchedulingPolicyModel(SemimajorAxisRatePolicy(targets=1))
        score_candidates!(model, first_target, state, mu)
        @test first_target.matching.score[1] < 0
        second_target = SchedulingPolicyModel(SemimajorAxisRatePolicy(targets=2))
        score_candidates!(model, second_target, state, mu)
        @test second_target.matching.score[1] > 0
        @test schedule_interlinks!(model, second_target, spacecraft, state, mu) == [key]
    end

    @testset "a run applies each link's own force" begin
        # One link from spacecraft 1, pushing it along its velocity: the scheduler
        # selects it (it raises the target's semimajor axis) and the integrated
        # interlink acceleration of spacecraft 1 grows as force / mass * time.
        args = gve_sma_configuration(duration=60.0)
        model = InterlinkModel()
        register_link!(model, args.dynamics_model.spacecraft, (1, 1), (2, 1), AlongTrackLink(2.0))
        args = SM.SimConfig._with_configuration(args; interlink_model=model)
        solution = run_simulation(args; return_solution=true)
        history = solution.prob.p.args.interlink_model.history
        @test all(sample -> sample.active == [((1, 1), (2, 1))], history)
        mass = solution.u[1].sc[1].mass
        @test norm(solution.u[end].sc[1].interlink_dv) ≈ 2.0 / mass * 60.0 rtol=1e-3
        @test norm(solution.u[end].sc[2].interlink_dv) > 0
        @test norm(solution.u[end].sc[3].interlink_dv) == 0                            # not linked
    end

    @testset "a kind without link_range_m is refused at setup" begin
        args = gve_sma_configuration(duration=10.0)
        model = InterlinkModel()
        register_link!(model, args.dynamics_model.spacecraft, (1, 1), (2, 1), RangelessLink())
        err = try
            run_simulation(SM.SimConfig._with_configuration(args; interlink_model=model))
            nothing
        catch caught
            caught
        end
        @test err isa ArgumentError && occursin("link_range_m", sprint(showerror, err))
    end

    @testset "the deprecated names" begin
        @test ConstellationInterLinkModel === InterlinkModel
        @test InterLinkParameters === LaserMomentumLink
        spacecraft = [SpacecraftModel(id=index) for index in 1:2]
        old_model = Ref{Any}(nothing)
        messages = deprecation_messages() do
            model = ConstellationInterLinkModel(constellation_link_type=:laser)
            register_allowed_link!(model, spacecraft, (1, 1), (2, 1);
                parameters=InterLinkParameters(P=12.0, B=3.0, range=100.0))
            old_model[] = model
        end
        if Base.JLOptions().depwarn == 2
            @test length(messages) == 1 && occursin("deprecated", messages[1])
        else
            @test any(m -> occursin("constellation_link_type keyword is deprecated", m), messages)
            @test any(m -> occursin("keyword P is deprecated", m), messages)
            @test any(m -> occursin("register_allowed_link! is deprecated", m), messages)
            new_model = InterlinkModel()
            register_link!(new_model, spacecraft, (1, 1), (2, 1),
                LaserMomentumLink(power_w=12.0, magnification=3.0, range_m=100.0))
            @test old_model[].graph.edge_properties == new_model.graph.edge_properties
            @test old_model[].graph.edges == new_model.graph.edges
            link = new_model.graph.edge_properties[1]
            state = pair_state()
            fields = deprecation_messages(() -> (link.P, link.B, link.range))
            @test length(fields) == 3 && all(m -> occursin("is deprecated", m), fields)
            @test (link.power_w, link.magnification, link.range_m) == (12.0, 3.0, 100.0)
            forces = Ref{Any}(nothing)
            calls = deprecation_messages() do
                new_model.graph.active[1] = true
                forces[] = (force_on_endpoint(state.sc[1], state.sc[2], link),
                    laser_force_on_spacecraft(new_model, state, 1), new_model.constellation_link_type)
            end
            @test length(calls) == 3
            a = SM.InterLinkModels.link_endpoint(state.sc[1], nothing, 1)
            b = SM.InterLinkModels.link_endpoint(state.sc[2], nothing, 2)
            @test forces[][1] == link_force(link, a, b, LinkContext())[1]
            @test forces[][2] == interlink_force_on_spacecraft(new_model, state, 1)
            @test forces[][3] === :laser
            radio = deprecation_messages(() -> @test_throws ArgumentError InterlinkModel(constellation_link_type=:radio))
            @test !isempty(radio)
        end
    end
end

end # module
