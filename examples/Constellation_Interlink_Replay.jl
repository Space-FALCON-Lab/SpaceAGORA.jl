using SpaceAGORA, LinearAlgebra, CSV, DataFrames, StaticArrays

# Native engine replay of the same terminal-gain plan, with no substitute rate optimizer.
function compare_native_replay(;full=false)
    case=walker27_constellation()
    times=collect(range(0.,(full ? 5. : .05)*case.period_s;length=full ? 1201 : 13))
    plan=plan_laser_schedule(case,times;objective=:a)
    off=TerminalGainPolicy(times,case.keys,zeros(size(plan.command)))
    function native(policy)
        args=constellation_replay_configuration(case,policy)
        run_simulation(args;return_solution=true,extra_callback_factories=(terminal_gain_schedule_callback,))
    end
    reference=native(off);controlled=native(plan.policy)
    replay=replay_laser_schedule(case,times,plan.command)
    OE=SpaceAGORA.SimulationModel.FrameTransforms
    elements(sol)=collect(OE.rvtoorbitalelement(SVector{3}(sol.u[end].sc[1].pos),
        SVector{3}(sol.u[end].sc[1].vel),case.config.planet)[1:5])
    delta=elements(controlled)-elements(reference)
    ks_delta=laser_element_history(case,[last(times)],replay)[:,1]-
        laser_element_history(case,[last(times)],plan.reference)[:,1]
    delta[3:5].=atan.(sin.(delta[3:5]),cos.(delta[3:5]))
    ks_delta[3:5].=atan.(sin.(ks_delta[3:5]),cos.(ks_delta[3:5]))
    scale=[1.,1.,180/pi,180/pi,180/pi]
    out=joinpath(@__DIR__,"..","outputs",full ? "constellation_native_interlinks_5orbits" : "constellation_native_interlinks_smoke")
    mkpath(out)
    table=DataFrame(element=["a_m","e","i_deg","RAAN_deg","omega_deg"],
        linear=scale.*plan.predicted[:,end],KS=scale.*ks_delta,native_Cartesian=scale.*delta,
        native_minus_KS=scale.*(delta-ks_delta))
    CSV.write(joinpath(out,"terminal_comparison.csv"),table)
    history=controlled.prob.p.args.interlink_model.history
    entries=[(;time_s=h.time,emitter=isempty(h.active) ? 0 : maximum(first.(only(h.active)))) for h in history]
    CSV.write(joinpath(out,"native_executed_schedule.csv"),DataFrame(entries))
    for h in history
        k=searchsortedlast(times,h.time)
        expected=k<=length(times)-1 ? case.keys[findall(!iszero,plan.command[:,k])] : InterLinkKey[]
        @assert h.active==expected "Native replay deviated from the selected schedule at $(h.time) s"
    end
    println(table)
    println("Native replay outputs: ",abspath(out))
end
if abspath(PROGRAM_FILE)==@__FILE__
    compare_native_replay(;full="--full" in ARGS)
end
