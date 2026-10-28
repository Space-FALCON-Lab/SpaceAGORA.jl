# Each invocation is a fresh process; checkpoint resume must survive that boundary.
ENV["SPACEAGORA_DENSITY_FREEZE_PER_STEP"]="false"
ENV["SPACEAGORA_GRAM_TRACK_CACHE_MODE"]="off"
include("lifecycle-fixtures.jl")
using KSPropagator
length(ARGS)>=3 || error("usage: lifecycle-worker.jl MODE PHASE DIRECTORY [coincident]")
mode=Symbol(ARGS[1]);phase=Symbol(ARGS[2]);directory=abspath(ARGS[3])
coincident=length(ARGS)>=4 && ARGS[4]=="coincident"
checkpoint=phase!==:direct
resume=phase===:resume
interrupt=phase===:interrupt
mkpath(directory)
cfg=lifecycle_config(mode;directory,checkpoint,resume,coincident)
factory=lifecycle_callbacks(;interrupt,coincident)
try
    out=run_simulation(cfg;return_results=true,return_solution=!checkpoint,
        priority_callback_factories=(factory,),isolate_state=!interrupt,save_fields=lifecycle_save_fields(cfg))
    c=out.configuration.control_model.control_effectors[1]
    final=checkpoint ? Vector(SE._load_checkpoint(out.configuration).u) : Vector(out.solution.u[end])
    record=(mode=mode,phase=phase,coincident=coincident,table=out.table,final=final,
        ticks=copy(c.ticks),actions=copy(c.actions),roots=deepcopy(c.roots),clocks=copy(c.clocks),trace=copy(c.trace),
        root_alias=c.root===only(out.configuration.dynamics_model.spacecraft).root,
        initial_audits=deepcopy(c.initial_audits),
        area=c.root.ref_area,held_force_N=c.held_force_N,
        checkpoint_time=checkpoint ? SE._load_checkpoint(out.configuration).t : nothing,
        returned_solution=!checkpoint,source_julia=string(VERSION))
    serialize(joinpath(directory,"$(mode)-$(phase)-result.bin"),record)
    println("LIFECYCLE_WORKER_SUCCESS mode=",mode," phase=",phase," rows=",size(out.table,1),
        " ticks=",length(c.ticks)," actions=",c.actions," roots=",length(c.roots),
        " final_mass_kg=",final[7]," final_heat_J_cm2=",final[8]);flush(stdout)
catch err
    if interrupt && occursin("deliberate lifecycle interruption",sprint(showerror,err))
        checkpoint_data=SE._load_checkpoint(cfg)
        c=cfg.control_model.control_effectors[1]
        serialize(joinpath(directory,"$(mode)-interrupted-audit.bin"),
            (error=sprint(showerror,err),checkpoint_time=checkpoint_data.t,
             observed_ticks=copy(c.ticks),observed_actions=copy(c.actions),observed_roots=deepcopy(c.roots),
             discarded_trace=copy(c.trace),checkpoint_snapshot=checkpoint_data,
             note="Discarded audit tail is not authoritative resumed history."))
        println("LIFECYCLE_EXPECTED_INTERRUPTION checkpoint_time=",checkpoint_data.t);flush(stdout)
    else
        rethrow()
    end
end
