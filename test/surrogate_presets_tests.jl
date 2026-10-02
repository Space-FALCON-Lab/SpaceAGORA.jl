using Test, SpaceAGORA, GRAMSuite, SHA, TOML, Serialization, Libdl, Artifacts, Pkg

const PRESET_ENV = SpaceAGORA.SimulationModel.EnvironmentModels
function preset_write_toml(path, record)
    open(path,"w") do io; TOML.print(io,record;sorted=true);end
end
preset_sha(path) = open(io->bytes2hex(sha256(io)),path)

@testset "Named surrogate preset integrity and native-free resolution" begin
    @test !any(x->occursin("libgram",lowercase(basename(x))),Libdl.dllist())
    catalog_template = TOML.parsefile(PRESET_ENV._SURROGATE_CATALOG)
    @test any(x->x["id"]=="odyssey_p20_frozen_v1" && x["version"]=="1.0.0",SpaceAGORA.available_surrogate_presets())
    mktempdir() do dir
        file=joinpath(dir,"fixture.jls"); catalog_file=joinpath(dir,"catalog.toml"); artifacts_file=joinpath(dir,"Artifacts.toml")
        catalog=deepcopy(catalog_template); entry=only(filter(p->p["id"]=="odyssey_p20_frozen_v1",catalog["presets"])); catalog["presets"]=[entry]
        entry["id"]="synthetic_preset";entry["release_enabled"]=false;entry["artifact_name"]="synthetic_preset_1_0_0"
        entry["axes"]=Dict("altitude"=>Dict("start"=>100.,"step"=>160.,"count"=>2),"latitude"=>Dict("start"=>40.,"step"=>50.,"count"=>2),"longitude"=>Dict("start"=>0.,"step"=>180.,"count"=>2))
        payload=deepcopy(entry["required_metadata"])
        payload["grid"]=Dict("alt_km"=>[100.,260.],"lat_deg"=>[40.,90.],"lon_deg"=>[0.,180.])
        payload["fields"]=Dict(k=>fill(v,2,2,2) for (k,v) in zip(("density_kgm3","temperature_K","wind_ew_ms","wind_ns_ms","wind_up_ms"),(1e-8,180.,2.,3.,4.)))
        function save_payload!()
            serialize(file,payload);entry["payload"]["sha256"]=preset_sha(file);entry["payload"]["bytes"]=filesize(file);preset_write_toml(catalog_file,catalog)
        end
        save_payload!();preset_write_toml(artifacts_file,Dict())
        defaults=(version="1.0.0",catalog_file=catalog_file,artifacts_file=artifacts_file)
        local_defaults=(;defaults...,file=file,allow_unreleased=true)
        @test_throws ArgumentError SpaceAGORA.resolve_surrogate_preset("missing";defaults...)
        @test_throws ArgumentError SpaceAGORA.resolve_surrogate_preset("synthetic_preset";merge(defaults,(version="latest",))...)
        @test_throws ArgumentError SpaceAGORA.resolve_surrogate_preset("synthetic_preset";defaults...)
        @test_throws ArgumentError SpaceAGORA.resolve_surrogate_preset("synthetic_preset";defaults...,allow_unreleased=true)
        @test_throws ArgumentError SpaceAGORA.resolve_surrogate_preset("synthetic_preset";defaults...,file=file)
        @test_throws ArgumentError SpaceAGORA.resolve_surrogate_preset("synthetic_preset";local_defaults...,planet="Earth")
        resolved=SpaceAGORA.resolve_surrogate_preset("synthetic_preset";local_defaults...,planet="MARS")
        @test resolved.file==file && resolved.planet=="mars"
        @test resolved.expected_sha256==preset_sha(file)
        @test resolved.provenance["local_development"]
        model=SpaceAGORA.surrogate_preset_model("synthetic_preset";local_defaults...)
        @test model isa SpaceAGORA.GRAMGridAtmosphereModel
        @test SpaceAGORA.getDensity(model,150000.,deg2rad(60.),deg2rad(359.),0.,true)[1]≈1e-8
        @test SpaceAGORA.getDensity(model,150000.,deg2rad(60.),0.,1e9,false)==SpaceAGORA.getDensity(model,150000.,deg2rad(60.),0.,0.,true)
        @test_throws DomainError SpaceAGORA.getDensity(model,99999.,deg2rad(60.),0.,0.,true)
        @test_throws DomainError SpaceAGORA.getDensity(model,150000.,deg2rad(39.),0.,0.,true)
        @test_throws DomainError SpaceAGORA.getDensity(model,260001.,deg2rad(60.),0.,0.,true)
        # A named preset fixes its domain policy: grid options raise a clear ArgumentError, not a MethodError.
        policy_error=try SpaceAGORA.surrogate_preset_model("synthetic_preset";local_defaults...,above_grid=:vacuum) catch err err end
        @test policy_error isa ArgumentError
        @test occursin("above_grid",policy_error.msg) && occursin("fixes its domain policy",policy_error.msg) && occursin("GRAMGridAtmosphereModel",policy_error.msg)
        @test Set(Base.kwarg_decl(only(methods(SpaceAGORA.resolve_surrogate_preset))))==Set(PRESET_ENV._PRESET_RESOLUTION_KEYWORDS)
        # The generic grid model accepts the vacuum policy above the ceiling; below the floor still fails.
        vacuum=SpaceAGORA.GRAMGridAtmosphereModel(planet="Mars",surrogate_file=resolved.file,above_grid=:vacuum)
        @test SpaceAGORA.getDensity(vacuum,260001.,deg2rad(60.),0.,0.,true)[1]==0
        @test_throws DomainError SpaceAGORA.getDensity(vacuum,99999.,deg2rad(60.),0.,0.,true)
        before=SpaceAGORA.atmosphere_provenance(model)
        @test before["catalog_sha256"]==preset_sha(catalog_file)
        before["domain"]["height_m"][1]=0
        @test SpaceAGORA.atmosphere_provenance(model)["domain"]["height_m"][1]==100000.
        @test SpaceAGORA.atmosphere_provenance(SpaceAGORA.NoAtmosphereModel())["backend"]=="configured_density_model"
        raw=SpaceAGORA.GRAMGridAtmosphereModel(planet="Mars",surrogate_file=file)
        @test SpaceAGORA.atmosphere_provenance(raw)["preset_status"]=="user_supplied_grid_without_named_preset_contract"
        pristine=read(file)
        for badfile in (joinpath(dir,"absent"),dir)
            @test_throws ArgumentError SpaceAGORA.resolve_surrogate_preset("synthetic_preset";merge(local_defaults,(file=badfile,))...)
        end
        for corrupt in (UInt8[],Vector{UInt8}(codeunits("version https://git-lfs.github.com/spec/v1\n")),vcat(pristine,0x00),vcat(pristine[1:end-1],pristine[end]⊻0x01))
            write(file,corrupt)
            @test_throws ArgumentError SpaceAGORA.resolve_surrogate_preset("synthetic_preset";local_defaults...)
            @test read(file)==corrupt
        end
        write(file,pristine)
        # Re-pin each mutated fixture so integrity passes and metadata/domain checks must reject it.
        for mutate in (p->p["generation_config"]["is_planetocentric"]=true,
            p->p["generation_config"]["mars_map_year"]=1,
            p->p["initial_time"]["day"]=8,p->delete!(p,"generation_config"),
            p->p["grid"]["lat_deg"][1]=39.,p->p["grid"]["alt_km"][2]=261.,
            p->p["planet"]="earth",p->p["format"]="unknown")
            original=deepcopy(payload);mutate(payload);save_payload!()
            @test_throws ArgumentError SpaceAGORA.surrogate_preset_model("synthetic_preset";local_defaults...)
            payload=original;save_payload!()
        end
        push!(catalog["presets"],deepcopy(entry));preset_write_toml(catalog_file,catalog)
        @test_throws ArgumentError SpaceAGORA.available_surrogate_presets(;catalog_file)
        pop!(catalog["presets"]);preset_write_toml(catalog_file,catalog)
        # Synthetic lazy artifacts exercise the real Julia artifact machinery, with no network server.
        Artifacts.with_artifacts_directory(joinpath(dir,"artifact_store")) do
            hash=Pkg.Artifacts.create_artifact() do target; cp(file,joinpath(target,"mars_grid.jls"));end
            archive=joinpath(dir,"fixture.tar.gz");archive_sha=Pkg.Artifacts.archive_artifact(hash,archive)
            entry["release_enabled"]=true
            entry["distribution"]=Dict("status"=>"synthetic_test_only","git_tree_sha1"=>string(hash),"archive_sha256"=>archive_sha,"urls"=>["file://"*archive])
            Pkg.Artifacts.bind_artifact!(artifacts_file,entry["artifact_name"],hash;lazy=true,download_info=[("file://"*archive,archive_sha)],force=true)
            preset_write_toml(catalog_file,catalog)
            installed=SpaceAGORA.surrogate_preset_model("synthetic_preset";defaults...,offline=true)
            @test SpaceAGORA.atmosphere_provenance(installed)["resolution"]=="julia_artifact"
            @test SpaceAGORA.atmosphere_provenance(installed)["artifact_git_tree_sha1"]==string(hash)
            # Explicit invalid file must not use the good installed artifact.
            @test_throws ArgumentError SpaceAGORA.resolve_surrogate_preset("synthetic_preset";defaults...,file=joinpath(dir,"bad"),offline=true)
            artifact_dir=Artifacts.artifact_path(hash);chmod(artifact_dir,0o555);chmod(joinpath(artifact_dir,"mars_grid.jls"),0o444)
            try
                @test SpaceAGORA.resolve_surrogate_preset("synthetic_preset";defaults...,offline=true).expected_sha256==preset_sha(file)
            finally
                chmod(artifact_dir,0o755);chmod(joinpath(artifact_dir,"mars_grid.jls"),0o644)
            end
            rm(artifact_dir;recursive=true)
            @test_throws ArgumentError SpaceAGORA.resolve_surrogate_preset("synthetic_preset";defaults...,offline=true)
            withenv("JULIA_PKG_SERVER"=>"") do
                # One line before and one after installation; Pkg's own status lines stay off the terminal.
                fetched=nothing
                terminal=mktemp() do path,io
                    redirect_stderr(io) do
                        fetched=@test_logs (:info,r"^Installing atmosphere preset synthetic_preset@1\.0\.0 .* from file://") (:info,r"^Installed atmosphere preset synthetic_preset@1\.0\.0 in .*SHA256") SpaceAGORA.resolve_surrogate_preset("synthetic_preset";defaults...)
                    end
                    close(io);read(path,String)
                end
                @test !occursin("artifact:",terminal)
                @test isfile(fetched.file) && preset_sha(fetched.file)==entry["payload"]["sha256"]
                @test !isdir(joinpath(dirname(artifact_dir),"unselected_artifact"))
                rm(artifact_dir;recursive=true)
                tasks=[Threads.@spawn SpaceAGORA.resolve_surrogate_preset("synthetic_preset";defaults...) for _ in 1:2]
                @test all(r->r.expected_sha256==entry["payload"]["sha256"],fetch.(tasks))
                rm(artifact_dir;recursive=true)
                archive_bytes=read(archive);write(archive,archive_bytes[1:10])
                install_error=try SpaceAGORA.resolve_surrogate_preset("synthetic_preset";defaults...) catch err err end
                @test install_error isa ArgumentError
                # The raised error keeps Pkg's diagnostics: the attempted source and its captured status lines.
                @test occursin("file://"*archive,install_error.msg) && occursin("Pkg output:",install_error.msg)
                @test !Artifacts.artifact_exists(hash)
                write(archive,archive_bytes)
            end
            # A Julia artifact override is honored but cannot bypass payload identity.
            override=joinpath(dir,"override");mkpath(override);write(joinpath(override,"mars_grid.jls"),"wrong")
            overrides=joinpath(dirname(artifact_dir),"Overrides.toml")
            preset_write_toml(overrides,Dict(string(hash)=>override));Artifacts.load_overrides(;force=true)
            try
                @test_throws ArgumentError SpaceAGORA.resolve_surrogate_preset("synthetic_preset";defaults...,offline=true)
                cp(file,joinpath(override,"mars_grid.jls");force=true)
                @test SpaceAGORA.resolve_surrogate_preset("synthetic_preset";defaults...,offline=true).file==joinpath(override,"mars_grid.jls")
            finally
                rm(overrides);Artifacts.load_overrides(;force=true)
            end
        end
        @test !any(x->occursin("libgram",lowercase(basename(x))),Libdl.dllist())
    end
end

@testset "Named presets with explicit non-uniform axis nodes" begin
    catalog_template = TOML.parsefile(PRESET_ENV._SURROGATE_CATALOG)
    @test any(x->x["id"]=="mars_global_upper_p20_frozen_v1" && x["version"]=="1.0.0" && x["release_enabled"],SpaceAGORA.available_surrogate_presets())
    shipped=only(filter(p->p["id"]=="mars_global_upper_p20_frozen_v1",catalog_template["presets"]))
    @test length(shipped["axes"]["altitude"]["nodes"])==253 && length(shipped["axes"]["latitude"]["nodes"])==229
    @test all(>(0),diff(shipped["axes"]["latitude"]["nodes"])) && shipped["axes"]["longitude"]["count"]==144
    mktempdir() do dir
        file=joinpath(dir,"fixture.jls"); catalog_file=joinpath(dir,"catalog.toml"); artifacts_file=joinpath(dir,"Artifacts.toml")
        catalog=deepcopy(catalog_template); entry=deepcopy(shipped); catalog["presets"]=[entry]
        entry["id"]="synthetic_nodes";entry["release_enabled"]=false;entry["artifact_name"]="synthetic_nodes_1_0_0"
        alt=[80.,95.5,365.]; lat=[-90.,-82.36,0.,90.]
        entry["axes"]=Dict("altitude"=>Dict("nodes"=>alt),"latitude"=>Dict("nodes"=>lat),"longitude"=>Dict("start"=>0.,"step"=>180.,"count"=>2))
        payload=deepcopy(entry["required_metadata"])
        payload["grid"]=Dict("alt_km"=>copy(alt),"lat_deg"=>copy(lat),"lon_deg"=>[0.,180.])
        payload["fields"]=Dict(k=>fill(v,3,4,2) for (k,v) in zip(("density_kgm3","temperature_K","wind_ew_ms","wind_ns_ms","wind_up_ms"),(1e-8,180.,2.,3.,4.)))
        function save_payload!()
            serialize(file,payload);entry["payload"]["sha256"]=preset_sha(file);entry["payload"]["bytes"]=filesize(file);preset_write_toml(catalog_file,catalog)
        end
        save_payload!();preset_write_toml(artifacts_file,Dict())
        opts=(version="1.0.0",catalog_file=catalog_file,artifacts_file=artifacts_file,file=file,allow_unreleased=true)
        model=SpaceAGORA.surrogate_preset_model("synthetic_nodes";opts...)
        @test SpaceAGORA.getDensity(model,90000.,deg2rad(-85.),deg2rad(90.),0.,true)[1]≈1e-8
        @test_throws DomainError SpaceAGORA.getDensity(model,79999.,deg2rad(10.),0.,0.,true)
        @test_throws DomainError SpaceAGORA.getDensity(model,365001.,deg2rad(10.),0.,0.,true)
        # A payload whose nodes differ from the declared nodes is rejected.
        payload["grid"]["lat_deg"][2]=-82.35;save_payload!()
        @test_throws ArgumentError SpaceAGORA.surrogate_preset_model("synthetic_nodes";opts...)
        payload["grid"]["lat_deg"][2]=-82.36;save_payload!()
        @test SpaceAGORA.surrogate_preset_model("synthetic_nodes";opts...) isa SpaceAGORA.GRAMGridAtmosphereModel
        # Invalid node declarations fail when the catalog is read.
        for (axis,bad) in (("altitude",Dict("nodes"=>[80.,80.,365.])),("altitude",Dict("nodes"=>[365.,80.])),("altitude",Dict("nodes"=>[80.])),
                           ("altitude",Dict("nodes"=>[80.,365.],"step"=>1.)),("latitude",Dict("nodes"=>[-90.,NaN,90.])),
                           ("longitude",Dict("nodes"=>[0.,180.])))
            broken=deepcopy(catalog);only(broken["presets"])["axes"][axis]=bad;preset_write_toml(catalog_file,broken)
            @test_throws ArgumentError SpaceAGORA.available_surrogate_presets(;catalog_file)
        end
        # The declared domain must equal the first and last nodes.
        broken=deepcopy(catalog);only(broken["presets"])["domain"]["height_m"]=[80000.,366000.];preset_write_toml(catalog_file,broken)
        @test_throws ArgumentError SpaceAGORA.available_surrogate_presets(;catalog_file)
    end
end

@testset "Named near-surface presets" begin
    catalog_template = TOML.parsefile(PRESET_ENV._SURROGATE_CATALOG)
    shipped=only(filter(p->p["id"]=="mars_global_near_surface_p20_frozen_v1" && p["version"]=="1.0.0",catalog_template["presets"]))
    @test shipped["kind"]=="gram_near_surface_scalars" && shipped["release_enabled"] && !haskey(shipped,"axes")
    @test shipped["atmosphere"]["winds_available"]===false && shipped["domain"]["top_areoid_height_m"]==75000.0
    @test any(x->x["id"]=="mars_global_near_surface_p20_frozen_v1" && x["version"]=="1.0.0",SpaceAGORA.available_surrogate_presets())
    mktempdir() do dir
        file=joinpath(dir,"fixture.jls"); catalog_file=joinpath(dir,"catalog.toml"); artifacts_file=joinpath(dir,"Artifacts.toml")
        catalog=deepcopy(catalog_template); entry=deepcopy(shipped); catalog["presets"]=[entry]
        entry["id"]="synthetic_near_surface";entry["release_enabled"]=false;entry["artifact_name"]="synthetic_near_surface_1_0_0"
        # Synthetic analytic payload on the shipped lattice (not Mars data): flat 0.5 km terrain, linear level states.
        payload=deepcopy(entry["required_metadata"]); lev=payload["levels_km"]; g=payload["lattice"]
        nl,ni,nj=length(lev),g["nlat"],g["nlon"]
        payload["radii_km"]=(payload["generation_config"]["equatorial_radius_km"],payload["generation_config"]["polar_radius_km"])
        payload["terrain"]=Dict{String,Any}("lat0_deg"=>-86.25,"lon0_deg"=>0.488,"step_deg"=>0.5,
            "surface_height_km"=>fill(0.5,346,720),"areoid_radius_km"=>fill(3390.0,346,720))
        payload["level_T_K"]=[220.0-2lev[a] for a in 1:nl, i in 1:ni, j in 1:nj]
        payload["level_R"]=fill(191.0,nl,ni,nj); payload["level_lnp"]=[log(700.0)-lev[a]/11 for a in 1:nl, i in 1:ni, j in 1:nj]
        payload["level_source"]=ones(UInt8,nl,ni,nj); payload["surface_T30_K"]=fill(214.0,ni,nj); payload["surface_T5_K"]=fill(216.0,ni,nj)
        keys_=[(b,c) for b in -12:11 for c in 0:39]; n=length(keys_)
        payload["q_models"]=Dict{String,Any}("band"=>first.(keys_),"cell"=>last.(keys_),"L"=>fill(1,n),"order"=>fill(1,n),
            "phic_center"=>[7.5b+3.75 for (b,_) in keys_],"lam_center"=>[9.0c+4.5 for (_,c) in keys_],
            "coef"=>hcat(fill(20.0,n),zeros(n,5)),"n_points"=>fill(1,n),"status"=>fill("qualified",n))
        function save_payload!()
            serialize(file,payload);entry["payload"]["sha256"]=preset_sha(file);entry["payload"]["bytes"]=filesize(file);preset_write_toml(catalog_file,catalog)
        end
        save_payload!();preset_write_toml(artifacts_file,Dict())
        opts=(version="1.0.0",catalog_file=catalog_file,artifacts_file=artifacts_file,file=file,allow_unreleased=true)
        # Invalid near-surface declarations fail when the catalog is read.
        for breakit in (e->e["kind"]="gram_mystery", e->e["axes"]=deepcopy(only(filter(p->haskey(p,"axes"),catalog_template["presets"][1:1]))["axes"]),
                        e->e["domain"]["top_areoid_height_m"]=80000.0, e->e["domain"]["outside"]="clamp",
                        e->e["atmosphere"]["winds_available"]=true, e->e["domain"]["minimum_clearance_m"]=NaN,
                        e->e["required_metadata"]["format"]="spaceagora_gram_static_grid_v1")
            broken=deepcopy(catalog);breakit(only(broken["presets"]));preset_write_toml(catalog_file,broken)
            @test_throws ArgumentError SpaceAGORA.available_surrogate_presets(;catalog_file)
        end
        preset_write_toml(catalog_file,catalog)
        # The wrapper's bypass membership and thread-safety trait do not need GRAMSuite's near-surface API, so they are
        # checked on a stub core in every environment.
        wrapper=SpaceAGORA.GRAMNearSurfaceAtmosphereModel(:stub_core)
        @test wrapper isa PRESET_ENV._NativeFreeSnapshotModel && wrapper.core===:stub_core
        @test SpaceAGORA.SimulationModel.SimulationCallbacks.density_model_threadsafe(wrapper)
        # The model needs a GRAMSuite with the near-surface API. With an older GRAMSuite (for example a CI job pinned to a
        # revision without it), the named preset must fail with a clear error, and the model checks are reported skipped.
        if !isdefined(GRAMSuite, :GRAMNearSurfaceAtmosphereModel)
            missing_api=try SpaceAGORA.surrogate_preset_model("synthetic_near_surface";opts...) catch err err end
            @test missing_api isa ArgumentError && occursin("near-surface API",missing_api.msg)
            @test_skip "near-surface model checks need a GRAMSuite with GRAMNearSurfaceAtmosphereModel"
            return
        end
        model=SpaceAGORA.surrogate_preset_model("synthetic_near_surface";opts...)
        @test model isa SpaceAGORA.GRAMNearSurfaceAtmosphereModel
        @test model isa PRESET_ENV._NativeFreeSnapshotModel
        @test SpaceAGORA.SimulationModel.SimulationCallbacks.density_model_threadsafe(model)
        # 3 km above the ellipsoid at 10 N is well inside the synthetic domain (areoid 3390 km, terrain 0.5 km)
        rho,T,wind=SpaceAGORA.getDensity(model,3000.,deg2rad(10.),deg2rad(40.),0.,true)
        @test (rho,T,wind)==GRAMSuite.density_state(model.core,3000.,deg2rad(10.),deg2rad(40.),0.,true)
        @test wind==zeros(3) && rho>0 && T>0
        @test SpaceAGORA.getDensity(model,3000.,deg2rad(10.),deg2rad(40.),1e9,false)==(rho,T,wind)
        s=GRAMSuite.near_surface_state(model.core,10.,40.,3000.)
        @test s.density_kgm3===rho && s.pressure_Pa>0
        @test_throws DomainError SpaceAGORA.getDensity(model,90000.,deg2rad(10.),deg2rad(40.),0.,true)
        @test_throws DomainError SpaceAGORA.getDensity(model,3000.,deg2rad(86.),deg2rad(40.),0.,true)
        @test_throws DomainError SpaceAGORA.getDensity(model,-20000.,deg2rad(10.),deg2rad(40.),0.,true)
        provenance=SpaceAGORA.atmosphere_provenance(model)
        @test provenance["backend"]=="gram_near_surface_surrogate" && provenance["preset_kind"]=="gram_near_surface_scalars"
        @test provenance["domain"]["minimum_clearance_m"]==5.0 && provenance["atmosphere"]["wind_returned"]=="zero_vector"
        raw=SpaceAGORA.GRAMNearSurfaceAtmosphereModel(planet="Mars",surrogate_file=file)
        @test SpaceAGORA.atmosphere_provenance(raw)["preset_status"]=="user_supplied_payload_without_named_preset_contract"
        policy_error=try SpaceAGORA.surrogate_preset_model("synthetic_near_surface";opts...,above_grid=:vacuum) catch err err end
        @test policy_error isa ArgumentError && occursin("above_grid",policy_error.msg)
        # Payload metadata that differs from the catalog contract is rejected after its identity is re-pinned.
        for mutate in (p->p["support"]["top_areoid_km"]=70.0, p->p["generation_config"]["mars_map_year"]=1,
                       p->p["format"]="unknown", p->p["planet"]="earth", p->p["epoch_utc"]="2001-11-07T11:51:05Z",
                       p->p["distribution"]["preset_id"]="another")
            original=deepcopy(payload);mutate(payload);save_payload!()
            @test_throws ArgumentError SpaceAGORA.surrogate_preset_model("synthetic_near_surface";opts...)
            payload=original;save_payload!()
        end
        # The released path installs a lazy artifact and builds the near-surface model from it.
        Artifacts.with_artifacts_directory(joinpath(dir,"artifact_store")) do
            hash=Pkg.Artifacts.create_artifact() do target; cp(file,joinpath(target,"mars_near_surface.jls"));end
            archive=joinpath(dir,"fixture.tar.gz");archive_sha=Pkg.Artifacts.archive_artifact(hash,archive)
            entry["release_enabled"]=true
            entry["distribution"]=Dict("status"=>"synthetic_test_only","git_tree_sha1"=>string(hash),"archive_sha256"=>archive_sha,"urls"=>["file://"*archive])
            Pkg.Artifacts.bind_artifact!(artifacts_file,entry["artifact_name"],hash;lazy=true,download_info=[("file://"*archive,archive_sha)],force=true)
            preset_write_toml(catalog_file,catalog)
            installed=SpaceAGORA.surrogate_preset_model("synthetic_near_surface";version="1.0.0",catalog_file,artifacts_file,offline=true)
            @test installed isa SpaceAGORA.GRAMNearSurfaceAtmosphereModel
            @test SpaceAGORA.atmosphere_provenance(installed)["resolution"]=="julia_artifact"
            @test SpaceAGORA.getDensity(installed,3000.,deg2rad(10.),deg2rad(40.),0.,true)==(rho,T,wind)
        end
        @test !any(x->occursin("libgram",lowercase(basename(x))),Libdl.dllist())
    end
end

@testset "Named near-surface preset versions" begin
    catalog_template = TOML.parsefile(PRESET_ENV._SURROGATE_CATALOG)
    near=filter(p->p["id"]=="mars_global_near_surface_p20_frozen_v1",catalog_template["presets"])
    @test sort([p["version"] for p in near])==["1.0.0","1.1.0"]
    v100=only(filter(p->p["version"]=="1.0.0",near)); v110=only(filter(p->p["version"]=="1.1.0",near))
    listed=SpaceAGORA.available_surrogate_presets()
    @test all(v->any(x->x["id"]=="mars_global_near_surface_p20_frozen_v1" && x["version"]==v,listed),("1.0.0","1.1.0"))
    # 1.0.0 keeps its published identities; 1.1.0 extends it to 81 km with two added levels and has its own identities.
    @test v100["domain"]["top_areoid_height_m"]==75000.0 && length(v100["required_metadata"]["levels_km"])==29
    @test v100["distribution"]["archive_sha256"]=="f0fa4a4cd05b747cef77f723bdb2a1a9697d8b824ee8b6159095c0bcb9212b37"
    @test v110["kind"]=="gram_near_surface_scalars" && v110["release_enabled"] && !haskey(v110,"axes")
    @test v110["atmosphere"]["winds_available"]===false && v110["atmosphere"]["wind_returned"]=="zero_vector"
    @test v110["domain"]["top_areoid_height_m"]==81000.0 && v110["required_metadata"]["support"]["top_areoid_km"]==81.0
    lev0,lev1=v100["required_metadata"]["levels_km"],v110["required_metadata"]["levels_km"]
    @test lev1[1:length(lev0)]==lev0 && lev1[end-1:end]==[80.03231545290365,85.03231545290365]
    for k in ("lattice","generation_config","epoch_utc","format","planet")
        @test v110["required_metadata"][k]==v100["required_metadata"][k]
    end
    @test v110["artifact_name"]=="mars_global_near_surface_p20_frozen_v1_1_1_0" && v110["artifact_name"]!=v100["artifact_name"]
    @test v110["payload"]["sha256"]=="5abb86083c0d9e824bf7a2f31203426834bb9c1f98bbafa332693cf0cc05a6e9"
    @test v110["distribution"]["archive_sha256"]=="3b3a7832cab9941d75104f7e797686c53d5d45b3bf1af33bdcbc2add6d3e0f15"
    @test v110["distribution"]["git_tree_sha1"]=="19bd3b55915f05b4e437510d1d66dc751852fb07"
    @test all(u->occursin("/surrogate-mars-global-near-surface-p20-1.1.0/",u),v110["distribution"]["urls"])
    @test v110["generation_provenance"]["runtime_evaluator_sha256"]=="1aa2c8f7244d7459208f301cdd7dd39f4f05723559ef91a4d4a0ac3941d8afb5"
    # Each version's Artifacts.toml binding matches its catalog entry.
    bindings=TOML.parsefile(PRESET_ENV._SURROGATE_ARTIFACTS)
    for e in (v100,v110)
        b=bindings[e["artifact_name"]]; d=only(b["download"])
        @test b["git-tree-sha1"]==e["distribution"]["git_tree_sha1"] && b["lazy"]===true
        @test d["sha256"]==e["distribution"]["archive_sha256"] && d["url"] in e["distribution"]["urls"]
    end
    mktempdir() do dir
        file=joinpath(dir,"fixture.jls"); catalog_file=joinpath(dir,"catalog.toml"); artifacts_file=joinpath(dir,"Artifacts.toml")
        catalog=deepcopy(catalog_template); new=deepcopy(v110); old=deepcopy(v100); catalog["presets"]=[new,old]
        for (e,v) in ((new,"1_1_0"),(old,"1_0_0"))
            e["id"]="synthetic_near_surface";e["release_enabled"]=false;e["artifact_name"]="synthetic_near_surface_"*v
        end
        # Synthetic analytic payload to the 1.1.0 contract (not Mars data): flat 0.5 km terrain, linear level states.
        payload=deepcopy(new["required_metadata"]); lev=payload["levels_km"]; g=payload["lattice"]
        nl,ni,nj=length(lev),g["nlat"],g["nlon"]
        payload["radii_km"]=(payload["generation_config"]["equatorial_radius_km"],payload["generation_config"]["polar_radius_km"])
        payload["terrain"]=Dict{String,Any}("lat0_deg"=>-86.25,"lon0_deg"=>0.488,"step_deg"=>0.5,
            "surface_height_km"=>fill(0.5,346,720),"areoid_radius_km"=>fill(3390.0,346,720))
        payload["level_T_K"]=[220.0-2lev[a] for a in 1:nl, i in 1:ni, j in 1:nj]
        payload["level_R"]=fill(191.0,nl,ni,nj); payload["level_lnp"]=[log(700.0)-lev[a]/11 for a in 1:nl, i in 1:ni, j in 1:nj]
        payload["level_source"]=ones(UInt8,nl,ni,nj); payload["surface_T30_K"]=fill(214.0,ni,nj); payload["surface_T5_K"]=fill(216.0,ni,nj)
        keys_=[(b,c) for b in -12:11 for c in 0:39]; n=length(keys_)
        payload["q_models"]=Dict{String,Any}("band"=>first.(keys_),"cell"=>last.(keys_),"L"=>fill(1,n),"order"=>fill(1,n),
            "phic_center"=>[7.5b+3.75 for (b,_) in keys_],"lam_center"=>[9.0c+4.5 for (_,c) in keys_],
            "coef"=>hcat(fill(20.0,n),zeros(n,5)),"n_points"=>fill(1,n),"status"=>fill("qualified",n))
        serialize(file,payload)
        for e in (new,old); e["payload"]["sha256"]=preset_sha(file); e["payload"]["bytes"]=filesize(file); end
        preset_write_toml(catalog_file,catalog); preset_write_toml(artifacts_file,Dict())
        opts(v)=(version=v,catalog_file=catalog_file,artifacts_file=artifacts_file,file=file,allow_unreleased=true)
        if !isdefined(GRAMSuite, :GRAMNearSurfaceAtmosphereModel)
            @test_skip "near-surface version checks need a GRAMSuite with GRAMNearSurfaceAtmosphereModel"
            return
        end
        # The 1.1.0 contract serves between the old 75 km top and 81 km, and refuses above 81 km.
        model=SpaceAGORA.surrogate_preset_model("synthetic_near_surface";opts("1.1.0")...)
        @test model isa SpaceAGORA.GRAMNearSurfaceAtmosphereModel
        s=GRAMSuite.near_surface_state(model.core,10.,40.,72000.)
        @test 75.0<s.areoid_height_km<81.0 && SpaceAGORA.getDensity(model,72000.,deg2rad(10.),deg2rad(40.),0.,true)[1]===s.density_kgm3
        @test_throws DomainError SpaceAGORA.getDensity(model,80000.,deg2rad(10.),deg2rad(40.),0.,true)
        @test SpaceAGORA.atmosphere_provenance(model)["preset_version"]=="1.1.0"
        # The same payload does not satisfy the 1.0.0 contract (levels, top and distribution version differ).
        @test_throws ArgumentError SpaceAGORA.surrogate_preset_model("synthetic_near_surface";opts("1.0.0")...)
    end
end

# A control effector that only requests a touchdown event, for the composition's synthetic descent.
mutable struct CombinedDescentTouchdown <: SpaceAGORA.AbstractControlEffectorModel
    radius_m::Float64
    clearance_m::Float64
    events::Vector{Any}
end
SpaceAGORA.SimulationModel.ControlHooks.touchdown_spec(m::CombinedDescentTouchdown, i::Int) =
    (terrain=SpaceAGORA.SimulationModel.NoTerrainModel(), reference_radius_m=m.radius_m, height_m=m.clearance_m,
     on_touchdown=(t, r_p, v_p, idx) -> push!(m.events, (t=t, r=sqrt(sum(abs2, r_p)), v=sqrt(sum(abs2, v_p)))))
SpaceAGORA.SimulationModel.ControlHooks.calcControlEffect!(::CombinedDescentTouchdown, u, p, t, i) = nothing
SpaceAGORA.SimulationModel.ControlHooks.calcControlForceTorque(::CombinedDescentTouchdown, u, p, i, t) =
    (SpaceAGORA.SimulationModel.StaticArrays.SVector(0.0, 0.0, 0.0), SpaceAGORA.SimulationModel.StaticArrays.SVector(0.0, 0.0, 0.0))

@testset "Combined descent atmosphere" begin
    callbacks=SpaceAGORA.SimulationModel.SimulationCallbacks
    # Construction rules, snapshot membership and the thread-safety trait need no GRAMSuite evaluation, so stub cores
    # check them in every environment.
    lower_stub=SpaceAGORA.GRAMNearSurfaceAtmosphereModel(:stub_lower); upper_stub=SpaceAGORA.GRAMGridAtmosphereModel(:stub_upper)
    stub=SpaceAGORA.CombinedAtmosphereModel(lower_stub,upper_stub;handover_height_m=60e3)
    @test stub isa PRESET_ENV._NativeFreeSnapshotModel && stub.lower===lower_stub && stub.upper===upper_stub && stub.handover_height_m===60000.0
    @test callbacks.density_model_threadsafe(stub)
    for bad in (NaN,Inf,true)
        @test_throws ArgumentError SpaceAGORA.CombinedAtmosphereModel(lower_stub,upper_stub;handover_height_m=bad)
    end
    for (lo,up) in ((SpaceAGORA.NoAtmosphereModel(),upper_stub),(lower_stub,SpaceAGORA.NoAtmosphereModel()))
        err=try SpaceAGORA.CombinedAtmosphereModel(lo,up;handover_height_m=60e3) catch e e end
        @test err isa ArgumentError && occursin("native-free snapshot",err.msg)
    end
    if !(isdefined(GRAMSuite,:GRAMNearSurfaceAtmosphereModel) && isdefined(GRAMSuite,:GRAMGridAtmosphereModel))
        @test_skip "combined-model queries need GRAMSuite's native-free grid and near-surface APIs"
    else
        catalog_template=TOML.parsefile(PRESET_ENV._SURROGATE_CATALOG)
        mktempdir() do dir
            # Synthetic analytic components (not Mars data), built as in the preset test sets above: a near-surface payload
            # with flat 0.5 km terrain and a flat 3390 km areoid, and a constant global grid from 40 to 260 km.
            near=only(filter(p->p["id"]=="mars_global_near_surface_p20_frozen_v1" && p["version"]=="1.0.0",catalog_template["presets"]))
            payload=deepcopy(near["required_metadata"]); lev=payload["levels_km"]; g=payload["lattice"]
            nl,ni,nj=length(lev),g["nlat"],g["nlon"]
            payload["radii_km"]=(payload["generation_config"]["equatorial_radius_km"],payload["generation_config"]["polar_radius_km"])
            payload["terrain"]=Dict{String,Any}("lat0_deg"=>-86.25,"lon0_deg"=>0.488,"step_deg"=>0.5,
                "surface_height_km"=>fill(0.5,346,720),"areoid_radius_km"=>fill(3390.0,346,720))
            payload["level_T_K"]=[220.0-2lev[a] for a in 1:nl, i in 1:ni, j in 1:nj]
            payload["level_R"]=fill(191.0,nl,ni,nj); payload["level_lnp"]=[log(700.0)-lev[a]/11 for a in 1:nl, i in 1:ni, j in 1:nj]
            payload["level_source"]=ones(UInt8,nl,ni,nj); payload["surface_T30_K"]=fill(214.0,ni,nj); payload["surface_T5_K"]=fill(216.0,ni,nj)
            keys_=[(b,c) for b in -12:11 for c in 0:39]; n=length(keys_)
            payload["q_models"]=Dict{String,Any}("band"=>first.(keys_),"cell"=>last.(keys_),"L"=>fill(1,n),"order"=>fill(1,n),
                "phic_center"=>[7.5b+3.75 for (b,_) in keys_],"lam_center"=>[9.0c+4.5 for (_,c) in keys_],
                "coef"=>hcat(fill(20.0,n),zeros(n,5)),"n_points"=>fill(1,n),"status"=>fill("qualified",n))
            lower_file=joinpath(dir,"near_surface.jls"); serialize(lower_file,payload)
            grid=deepcopy(only(filter(p->p["id"]=="odyssey_p20_frozen_v1",catalog_template["presets"]))["required_metadata"])
            grid["grid"]=Dict("alt_km"=>[40.,260.],"lat_deg"=>[-90.,90.],"lon_deg"=>[0.,180.])
            grid["fields"]=Dict(k=>fill(v,2,2,2) for (k,v) in zip(("density_kgm3","temperature_K","wind_ew_ms","wind_ns_ms","wind_up_ms"),(2e-9,150.,2.,3.,4.)))
            upper_file=joinpath(dir,"grid.jls"); serialize(upper_file,grid)
            lower=SpaceAGORA.GRAMNearSurfaceAtmosphereModel(planet="Mars",surrogate_file=lower_file)
            upper=SpaceAGORA.GRAMGridAtmosphereModel(planet="Mars",surrogate_file=upper_file)
            model=SpaceAGORA.CombinedAtmosphereModel(lower,upper;handover_height_m=60e3)
            @test model isa PRESET_ENV._NativeFreeSnapshotModel && callbacks.density_model_threadsafe(model)
            lat,lon=deg2rad(10.),deg2rad(40.)
            # Below the handover the lower component answers; at and above it, the upper one does.
            for h in (3000.,30000.,59999.)
                @test SpaceAGORA.getDensity(model,h,lat,lon,0.,true)===SpaceAGORA.getDensity(lower,h,lat,lon,0.,true)
            end
            for h in (60000.,150000.)
                @test SpaceAGORA.getDensity(model,h,lat,lon,0.,true)===SpaceAGORA.getDensity(upper,h,lat,lon,0.,true)
            end
            # Each component's own wind: none below the handover, the grid's stored wind at it.
            below=SpaceAGORA.getDensity(model,59999.,lat,lon,0.,true); at=SpaceAGORA.getDensity(model,60000.,lat,lon,0.,true)
            @test below[3]==zeros(3) && below[1]>0 && at[3]==[2.,3.,4.] && at[1]==2e-9
            @test SpaceAGORA.getDensity(model,3000.,lat,lon,1e9,false,nothing)===SpaceAGORA.getDensity(model,3000.,lat,lon,0.,true)
            # Refusals come from the selected component, with no fallback to the other one.
            @test_throws DomainError SpaceAGORA.getDensity(model,30000.,deg2rad(86.),lon,0.,true)
            @test_throws DomainError SpaceAGORA.getDensity(upper,30000.,lat,lon,0.,true)
            @test_throws DomainError SpaceAGORA.getDensity(model,261000.,lat,lon,0.,true)
            provenance=SpaceAGORA.atmosphere_provenance(model)
            @test provenance["backend"]=="combined_native_free_snapshot" && provenance["handover_height_m"]==60000.0
            @test provenance["lower"]==SpaceAGORA.atmosphere_provenance(lower) && provenance["upper"]==SpaceAGORA.atmosphere_provenance(upper)
            # Concurrent queries across the handover match serial ones.
            hs=collect(range(1000.,200000.;length=64))
            serial=[SpaceAGORA.getDensity(model,h,lat,lon,0.,true) for h in hs]
            @test fetch.([Threads.@spawn SpaceAGORA.getDensity(model,h,lat,lon,0.,true) for h in hs])==serial
            # Compatibility checks read the components' recorded metadata.
            shifted=deepcopy(grid); shifted["initial_time"]["minute"]=52; shifted["generation_config"]["initial_time"]["minute"]=52
            serialize(joinpath(dir,"grid_shifted.jls"),shifted)
            err=try SpaceAGORA.CombinedAtmosphereModel(lower,SpaceAGORA.GRAMGridAtmosphereModel(planet="Mars",surrogate_file=joinpath(dir,"grid_shifted.jls"));handover_height_m=60e3) catch e e end
            @test err isa ArgumentError && occursin("different instants",err.msg)
            other=deepcopy(grid); other["generation_config"]["equatorial_radius_km"]=3397.0
            serialize(joinpath(dir,"grid_ellipsoid.jls"),other)
            err=try SpaceAGORA.CombinedAtmosphereModel(lower,SpaceAGORA.GRAMGridAtmosphereModel(planet="Mars",surrogate_file=joinpath(dir,"grid_ellipsoid.jls"));handover_height_m=60e3) catch e e end
            @test err isa ArgumentError && occursin("reference ellipsoids",err.msg)
            for h in (30e3,260e3)
                err=try SpaceAGORA.CombinedAtmosphereModel(lower,upper;handover_height_m=h) catch e e end
                @test err isa ArgumentError && occursin("upper component's grid",err.msg)
            end
            # A synthetic descent that crosses the handover, with a touchdown 2 km above the flat synthetic terrain.
            SM=SpaceAGORA.SimulationModel; SV=SM.StaticArrays.SVector
            profile=deepcopy(grid); column(v)=[v[a] for a in 1:5, i in 1:2, j in 1:2]
            profile["grid"]=Dict("alt_km"=>[40.,60.,80.,120.,260.],"lat_deg"=>[-90.,90.],"lon_deg"=>[0.,180.])
            profile["fields"]=Dict("density_kgm3"=>column([2e-3,1.4e-4,2.5e-5,4e-7,1e-11]),"temperature_K"=>column(fill(150.,5)),
                "wind_ew_ms"=>column(fill(30.,5)),"wind_ns_ms"=>column(fill(10.,5)),"wind_up_ms"=>column(zeros(5)))
            serialize(joinpath(dir,"grid_profile.jls"),profile)
            descent_model=SpaceAGORA.CombinedAtmosphereModel(lower,
                SpaceAGORA.GRAMGridAtmosphereModel(planet="Mars",surrogate_file=joinpath(dir,"grid_profile.jls"));handover_height_m=60e3)
            planet=SM.Mars(); direction=SV(cosd(30.0),0.0,sind(30.0))   # the ellipsoid radius is about 3391.2 km here
            function descent(; density_model=descent_model, clearance=2000.0, tol=(1e-9,1e-11,0.5), flags=("0","0"))
                root=SM.Link(root=true, m=1500.0, ref_area=1.0)
                ic=SM.CartesianInitialCondition(3_491_200.0*direction, -1500.0*direction+SV(0.0,2500.0,0.0))
                sc=SM.SpacecraftModel(SM.Joint[], [root], root, true, 1500.0, 0.0, root.inertia, 0, 0, ic, 1)
                args=SM.SimulationConfiguration(
                    simulation_settings=SM.SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false),
                    mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, mission_time=3000.0, number_of_orbits=1,
                        keplerian=false, orientation_sim=false, num_steps_to_save=200),
                    environment_model=SM.EnvironmentModel(planet=planet, EI=250.0, density_model=density_model, topography=false, wind=true,
                        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
                        ephemerides_model=SM.SimpleEphemeridesModel(prime_meridian_at_reference_rad=0.0)),
                    dynamics_model=SM.DynamicsModel([sc], (SM.InverseSquaredJ2GravityModel(), SM.AerodynamicCoefficientfM())),
                    guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
                    navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
                    control_model=SM.ControlModel(control_effectors=(CombinedDescentTouchdown(3_390_500.0, clearance, Any[]),), control_rates=[1.0]),
                    initial_time=SM.InitialTime(year=2001, month=11, day=7, hour=11, minute=51, second=4.794789),
                    integration_tolerances=SM.IntegrationTolerances(reltol_orbit=tol[1], abstol_orbit=tol[2], reltol_atmosphere=tol[1],
                        abstol_atmosphere=tol[2], dt_max_orbit=tol[3], dt_max_atmosphere=tol[3]),
                    solver_config=SM.SolverConfig(solver_mode=:tsit5))
                sol=withenv("SPACEAGORA_DENSITY_FREEZE_PER_STEP"=>flags[1], "SPACEAGORA_VACUUM_GRAM_CACHE"=>flags[2]) do
                    SpaceAGORA.run_simulation(args; return_solution=true)
                end
                return sol, only(sol.prob.p.args.control_model.control_effectors).events
            end
            radius(u)=sqrt(u[1]^2+u[2]^2+u[3]^2)
            sol,events=descent()
            @test string(sol.retcode)=="Terminated" && length(events)==1 && all(!,sol.prob.p.is_active)
            @test radius(sol.u[1])>3_391_200.0+60e3 && isapprox(events[1].r,3_392_500.0;atol=1e-3)
            # Trajectory density caches and per-step freezing stay bypassed: enabling them changes nothing.
            on,events_on=descent(flags=("1","1"))
            @test on.t==sol.t && on.u==sol.u && events_on==events
            # Touchdown time and speed do not depend on the step-size controls.
            for tol in ((1e-7,1e-9,2.0),(1e-11,1e-13,0.05))
                _,e=descent(;tol)
                @test isapprox(e[1].t,events[1].t;rtol=1e-6) && isapprox(e[1].v,events[1].v;rtol=1e-6)
            end
            # The atmosphere acts on the descent: without it the spacecraft reaches the surface faster.
            _,vacuum=descent(density_model=SpaceAGORA.NoAtmosphereModel())
            @test vacuum[1].v>events[1].v+1.0
            # Below the lower component's minimum clearance the run fails with that component's error; nothing is extrapolated.
            err=try descent(clearance=1.0) catch e e end
            @test err isa DomainError && occursin("5 m",sprint(showerror,err))
            @test !any(x->occursin("libgram",lowercase(basename(x))),Libdl.dllist())
        end
    end
end
