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
        catalog=deepcopy(catalog_template); entry=only(catalog["presets"])
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
