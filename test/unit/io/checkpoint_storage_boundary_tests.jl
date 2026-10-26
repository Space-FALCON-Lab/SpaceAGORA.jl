using Test
using Serialization
using TOML
using SHA

module CheckpointStorageBoundaryFixture
    const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
    include(joinpath(ROOT, "src", "interfaces", "io", "paths", "io_paths.jl"))
    include(joinpath(ROOT, "src", "interfaces", "io", "checkpoints", "checkpoint_storage.jl"))
end

@testset "I/O paths and checkpoint storage retain their schema" begin
    paths = CheckpointStorageBoundaryFixture.IOConfig
    storage = CheckpointStorageBoundaryFixture.IOSerialization
    @test nameof(paths) == :IOConfig
    @test nameof(storage) == :IOSerialization
    @test parentmodule(paths._checkpoint_paths) === paths
    @test parentmodule(storage._write_checkpoint!) === storage
    mktempdir() do root
        config = (simulation_settings=(results_directory=joinpath(root, "results"), checkpoint_directory=""),)
        p = paths._checkpoint_paths(config)
        @test p.data == joinpath(root, "results", "checkpoints", "simulation_checkpoint.bin")
        @test p.manifest == joinpath(root, "results", "checkpoints", "simulation_checkpoint.manifest.toml")
        @test paths._results_bundle_prefix(config) == joinpath(root, "results", "simulation_results")
        @test paths._results_csv_path(config) == joinpath(root, "results", "simulation_results.csv")
        explicit = (simulation_settings=(results_directory=root, checkpoint_directory=joinpath(root, "custom")),)
        @test paths._checkpoint_paths(explicit).data == joinpath(root, "custom", "simulation_checkpoint.bin")
        @test storage._load_checkpoint(config) === nothing
        state = [1.0, -2.5, 3.0]
        @test storage._write_checkpoint!(config, 12.5, state, 7; solver_mode="tsit5") === nothing
        state[1] = 99.0
        payload = deserialize(p.data)
        @test keys(payload) == (:schema_version, :created_utc, :t, :solver_mode, :u)
        @test payload.schema_version == 7
        @test payload.t == 12.5
        @test payload.u == [1.0, -2.5, 3.0]
        @test payload.solver_mode == "tsit5"
        manifest = TOML.parsefile(p.manifest)
        @test Set(keys(manifest)) == Set(["schema_version", "created_utc", "time_s", "solver_mode", "data_path", "data_size_bytes", "data_sha256"])
        @test manifest["schema_version"] == 7
        @test manifest["time_s"] == 12.5
        @test manifest["solver_mode"] == "tsit5"
        @test manifest["data_path"] == p.data
        @test manifest["data_size_bytes"] == filesize(p.data)
        @test manifest["data_sha256"] == bytes2hex(sha256(read(p.data)))
        restored = storage._load_checkpoint(config)
        @test restored == (t=12.5, u=[1.0, -2.5, 3.0], solver_mode="tsit5", data_path=p.data, manifest_path=p.manifest)
        serialize(p.data, (t=1, u=[2.0]))
        @test storage._load_checkpoint(config).solver_mode === nothing
        serialize(p.data, (t=1,))
        @test_throws ArgumentError storage._load_checkpoint(config)
        write(p.interlinks, "active history")
        archived = p.interlinks * ".before_resume"
        write(archived, "diagnostic history")
        @test storage._clear_checkpoint!(config) === nothing
        @test !ispath(p.interlinks)
        @test read(archived, String) == "diagnostic history"
        @test storage._clear_checkpoint!(config) === nothing
        @test !ispath(p.data) && !ispath(p.manifest)
        @test storage._load_checkpoint(config) === nothing

        target = joinpath(root, "atomic.txt")
        write(target, "original")
        before = sort(readdir(root))
        writer = tmp -> begin
            write(tmp, "partial")
            error("fixture writer failure")
        end
        @test_throws ErrorException storage._atomic_write_file(target, writer)
        @test read(target, String) == "original"
        @test sort(readdir(root)) == before
    end
end

@testset "Checkpoint history hashes only the requested byte range" begin
    storage = CheckpointStorageBoundaryFixture.IOSerialization
    mktempdir() do root
        path = joinpath(root, "history.bin")
        bytes = UInt8[mod(i, 251) for i in 1:200_000]
        write(path, bytes)
        context = SHA256_CTX()
        storage._hash_checkpoint_rows!(context, path, 0, 70_001)
        @test digest!(deepcopy(context)) == sha256(bytes[1:70_001])
        # Mutating the already-hashed prefix demonstrates that an append hashes
        # only the new range; the saved context retains the original prefix.
        open(path, "r+") do io
            write(io, UInt8(0xff))
        end
        storage._hash_checkpoint_rows!(context, path, 70_001, length(bytes))
        @test digest!(deepcopy(context)) == sha256(bytes)
        @test digest!(deepcopy(context)) != sha256(read(path))
        storage._hash_checkpoint_rows!(context, path, length(bytes), length(bytes))
        @test digest!(deepcopy(context)) == sha256(bytes)
        @test_throws ArgumentError storage._hash_checkpoint_rows!(context, path, -1, 20)
        @test_throws ArgumentError storage._hash_checkpoint_rows!(context, path, 20, 19)
        @test_throws EOFError storage._hash_checkpoint_rows!(context, path, 0, length(bytes) + 1)
    end
end
