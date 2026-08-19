#!/usr/bin/env julia

using Dates
using SpaceAGORA_RL

const TERMINAL_LOG_FILENAME = "terminal_output.txt"

function _usage(io::IO=stdout)
    println(io, """
Usage:
  julia --project=SpaceAGORA_RL.jl \\
    SpaceAGORA_RL.jl/scripts/train.jl [CONFIG] [--wind-mode MODE]

Options:
  --wind-mode MODE   MarsGRAM wind mode: zero, nominal, or perturbed.
                     Overrides spaceagora_physics.gram_wind_mode.
  --help             Show this message.
""")
end

function _parse_cli(args)
    any(==("--help"), args) && return (help=true,)
    config_path = nothing
    wind_mode = nothing
    index = 1
    while index <= length(args)
        arg = args[index]
        if arg == "--wind-mode"
            index == length(args) && throw(ArgumentError("missing value for --wind-mode"))
            wind_mode = canonical_gram_wind_mode(args[index + 1])
            index += 2
        elseif startswith(arg, "--")
            throw(ArgumentError("unknown option: $arg"))
        elseif config_path === nothing
            config_path = arg
            index += 1
        else
            throw(ArgumentError("unexpected positional argument: $arg"))
        end
    end
    return (
        help=false,
        config_path=something(config_path, default_config_path()),
        wind_mode=wind_mode,
    )
end

function _terminal_log_timestamp()
    return Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ")
end

function _write_locked(f::Function, io_lock::ReentrantLock)
    lock(io_lock)
    try
        return f()
    finally
        unlock(io_lock)
    end
end

function _tee_stream!(reader::IO, terminal::IO, log_io::IO, lock::ReentrantLock)
    try
        while true
            bytes = readavailable(reader)
            isempty(bytes) && break
            _write_locked(lock) do
                write(terminal, bytes)
                flush(terminal)
                write(log_io, bytes)
                flush(log_io)
            end
        end
    catch
        # The pipe is expected to close during shutdown; do not let tee cleanup
        # mask the original training error.
    finally
        close(reader)
    end
    return nothing
end

function _with_terminal_log(f::Function, output_dir::AbstractString)
    mkpath(output_dir)
    log_path = joinpath(output_dir, TERMINAL_LOG_FILENAME)
    open(log_path, "a") do log_io
        log_lock = ReentrantLock()
        original_stdout = stdout
        original_stderr = stderr
        stdout_reader, stdout_writer = redirect_stdout()
        stderr_reader, stderr_writer = redirect_stderr()
        stdout_task = @async _tee_stream!(stdout_reader, original_stdout, log_io, log_lock)
        stderr_task = @async _tee_stream!(stderr_reader, original_stderr, log_io, log_lock)
        try
            println("saving terminal output to ", log_path)
            println("=== terminal log started ", _terminal_log_timestamp(), " ===")
            result = f(log_path)
            println("=== terminal log completed ", _terminal_log_timestamp(), " ===")
            return result
        catch err
            bt = catch_backtrace()
            _write_locked(log_lock) do
                println(log_io)
                println(log_io, "=== terminal log captured error ", _terminal_log_timestamp(), " ===")
                showerror(log_io, err, bt)
                println(log_io)
                flush(log_io)
            end
            rethrow()
        finally
            flush(stdout)
            flush(stderr)
            redirect_stdout(original_stdout)
            redirect_stderr(original_stderr)
            close(stdout_writer)
            close(stderr_writer)
            wait(stdout_task)
            wait(stderr_task)
        end
    end
end

function _train_session(session)
    if session.learner.device isa CUDATrainingDevice
        return Base.invokelatest(train_parallel!, session)
    end
    return train_parallel!(session)
end

function _validate_checkpoints(session)
    config = session.config
    config.training.validate_checkpoints || return nothing
    config.training.algorithm in (:pr_drl, :ddqn, :a2c) || return nothing
    output_dir = joinpath(session.output_dir, "checkpoint_validation")
    return validate_frozen_checkpoints(
        session.output_dir,
        config.scenario;
        episodes=config.training.validation_episodes,
        seed=config.training.validation_seed,
        output_dir=output_dir,
        protected_initialization=protected_initialization_config(config.training),
        checkpoint_stride=config.training.validation_checkpoint_stride,
        write_plots=config.reports.write_plots,
    )
end

function main(args=ARGS)
    options = _parse_cli(args)
    if options.help
        _usage()
        return nothing
    end
    config = resolve_config(options.config_path; gram_wind_mode=options.wind_mode)
    session = build_training_session(config)
    return _with_terminal_log(session.output_dir) do _
        println("MarsGRAM wind mode: ", config.scenario.spaceagora_gram_wind_mode)
        result = _train_session(session)
        validation = _validate_checkpoints(session)
        println("wrote run artifacts to ", result.output_dir)
        println("episodes: ", length(result.summaries),
                " transitions: ", length(result.transitions),
                " global_step: ", result.global_step)
        return merge(result, (validation=validation,))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main()
    catch error
        showerror(stderr, error)
        println(stderr)
        println(stderr, "Use --help for usage.")
        exit(1)
    end
end
