@testset "Aerodynamic Helper Branch Coverage" begin
    dynamic_effectors = SimulationModel.DynamicEffectors

    @test dynamic_effectors._parse_bool_env("SPACEAGORA_TEST_BOOL_PARSE", false) == false
    withenv("SPACEAGORA_TEST_BOOL_PARSE" => "yes") do
        @test dynamic_effectors._parse_bool_env("SPACEAGORA_TEST_BOOL_PARSE", false) == true
    end
    withenv("SPACEAGORA_TEST_BOOL_PARSE" => "off") do
        @test dynamic_effectors._parse_bool_env("SPACEAGORA_TEST_BOOL_PARSE", true) == false
    end
    withenv("SPACEAGORA_TEST_BOOL_PARSE" => "invalid") do
        @test_throws ArgumentError dynamic_effectors._parse_bool_env("SPACEAGORA_TEST_BOOL_PARSE", false)
    end

    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "off") do
        @test dynamic_effectors._multibody_parallel_mode() == :off
    end
    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "on") do
        @test dynamic_effectors._multibody_parallel_mode() == :on
    end
    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "auto") do
        @test dynamic_effectors._multibody_parallel_mode() == :auto
    end
    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "invalid") do
        @test_throws ArgumentError dynamic_effectors._multibody_parallel_mode()
    end

    withenv("SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "4") do
        @test dynamic_effectors._multibody_thread_threshold() == 4
    end
    withenv("SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "0") do
        @test dynamic_effectors._multibody_thread_threshold() == 1
    end
    withenv("SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "oops") do
        @test_throws ArgumentError dynamic_effectors._multibody_thread_threshold()
    end
    withenv("SPACEAGORA_MULTIBODY_MAX_THREADS" => "3") do
        @test dynamic_effectors._multibody_max_threads() == 3
    end
    withenv("SPACEAGORA_MULTIBODY_MAX_THREADS" => "oops") do
        @test_throws ArgumentError dynamic_effectors._multibody_max_threads()
    end

    withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1") do
        @test dynamic_effectors._multibody_outer_parallel_hint() == true
    end
    withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0") do
        @test dynamic_effectors._multibody_outer_parallel_hint() == false
    end
    withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "oops") do
        @test_throws ArgumentError dynamic_effectors._multibody_outer_parallel_hint()
    end

    @test dynamic_effectors._multibody_use_threads(1) == false
    if Threads.nthreads() > 1
        withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "on") do
            @test dynamic_effectors._multibody_use_threads(64) == true
        end
        withenv(
            "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
            "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "2",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
            "SPACEAGORA_MULTIBODY_PARALLEL_ALLOW_WITH_OUTER" => "0"
        ) do
            @test dynamic_effectors._multibody_use_threads(64) == false
        end
        withenv(
            "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
            "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "2",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
            "SPACEAGORA_MULTIBODY_PARALLEL_HEAVY_ONLY" => "1"
        ) do
            @test dynamic_effectors._multibody_use_threads(64; heavy_work=false) == false
        end
    end

    @test dynamic_effectors._threadid_capacity() >= Threads.maxthreadid()

    body_a = Link{0}(root=true)
    body_b = Link{0}(root=false)
    body_a.net_force .= SVector{3, Float64}(1.0, 2.0, 3.0)
    body_a.net_torque .= SVector{3, Float64}(4.0, 5.0, 6.0)
    body_b.net_force .= SVector{3, Float64}(-0.5, 0.0, 0.5)
    body_b.net_torque .= SVector{3, Float64}(1.0, -1.0, 0.0)
    force_sum, torque_sum = dynamic_effectors.collect_and_reset_link_wrenches!([body_a, body_b])

    @test force_sum == SVector{3, Float64}(0.5, 2.0, 3.5)
    @test torque_sum == SVector{3, Float64}(5.0, 4.0, 6.0)
    @test body_a.net_force == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test body_a.net_torque == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test body_b.net_force == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test body_b.net_torque == SVector{3, Float64}(0.0, 0.0, 0.0)
end




