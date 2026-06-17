#!/usr/bin/env julia

using CSV
using DataFrames
using Printf
using Plots
using SpaceAGORA_RL

function _example_output_dir()
    return joinpath(dirname(dirname(@__DIR__)), "outputs", "examples", "mars_odyssey_rl")
end

function _plot_pass_logs(pass_logs::DataFrame, output_dir::AbstractString)
    mkpath(output_dir)
    grouped = groupby(pass_logs, [:policy, :episode])

    p1 = plot(xlabel="Pass", ylabel="Apoapsis radius (km)", legend=:outerright)
    p2 = plot(xlabel="Pass", ylabel="Heat rate (W/cm^2)", legend=:outerright)
    p3 = plot(xlabel="Pass", ylabel="Action delta-V (m/s)", legend=:outerright)

    for group in grouped
        label = string(group.policy[1], " ep", group.episode[1])
        plot!(p1, group.pass, group.apoapsis_radius_km; label=label)
        plot!(p2, group.pass, group.heat_rate_w_cm2; label=label)
        plot!(p3, group.pass, group.action_delta_v_mps; label=label, seriestype=:steppre)
    end

    hline!(p2, [0.05, 0.25]; label=["low corridor" "high corridor"], linestyle=:dash, color=[:gray :red])
    savefig(p1, joinpath(output_dir, "apoapsis_history.png"))
    savefig(p2, joinpath(output_dir, "heat_rate_history.png"))
    savefig(p3, joinpath(output_dir, "action_history.png"))
end

function _plot_episode_metrics(metrics::DataFrame, output_dir::AbstractString)
    labels = string.(metrics.policy, " ep", metrics.episode)
    p = bar(labels, abs.(metrics.target_error_km);
            xlabel="Policy",
            ylabel="Final target error (km)",
            legend=false,
            xrotation=20)
    savefig(p, joinpath(output_dir, "target_error_by_policy.png"))
end

function _js_number_array(values)
    return "[" * join((@sprintf("%.6f", Float64(value)) for value in values), ",") * "]"
end

function _js_string(value)
    escaped = replace(String(value), "\\" => "\\\\", "\"" => "\\\"")
    return "\"" * escaped * "\""
end

function _perifocal_to_inertial(xp, yp, omega, raan, inc)
    co = cos(omega)
    so = sin(omega)
    cO = cos(raan)
    sO = sin(raan)
    ci = cos(inc)
    si = sin(inc)

    x = (cO * co - sO * so * ci) * xp + (-cO * so - sO * co * ci) * yp
    y = (sO * co + cO * so * ci) * xp + (-sO * so + cO * co * ci) * yp
    z = (so * si) * xp + (co * si) * yp
    return x, y, z
end

function _orbit_series(pass_logs::DataFrame, planet_radius_km::Real; samples_per_pass::Int=96)
    series = String[]
    for group in groupby(pass_logs, [:policy, :episode])
        x = Float64[]
        y = Float64[]
        z = Float64[]
        pass_axis = Float64[]
        for row in eachrow(group)
            ra = Float64(row.apoapsis_radius_km)
            rp = Float64(planet_radius_km + row.periapsis_altitude_km)
            a = (ra + rp) / 2
            e = clamp((ra - rp) / (ra + rp), 0.0, 0.999)
            p = a * (1 - e^2)
            for sample in 0:(samples_per_pass - 1)
                theta = 2pi * sample / samples_per_pass
                r = p / (1 + e * cos(theta))
                xi, yi, zi = _perifocal_to_inertial(
                    r * cos(theta),
                    r * sin(theta),
                    Float64(row.argument_of_periapsis_rad),
                    Float64(row.raan_rad),
                    Float64(row.inclination_rad),
                )
                push!(x, xi)
                push!(y, yi)
                push!(z, zi)
                push!(pass_axis, Float64(row.pass) + sample / samples_per_pass)
            end
        end
        label = string(group.policy[1], " episode ", group.episode[1])
        push!(series, string(
            "{label:", _js_string(label),
            ",x:", _js_number_array(x),
            ",y:", _js_number_array(y),
            ",z:", _js_number_array(z),
            ",pass:", _js_number_array(pass_axis), "}",
        ))
    end
    return "[" * join(series, ",") * "]"
end

function _write_orbit_html(pass_logs::DataFrame, output_dir::AbstractString, config)
    planet_radius_km = config.mars_radius_m / 1000
    series_json = _orbit_series(pass_logs, planet_radius_km)
    html_path = joinpath(output_dir, "earth_orbit_trajectory_3d.html")
    html = """
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>3D Earth Orbit Trajectory</title>
  <style>
    :root { color-scheme: dark; font-family: Inter, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }
    body { margin: 0; background: #08111f; color: #eef5ff; overflow: hidden; }
    #stage { width: 100vw; height: 100vh; display: block; background: radial-gradient(circle at 50% 42%, #193653 0%, #08111f 55%, #030712 100%); }
    #panel { position: fixed; left: 18px; right: 18px; bottom: 16px; display: grid; grid-template-columns: minmax(170px, 260px) 90px 1fr 120px; gap: 12px; align-items: center; padding: 12px; border: 1px solid rgba(180,210,255,0.2); background: rgba(4,10,22,0.78); backdrop-filter: blur(10px); border-radius: 8px; }
    select, button, input { width: 100%; }
    select, button { color: #eef5ff; background: #10223b; border: 1px solid rgba(180,210,255,0.28); border-radius: 6px; height: 34px; }
    input[type=range] { accent-color: #f7a84b; }
    #readout { font-size: 13px; color: #b8c8dd; text-align: right; }
    #title { position: fixed; left: 20px; top: 16px; font-size: 17px; letter-spacing: 0; color: #f7fbff; }
    #hint { position: fixed; right: 20px; top: 18px; font-size: 12px; color: #a8bad3; }
  </style>
</head>
<body>
  <canvas id="stage"></canvas>
  <div id="title">3D Earth-style orbit trajectory</div>
  <div id="hint">Drag to rotate, scroll to zoom</div>
  <div id="panel">
    <select id="series"></select>
    <button id="play">Play</button>
    <input id="time" type="range" min="0" max="0" value="0" step="1">
    <div id="readout"></div>
  </div>
  <script>
    const PLANET_RADIUS = $(@sprintf("%.6f", planet_radius_km));
    const SERIES = $series_json;
    const canvas = document.getElementById("stage");
    const ctx = canvas.getContext("2d");
    const selector = document.getElementById("series");
    const slider = document.getElementById("time");
    const playButton = document.getElementById("play");
    const readout = document.getElementById("readout");
    let yaw = -0.72;
    let pitch = 0.42;
    let zoom = 1.0;
    let playing = false;
    let drag = null;

    SERIES.forEach((series, index) => {
      const option = document.createElement("option");
      option.value = index;
      option.textContent = series.label;
      selector.appendChild(option);
    });

    function currentSeries() { return SERIES[Number(selector.value || 0)] || SERIES[0]; }
    function resize() {
      const dpr = window.devicePixelRatio || 1;
      canvas.width = Math.round(innerWidth * dpr);
      canvas.height = Math.round(innerHeight * dpr);
      canvas.style.width = innerWidth + "px";
      canvas.style.height = innerHeight + "px";
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      draw();
    }

    function setSliderBounds() {
      const series = currentSeries();
      slider.max = Math.max(0, series.x.length - 1);
      slider.value = 0;
      draw();
    }

    function rotatePoint(x, y, z) {
      const cy = Math.cos(yaw), sy = Math.sin(yaw);
      const cp = Math.cos(pitch), sp = Math.sin(pitch);
      const x1 = cy * x + sy * z;
      const z1 = -sy * x + cy * z;
      const y1 = cp * y - sp * z1;
      const z2 = sp * y + cp * z1;
      return [x1, y1, z2];
    }

    function project(x, y, z, scale, cx, cy) {
      const p = rotatePoint(x, y, z);
      return [cx + p[0] * scale, cy - p[1] * scale, p[2]];
    }

    function drawGrid(scale, cx, cy) {
      ctx.lineWidth = 1;
      for (let lat = -60; lat <= 60; lat += 30) {
        ctx.beginPath();
        let started = false;
        for (let lon = 0; lon <= 360; lon += 4) {
          const phi = lat * Math.PI / 180;
          const lam = lon * Math.PI / 180;
          const r = PLANET_RADIUS * Math.cos(phi);
          const x = r * Math.cos(lam);
          const y = PLANET_RADIUS * Math.sin(phi);
          const z = r * Math.sin(lam);
          const pt = project(x, y, z, scale, cx, cy);
          if (pt[2] < 0) { started = false; continue; }
          if (!started) { ctx.moveTo(pt[0], pt[1]); started = true; } else { ctx.lineTo(pt[0], pt[1]); }
        }
        ctx.strokeStyle = "rgba(202,230,255,0.15)";
        ctx.stroke();
      }
      for (let lon = 0; lon < 180; lon += 30) {
        ctx.beginPath();
        let started = false;
        for (let lat = -90; lat <= 90; lat += 4) {
          const phi = lat * Math.PI / 180;
          const lam = lon * Math.PI / 180;
          const r = PLANET_RADIUS * Math.cos(phi);
          const x = r * Math.cos(lam);
          const y = PLANET_RADIUS * Math.sin(phi);
          const z = r * Math.sin(lam);
          const pt = project(x, y, z, scale, cx, cy);
          if (pt[2] < 0) { started = false; continue; }
          if (!started) { ctx.moveTo(pt[0], pt[1]); started = true; } else { ctx.lineTo(pt[0], pt[1]); }
        }
        ctx.strokeStyle = "rgba(202,230,255,0.12)";
        ctx.stroke();
      }
    }

    function drawPath(series, endIndex, scale, cx, cy) {
      function strokeRange(start, stop, color, width) {
        ctx.beginPath();
        let started = false;
        for (let i = start; i <= stop; i++) {
          const pt = project(series.x[i], series.y[i], series.z[i], scale, cx, cy);
          if (!started) { ctx.moveTo(pt[0], pt[1]); started = true; } else { ctx.lineTo(pt[0], pt[1]); }
        }
        ctx.strokeStyle = color;
        ctx.lineWidth = width;
        ctx.stroke();
      }
      if (series.x.length < 2) return;
      strokeRange(0, series.x.length - 1, "rgba(124,188,255,0.22)", 1.2);
      strokeRange(0, Math.max(1, endIndex), "rgba(255,177,84,0.95)", 2.6);
      const pt = project(series.x[endIndex], series.y[endIndex], series.z[endIndex], scale, cx, cy);
      ctx.beginPath();
      ctx.arc(pt[0], pt[1], 5.5, 0, 2 * Math.PI);
      ctx.fillStyle = "#ffffff";
      ctx.fill();
      ctx.strokeStyle = "#f7a84b";
      ctx.lineWidth = 2;
      ctx.stroke();
    }

    function drawPlanet(scale, cx, cy) {
      const radiusPx = PLANET_RADIUS * scale;
      const gradient = ctx.createRadialGradient(cx - radiusPx * 0.35, cy - radiusPx * 0.42, radiusPx * 0.12, cx, cy, radiusPx);
      gradient.addColorStop(0, "#9de3ff");
      gradient.addColorStop(0.42, "#2879d0");
      gradient.addColorStop(1, "#062a5c");
      ctx.beginPath();
      ctx.arc(cx, cy, radiusPx, 0, 2 * Math.PI);
      ctx.fillStyle = gradient;
      ctx.fill();
      ctx.strokeStyle = "rgba(215,238,255,0.45)";
      ctx.lineWidth = 1.4;
      ctx.stroke();
      drawGrid(scale, cx, cy);
    }

    function draw() {
      ctx.clearRect(0, 0, innerWidth, innerHeight);
      const series = currentSeries();
      if (!series) return;
      const maxRadius = Math.max(...series.x.map((x, i) => Math.hypot(x, series.y[i], series.z[i])), PLANET_RADIUS);
      const scale = 0.42 * Math.min(innerWidth, innerHeight) / maxRadius * zoom;
      const cx = innerWidth / 2;
      const cy = innerHeight / 2 - 18;
      const endIndex = Math.min(Number(slider.value || 0), series.x.length - 1);
      drawPlanet(scale, cx, cy);
      drawPath(series, endIndex, scale, cx, cy);
      const passValue = series.pass[endIndex] || 0;
      readout.textContent = "pass " + passValue.toFixed(2) + " / " + Math.ceil(series.pass[series.pass.length - 1]);
    }

    selector.addEventListener("change", setSliderBounds);
    slider.addEventListener("input", draw);
    playButton.addEventListener("click", () => { playing = !playing; playButton.textContent = playing ? "Pause" : "Play"; });
    canvas.addEventListener("mousedown", event => { drag = {x: event.clientX, y: event.clientY, yaw, pitch}; });
    window.addEventListener("mouseup", () => { drag = null; });
    window.addEventListener("mousemove", event => {
      if (!drag) return;
      yaw = drag.yaw + (event.clientX - drag.x) * 0.006;
      pitch = Math.max(-1.35, Math.min(1.35, drag.pitch + (event.clientY - drag.y) * 0.006));
      draw();
    });
    canvas.addEventListener("wheel", event => {
      event.preventDefault();
      zoom = Math.max(0.45, Math.min(2.8, zoom * Math.exp(-event.deltaY * 0.001)));
      draw();
    }, {passive: false});
    window.addEventListener("resize", resize);
    function tick() {
      if (playing) {
        const next = Number(slider.value) + 2;
        slider.value = next > Number(slider.max) ? 0 : next;
        draw();
      }
      requestAnimationFrame(tick);
    }
    setSliderBounds();
    resize();
    tick();
  </script>
</body>
</html>
"""
    open(html_path, "w") do io
        write(io, html)
    end
    return html_path
end

function _comparison_results(config, episodes::Int, seed::Int,
                             checkpoint_path::Union{Nothing,AbstractString};
                             include_all_baselines::Bool)
    if include_all_baselines
        results = evaluate_baselines(config; episodes=episodes, seed=seed)
    else
        results = Dict(
            "aads_heuristic" => evaluate_policy(AADSHeuristicPolicy(), config;
                                                episodes=episodes,
                                                seed=seed,
                                                policy_name="aads_heuristic"),
        )
    end

    if checkpoint_path !== nothing
        policy = load_trained_ddqn_policy(checkpoint_path)
        results["trained_ddqn"] = evaluate_policy(policy, config;
                                                  episodes=episodes,
                                                  seed=seed,
                                                  policy_name="trained_ddqn")
    end
    return results
end

function run_example(; episodes::Int=3,
                     seed::Int=7,
                     checkpoint_path::Union{Nothing,AbstractString}=nothing,
                     include_all_baselines::Bool=true,
                     output_dir::AbstractString=_example_output_dir())
    config = default_aerobraking_config(phase="Main", nominal=true, max_passes=80, training=false)
    results = _comparison_results(config, episodes, seed, checkpoint_path;
                                  include_all_baselines=include_all_baselines)
    paths = write_evaluation_artifacts(output_dir, results)

    metrics = DataFrame(CSV.File(paths.metrics))
    pass_logs = DataFrame(CSV.File(paths.pass_logs))
    _plot_pass_logs(pass_logs, output_dir)
    _plot_episode_metrics(metrics, output_dir)
    orbit_html = _write_orbit_html(pass_logs, output_dir, config)

    println("wrote CSV and plot artifacts to ", output_dir)
    println("wrote interactive 3D orbit HTML to ", orbit_html)
    return (results=results, paths=merge(paths, (; orbit_html=orbit_html)), output_dir=output_dir)
end

function run_trained_comparison(checkpoint_path::AbstractString;
                                episodes::Int=40,
                                seed::Int=42,
                                output_dir::AbstractString=joinpath(_example_output_dir(), "trained_vs_heuristic"))
    return run_example(episodes=episodes,
                       seed=seed,
                       checkpoint_path=checkpoint_path,
                       include_all_baselines=false,
                       output_dir=output_dir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_example()
end
