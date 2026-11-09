# Animations of the six backlog cases

[Open the video gallery](animations.html). These are the original **satellite-to-satellite** experiments, not the later forty-spacecraft ground-station run. No simulation was rerun.

| Case | MP4 | Preview | Playback (s) | Recorded set changes | Delivered (Mbit) | Overflow (Mbit) |
|---|---|---|---:|---:|---:|---:|
| producer_control | [Video](results/backlog/producer_control/animation/animation.mp4) | [Initial](results/backlog/producer_control/animation/preview.png) / [Final](results/backlog/producer_control/animation/final.png) | 23.467 | 0 | 4.000000 | 0.000000 |
| production_competition | [Video](results/backlog/production_competition/animation/animation.mp4) | [Initial](results/backlog/production_competition/animation/preview.png) / [Final](results/backlog/production_competition/animation/final.png) | 27.967 | 9 | 11.082066 | 0.000000 |
| busy_relay | [Video](results/backlog/busy_relay/animation/animation.mp4) | [Initial](results/backlog/busy_relay/animation/preview.png) / [Final](results/backlog/busy_relay/animation/final.png) | 23.967 | 1 | 20.000061 | 19.999879 |
| direct_route_control | [Video](results/backlog/direct_route_control/animation/animation.mp4) | [Initial](results/backlog/direct_route_control/animation/preview.png) / [Final](results/backlog/direct_route_control/animation/final.png) | 23.467 | 0 | 30.000000 | 0.000000 |
| production_competition_refined | [Video](results/backlog/production_competition_refined/animation/animation.mp4) | [Initial](results/backlog/production_competition_refined/animation/preview.png) / [Final](results/backlog/production_competition_refined/animation/final.png) | 32.967 | 19 | 11.082066 | 0.000000 |
| busy_relay_refined | [Video](results/backlog/busy_relay_refined/animation/animation.mp4) | [Initial](results/backlog/busy_relay_refined/animation/preview.png) / [Final](results/backlog/busy_relay_refined/animation/final.png) | 23.967 | 1 | 20.000061 | 19.999879 |

## Method and interpretation

Adapted from [the reference GLMakie renderer](<../../Summary_Table_for_Other_Labs/3_Kuang's Prototype Code_for_animation/functions/10_Animation_ver2.jl>) and its recorded-data replay approach: reactive `Observable` frame updates, `Axis3` orbital context, `GLMakie.record` MP4 export, PNG snapshots, and headless rendering under `xvfb-run`. The existing Earth texture and rendering environment are reused read-only. The reference renderer itself is not modified.

The original renderer is laser/helper/target-specific. This adaptation uses directed **data** links, actual terminal IDs, storage bars, and delivery/production/overflow counters. The right-hand network is explicitly schematic because the real spacecraft are only 10-20 km apart and overlap visually in a full Earth view. The dashed orbital circle is a reference orbit; colored trails use the saved 20-second trajectories.

- Every original selected-set change is represented, including the busy relay's initial direct link lasting only about 0.000121 s.
- A labelled review hold pauses simulation time for 1.5 s initially, 0.5 s at each selection change, and 2 s at final totals. These are **presentation pauses**, not longer physical contacts. Playback duration therefore differs between cases.
- The base display grid is 30 frames per simulated second, with exact change times inserted. Positions and cumulative counters use linear interpolation of saved values. The CSVs do not contain velocities, so the reference's cubic Hermite interpolation is not claimed. Selections, availability and scores are held from the latest accepted scheduler sample; they are never interpolated or recomputed.
- Orange means selected; dashed gray means available but idle; dotted light gray means excluded or unavailable. Arrows indicate direction, not individual packets or forces. The label's maximum rate and backlog score are **not** claims of actual carried rate.
- All numerical source CSVs remain unchanged. The per-video `frames.csv` maps media frames to simulation time, source samples, selections and displayed counters. `validation.toml` includes input/code/video hashes and codec checks.

## Validation

Input states are checked against the accepted-step history and link windows. Every encoded frame checks the rendered positions, connection visibility and data curves. FFprobe verifies H.264, yuv420p, 1600x1080, 30 fps, frame count and duration; FFmpeg decodes each complete video with errors treated as failures. See [animation_index.csv](results/backlog/animation_index.csv) and each case's validation record.

## Reproduce

From the repository root, using the existing reference rendering environment:

```sh
xvfb-run -a julia --startup-file=no \
  --project="../Summary_Table_for_Other_Labs/2_SpaceAGORA.jl" \
  kuang_later_tasks/animate_backlog.jl
```

Append one or more of the six case names to render only those cases. A full run refreshes this index and the gallery; partial runs update only the requested case outputs.

