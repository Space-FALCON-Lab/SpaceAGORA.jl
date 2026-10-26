"""Exact segment-to-point-cloud clearance, with the existing station/chaser buffer."""
function rpo_comparison_segment_clearance end

"""
Execute a retimed path using recovered nonlinear Earth two-body propagation.
The default six-axis actuator holds fixed inertial body axes and starts at the
first 0.1 s control tick, matching the full demo. Above-deadband commands use
interval-average throttle. `actuator_frame=:rtn` retains the rotating-axis
idealization for explicitly labelled diagnostics. This is a translational evaluator, not a resolved
valve pulse or attitude/wheel qualification. Fuel uses delivered force and mass
is depleted at every update. Every 0.05 s plant substep checks swept clearance.
"""
function rpo_track_two_body_path end

function rpo_execute_two_body_reference end
