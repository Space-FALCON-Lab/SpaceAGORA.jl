# Plots from saved Feather data

Source: `/home/space-falcon-3/Desktop/Summary_Table_for_Other_Labs/1_Kuang's Prototype Code/output/feather/prototype_N20_h1000km_t1000km_ih0_it0_e0_nu1_T315355.97s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s`. No orbit simulation was rerun.

Reuses the target-only result plots from test16_options.jl; target spacecraft is 21. Orbital-element differences retain the reference's spacecraft 1 versus 2 comparison.

Positions, velocities and masses come from saved samples. Laser force, acceleration, delta-v and impulse are reconstructed with the prototype diagnostic routines: power=10000.0 W, magnification=100.0, schedule=gve_sma, J2=true, drag=false. These parameters are supplied explicitly because the Feather metadata does not contain all of them.

Maximum reconstructed versus stored target RTN delta-v difference: 0.0 m/s. Sampling is 10.0 seconds plus the endpoint; plots do not restore unsaved solver steps.
