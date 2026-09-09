# KS Dynamics

SpaceAGORA owns reusable Kustaanheimo–Stiefel dynamics under
`src/dynamics/ks_dynamics`, alongside its translational and rotational dynamics
components. They are independent of any particular guidance or control
implementation. Aerobraking MPC consumes the KS dynamics through this public
interface; other mission and analysis code can do the same.

## Kustaanheimo–Stiefel dynamics

Construct the planet constants, transform a Cartesian state, advance it in
fictitious time, and convert it back:

```julia
using SpaceAGORA

params = KSPropagationParams(
    Re=planet.Rp_e,
    μ=planet.μ,
    J2=planet.J2,
    Ω=norm(planet.ω),
)
state = cartesian_to_ks_state(position_ii_m, velocity_ii_m, params)
state_next = ks_rk4_step(state, params, area_m2, delta_s;
    density_kg_m3=density,
    drag_coefficient=cd,
    mass_kg=mass,
    use_drag=true,
)
cartesian = ks_state_to_cartesian(state_next)
```

The ten-element state is `[u₁,u₂,u₃,u₄,u₁′,u₂′,u₃′,u₄′,h_KS,t]`. It uses the
paper convention

```math
h_{KS}=-\varepsilon, \qquad \omega_{KS}^{2}=\frac{h_{KS}}{2},
```

where `ε` is specific orbital energy in J/kg. Consequently,
`u″=-(h_KS/2)u+G_KS a_p` and `h_KS′=-R vᵀa_p`. The helpers
`ks_energy_parameter` and `specific_energy_from_ks` make the sign convention
explicit. The independent
step `delta_s` is fictitious time; the physical elapsed time is propagated in
the final state component. Inputs use SI units.

!!! warning "Legacy KS states"
    KS states saved by the pre-change implementation stored `h_old=-2ε`. They
must not be passed directly to the current dynamics implementation. Regenerate them from
    Cartesian position and velocity, or replace component 9 with `h_old/2`.

Public functions are `ks_position`, `ks_velocity`,
`cartesian_to_ks_state`, `ks_state_to_cartesian`,
`ks_j2_acceleration_si`, `ks_drag_acceleration_si`, `ks_rhs`, and
`ks_rk4_step`. The reusable variational interface includes
`ks_kinematics_jacobians`, `ks_j2_acceleration_jacobian_si`,
`ks_rhs_jacobian`, and `ks_step_jacobian`.

The current force model supports inverse-square central gravity through the KS
energy parameter, plus optional J2 and atmospheric drag perturbations. Generic
KS state and step linearizations live with the KS dynamics. The MPC directory
retains only its drag, heat-rate, energy-output, area-input, and condensed-model
linearization because those definitions depend on the controller configuration.

## Linearization convention

For the nonlinear state `x=[p;q;h_KS;t]`, the reusable Jacobian is formed from
the paper-convention equations

```math
\begin{aligned}
p' &= q,\\
q' &= -\frac{h_{KS}}{2}p + G_{KS}(p)a_p,\\
h_{KS}' &= -R v^T a_p,\\
t' &= R.
\end{aligned}
```

In particular, `∂q′/∂h_KS=-p/2`, `∂p′/∂q=I`, and
`∂t′/∂p=2pᵀ`. `ks_rhs_jacobian` linearizes the continuous
fictitious-time equations, while `ks_step_jacobian` differentiates the complete
nonlinear RK4 step.

Writing `g(p,q,A)=G_KS(p)a_p` and `ℓ(p,q,A)=-R vᵀa_p`, the continuous
nine-state control linearization has the block structure

```math
A_c = \begin{bmatrix}
0 & I & 0\\
-\tfrac{h_{KS}}{2}I + g_p & g_q & -\tfrac12 p\\
\ell_p & \ell_q & 0
\end{bmatrix}, \qquad
B_c = \begin{bmatrix}
0\\ G_{KS}a_{p,A}\\ -R v^T a_{p,A}
\end{bmatrix}.
```

Thus changing from the legacy variable is not accomplished by changing only
the oscillator frequency. The `-p/2` column, the complete energy-rate row, and
the energy component of the area-input column must all be retained.

The MPC prediction state is `δx=[δp;δq;δh_KS]`; physical time is supplied
by the known reference horizon. Its discrete `A` and area-input `B` matrices are
computed from the same nonlinear RK4 map, including density variation at the
RK stages. The energy-output row is exactly `E=-h_KS` in J/kg, or
`E_MJ/kg=-h_KS/1e6`. This replaces the former frozen-energy eight-state
approximation.
