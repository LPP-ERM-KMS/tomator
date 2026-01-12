# tomator_dolfinx

1D cylindrical plasma transport solver using FEniCS/dolfinx.

## Governing Equation

The solver solves the cylindrical transport equation for each species:

$$\frac{\partial n}{\partial t} = \frac{1}{r} \frac{\partial}{\partial r}\left(r D \frac{\partial n}{\partial r}\right) - \frac{1}{r} \frac{\partial}{\partial r}(r V n) + S$$

where:
- $n$ = species density [m⁻³]
- $D$ = diffusion coefficient [m²/s]
- $V$ = advection velocity [m/s]
- $S$ = source/sink terms [m⁻³s⁻¹]
- $r$ = radial coordinate [m]

## Physics Steps (per timestep)

The solver uses **operator splitting** with implicit parallel losses:

1. **Collision Sources** — Compute ionization, recombination, charge exchange rates (ADAS data) and collision frequencies (ν)
2. **RF Power** — Compute electron heating from IC/EC waves
3. **Collision Timestep** — Compute per-cell dt_collision from source magnitudes (ignoring vacuum regions where n < nevac)
4. **Parallel Loss Rates** — Compute k_n, k_E loss rates for implicit LHS treatment (limiter + Bpol losses)
5. **Transport Coefficients** — Calculate D (Bohm/gyro-geometric) and V for each species
6. **Diffusion Timestep** — Compute per-cell dt_diffusion from stability criterion
7. **Optimal Timestep** — Combine dt_collision and dt_diffusion; first step always uses dt_init (ignores neutrals)
8. **BDF2 Coefficients** — Update time discretization weights (variable for adaptive dt)
9. **Transport Solve** — Advance densities and energies via BDF2 FEM with implicit parallel losses:
   - Ion transport (n, E) with implicit k_n, k_E on LHS
   - Neutral transport (n, E) — no parallel losses
   - Electron energy transport with implicit k_E on LHS
10. **Reaction Solve** — Implicit Newton solve for all reactions at each mesh point
11. **Post-processing** — Apply density floors, quasi-neutrality, temperature clamps
12. **Timestep Adaptation** — Adjust dt based on solution change (grows toward dt_min with maxtstepincrement)

## Implicit Parallel Losses

Parallel losses (limiter and Bpol) are treated **implicitly** by adding destruction rate terms to the transport equation LHS:

$$c_1 n^{k+1} + k_n \Delta t \cdot n^{k+1} = c_2 n^k - c_3 n^{k-1} + \Delta t \cdot \text{RHS}$$

This allows much larger timesteps than explicit treatment since parallel losses don't restrict dt_collision.

## Timestep Control

- **dt_init**: Initial timestep (always used for first step)
- **dt_min**: Minimum timestep floor (only enforced after dt grows above it)
- **dt_max**: Maximum timestep ceiling
- **accur**: Target relative change per step
- **maxtstepincrement**: Maximum factor by which dt can increase per step
- **nevac**: Vacuum density threshold — grid points with n < nevac are ignored in dt_collision calculation

## Module Overview

| Module | Description |
|--------|-------------|
| `solver.py` | BDF2 time stepper, transport equation assembly, implicit parallel losses |
| `species.py` | Species definitions, PlasmaState container |
| `mesh.py` | 1D mesh generation |
| `boundary.py` | Robin/Dirichlet boundary conditions |
| `transport.py` | Diffusion models (Bohm, gyro-geometric) |
| `reactions/` | Collision rates, ADAS data, implicit Newton reaction solver |
| `parallel.py` | Limiter and Bpol loss rate calculations |
| `io/` | JSON input, CSV output |
| `gui/` | Bokeh-based interactive plotter |
