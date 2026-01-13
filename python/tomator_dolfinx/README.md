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

The solver uses **BDF2 transport with explicit reaction sources** and implicit parallel losses:

1. **Collision Sources** — Compute ionization, recombination, charge exchange rates (ADAS data) and collision frequencies (ν)
2. **RF Power** — Compute electron heating from IC/EC waves (merged into energy sources)
3. **Collision Timestep** — Compute per-cell dt_collision from source magnitudes (for analysis only)
4. **Parallel Loss Rates** — Compute k_n, k_E loss rates for implicit LHS treatment (limiter + Bpol losses)
5. **Transport Coefficients** — Calculate D (Bohm/gyro-geometric) and V for each species
6. **Diffusion Timestep** — Compute per-cell dt_diffusion from stability criterion (for analysis only)
7. **BDF2 Coefficients** — Update time discretization weights (variable for adaptive dt)
8. **Transport + Reactions Solve** — Advance densities and energies via BDF2 FEM:
   - Ion transport (n, E) with implicit k_n, k_E on LHS + explicit reaction sources
   - Neutral transport (n, E) with explicit reaction sources
   - Electron energy transport with implicit k_E on LHS + collision/RF sources
9. **Post-processing** — Apply density floors, quasi-neutrality, temperature clamps
10. **Timestep Adaptation** — Adjust dt for next step based on max_change

## Timestep Strategy

The solver uses a **simple adaptive timestep** scheme:

1. **Initialize**: `dt = dt_init` (can be smaller than dt_min)
2. **After each step**: Compute `max_change = max(|Δn|/n, |ΔE|/E)` across non-vacuum cells
3. **Adapt dt**: `dt_next = dt * (accur / max_change)`, clamped to dt_max
4. **dt_min enforcement**: Only enforced once dt has grown above it (allows starting small)

This way:
- If `max_change > accur`: dt shrinks → smaller changes next step
- If `max_change < accur`: dt grows → larger steps while staying accurate

No step rejection — just continuous adaptation.

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
| `reactions/` | Collision rates, ADAS data |
| `parallel.py` | Limiter and Bpol loss rate calculations |
| `io/` | JSON input, CSV output |
| `gui/` | Bokeh-based interactive plotter |
