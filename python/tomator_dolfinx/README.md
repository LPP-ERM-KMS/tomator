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

1. **BDF2 Coefficients** — Update time discretization weights
2. **Collision Rates** — Compute ionization, recombination, charge exchange rates (ADAS data)
3. **Transport Coefficients** — Calculate D (Bohm/gyro-geometric) and V for each species
4. **Parallel Losses** — Compute limiter and Bpol losses (particle sinks at boundaries)
5. **RF Power** — Add electron heating from IC/EC waves
6. **Transport Solve** — Advance densities via FEM weak form
7. **Reaction Solve** — Apply atomic/molecular reactions (see solver modes below)
8. **Energy Solve** — Advance electron energy equation
9. **Post-processing** — Apply density floors, quasi-neutrality, temperature clamps

## Solver Modes

Two approaches are available for handling reaction sources:

### 1. Explicit Mode (`operator_splitting: false`)
All sources (collisions + parallel losses) are combined and treated **explicitly** in the transport step:

$$n^{k+1} = n^k + \Delta t \left[ \nabla \cdot (D \nabla n - V n) + S^k \right]$$

- Simpler implementation
- May require smaller timesteps for stiff reactions

### 2. Operator Splitting Mode (`operator_splitting: true`, default)
Transport and reactions are solved **separately**:

1. **Transport step**: Solve diffusion/advection with parallel losses only
2. **Reaction step**: Solve reaction ODE system **implicitly** at each mesh point using Newton iteration

$$\frac{dn}{dt} = S_{\text{reactions}}(n, T) \quad \text{(implicit solve)}$$

- More robust for stiff reaction systems
- Better conservation properties
- Allows larger timesteps

## Module Overview

| Module | Description |
|--------|-------------|
| `solver.py` | BDF2 time stepper, transport equation assembly |
| `species.py` | Species definitions, PlasmaState container |
| `mesh.py` | 1D mesh generation |
| `boundary.py` | Robin/Dirichlet boundary conditions |
| `transport.py` | Diffusion models (Bohm, gyro-geometric) |
| `reactions/` | Collision rates, ADAS data, implicit reaction solver |
| `parallel.py` | Limiter and Bpol loss calculations |
| `io/` | JSON input, CSV output |
| `gui/` | Bokeh-based interactive plotter |
