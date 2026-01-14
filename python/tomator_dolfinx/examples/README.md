# Examples

## Quick Start

```bash
cd tomator/python
conda activate t1dl-env

# Run simulation until 0.1s
python -u tomator_dolfinx/examples/run_from_json.py -t 0.1 tomator_dolfinx/examples/TCV5151X_fixPDV.json

# Plot results from a running or completed simulation
python -m tomator_dolfinx.gui.plotter /TomatorResults/TCV5151X_fixPDV/Res_*.csv
```

## Input Files

| File | Description |
|------|-------------|
| `TCV5151X_fixPDV.json` | TCV #5151X with fixed D and V, fixed power fraction |
