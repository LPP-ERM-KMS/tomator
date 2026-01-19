"""
Tomator-dolfinx GUI module.

Provides:
- Interactive Bokeh-based visualization for simulation results
- Tkinter-based simulation parameter interface
"""

from .plotter import launch_plotter, stop_plotter
# from .simulation_interface import SimulationInterface, main as launch_interface

# __all__ = ['launch_plotter', 'stop_plotter', 'SimulationInterface', 'launch_interface']
__all__ = ['launch_plotter', 'stop_plotter'] #, 'SimulationInterface', 'launch_interface']
