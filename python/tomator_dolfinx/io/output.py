"""
Output file writing for Tomator simulations.

Writes simulation results to CSV files compatible with the
C++ Tomator1D output format for comparison and post-processing.
"""

import os
import csv
from datetime import datetime
from typing import Optional
import numpy as np

from ..species import PlasmaState


def write_csv_output(
    state: PlasmaState,
    t: float,
    radial_positions: np.ndarray,
    output_dir: str,
    filename: Optional[str] = None,
    append: bool = True
) -> str:
    """
    Write current state to CSV file.
    
    Creates output compatible with C++ Tomator1D format:
    tmain, RadialPositions, ne, Ee, Te, nH, EH, TH, nHi, EHi, THi, ...
    
    Parameters
    ----------
    state : PlasmaState
        Current plasma state.
    t : float
        Current time [s].
    radial_positions : np.ndarray
        Radial positions [m].
    output_dir : str
        Output directory path.
    filename : str, optional
        Specific filename. If None, auto-generated with timestamp.
    append : bool
        If True, append to existing file. If False, overwrite.
        
    Returns
    -------
    filepath : str
        Path to written file.
    """
    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate filename if not provided
    if filename is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"Res_{timestamp}.csv"
    
    filepath = os.path.join(output_dir, filename)
    
    # Build header and data
    header = ['t', 'R']
    data_columns = [
        np.full(len(radial_positions), t),
        radial_positions
    ]
    
    # Add species data
    for name in ['e', 'H', 'Hi', 'HeI', 'HeII', 'HeIII', 'H2', 'H2i', 'H3i']:
        if name in state.species:
            species = state.species[name]
            n = species.n.x.array
            E = species.E.x.array
            T = species.T
            
            header.extend([f'n{name}', f'E{name}', f'T{name}'])
            data_columns.extend([n, E, T])
    
    # Stack data
    data = np.column_stack(data_columns)
    
    # Determine write mode
    file_exists = os.path.exists(filepath)
    mode = 'a' if append and file_exists else 'w'
    write_header = (mode == 'w')
    
    # Write CSV
    with open(filepath, mode, newline='') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(header)
        writer.writerows(data)
    
    return filepath


def write_time_series(
    state: PlasmaState,
    t: float,
    output_file: str,
    radial_index: int = None,
    append: bool = True
) -> None:
    """
    Write time series data at a specific radial location.
    
    Parameters
    ----------
    state : PlasmaState
        Current plasma state.
    t : float
        Current time [s].
    output_file : str
        Output file path.
    radial_index : int, optional
        Radial mesh index to extract. If None, uses center.
    append : bool
        Append to existing file (True) or overwrite (False).
    """
    # Get radial index (default: center)
    nmesh = len(state.electrons.n.x.array)
    if radial_index is None:
        radial_index = nmesh // 2
    
    # Build header
    header = ['t']
    values = [t]
    
    for name in ['e', 'H', 'Hi', 'HeI', 'HeII', 'HeIII']:
        if name in state.species:
            species = state.species[name]
            header.extend([f'n{name}', f'T{name}'])
            values.extend([
                species.n.x.array[radial_index],
                species.T[radial_index]
            ])
    
    # Write
    mode = 'a' if append and os.path.exists(output_file) else 'w'
    write_header = (mode == 'w')
    
    with open(output_file, mode, newline='') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(header)
        writer.writerow(values)


def write_snapshot(
    state: PlasmaState,
    t: float,
    radial_positions: np.ndarray,
    output_dir: str
) -> str:
    """
    Write a complete snapshot of the simulation state.
    
    Creates a detailed output file with all species data.
    
    Parameters
    ----------
    state : PlasmaState
        Current plasma state.
    t : float
        Current time [s].
    radial_positions : np.ndarray
        Radial positions [m].
    output_dir : str
        Output directory path.
        
    Returns
    -------
    filepath : str
        Path to written file.
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"snapshot_t{t:.6e}_{timestamp}.csv"
    
    return write_csv_output(state, t, radial_positions, output_dir, filename)


class OutputManager:
    """
    Manager for simulation output with configurable intervals.
    
    Attributes
    ----------
    output_dir : str
        Output directory.
    save_interval : float
        Time interval between saves [s].
    last_save_time : float
        Time of last save.
    """
    
    def __init__(
        self,
        output_dir: str,
        save_interval: float = 1e-4,
        radial_positions: np.ndarray = None
    ):
        """
        Initialize output manager.
        
        Parameters
        ----------
        output_dir : str
            Output directory path.
        save_interval : float
            Time interval between saves [s].
        radial_positions : np.ndarray, optional
            Radial positions array.
        """
        self.output_dir = output_dir
        self.save_interval = save_interval
        self.radial_positions = radial_positions
        self.last_save_time = -save_interval  # Ensure first step is saved
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Time series file
        self.time_series_file = os.path.join(output_dir, "time_series.csv")
        self.time_series_initialized = False
    
    def should_save(self, t: float) -> bool:
        """Check if we should save at current time."""
        return (t - self.last_save_time) >= self.save_interval
    
    def save(self, state: PlasmaState, t: float) -> None:
        """
        Save current state if interval has passed.
        
        Parameters
        ----------
        state : PlasmaState
            Current plasma state.
        t : float
            Current time [s].
        """
        if not self.should_save(t):
            return
        
        # Save full profile
        if self.radial_positions is not None:
            write_csv_output(state, t, self.radial_positions, self.output_dir)
        
        # Save time series
        write_time_series(
            state, t, self.time_series_file,
            append=self.time_series_initialized
        )
        self.time_series_initialized = True
        
        self.last_save_time = t
    
    def save_final(self, state: PlasmaState, t: float) -> None:
        """
        Save final state regardless of interval.
        
        Parameters
        ----------
        state : PlasmaState
            Final plasma state.
        t : float
            Final time [s].
        """
        if self.radial_positions is not None:
            filepath = write_snapshot(state, t, self.radial_positions, self.output_dir)
            print(f"Final snapshot saved: {filepath}")
