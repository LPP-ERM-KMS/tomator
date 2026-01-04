#!/usr/bin/env python3
"""
Bokeh plotting application for Tomator simulation results.

This is the Bokeh server application that displays interactive plots
of simulation results. It monitors the CSV file for changes and updates
the plots automatically.

Adapted from the original C++ Tomator1D GUI.
"""

import os
import time
import pandas as pd
import numpy as np
from bokeh.models import (
    ColumnDataSource,
    LinearAxis,
    LogScale,
    Range1d,
    TapTool,
    CustomJS,
    Span,
    Label,
    PrintfTickFormatter,
)
from bokeh.plotting import figure, curdoc
from bokeh.layouts import gridplot

# =============================================================================
# Configuration
# =============================================================================

# Get CSV file path from environment variable
DATA_FILE = os.environ.get('TOMATOR_CSV_FILE', None)

# Plot dimensions
SCALE_FACTOR = 0.7
PLOT_WIDTH = round((960 * 2 / 3) * SCALE_FACTOR)
PLOT_HEIGHT = round(445 * SCALE_FACTOR)

# Color palette for multiple species
COLORS = [
    "blue", "magenta", "green", "red", "cyan",
    "yellow", "black", "orange", "purple", "brown",
]

# Column name mappings (Python output uses different names than C++)
# Python: 't', 'R'
# C++: 'tmain', 'RadialPositions'
TIME_COL = 't'
RADIAL_COL = 'R'

# =============================================================================
# Global state
# =============================================================================

selected_radius = None
selected_time = None
last_modified_time = None
sources_rad = {}
sources_time = {}
timestamp_trigger_source = ColumnDataSource(data={'last_time': []})
df = None

# Track if user has explicitly selected a time (vs auto-following latest)
user_selected_time = False


# =============================================================================
# Utility functions
# =============================================================================

def print_blue(text):
    """Print text in blue for visibility."""
    print(f"\033[94m{text}\033[00m")


def get_closest_value(array, value):
    """Find the closest value in an array to a target value."""
    if len(array) == 0:
        return value
    closest_distance = float("inf")
    closest_value = array[0]
    for v in array:
        distance = abs(v - value)
        if distance < closest_distance:
            closest_distance = distance
            closest_value = v
    return closest_value


def should_plot_concentration(column):
    """Check if a column should be plotted as a concentration."""
    undesired_subs = ["d", "u", "C", "tnew", "minstopcrit"]
    return column.startswith('n') and not any(sub in column for sub in undesired_subs)


def should_plot_temperature(column):
    """Check if a column should be plotted as a temperature."""
    return column.startswith('T') and column != 'Te'


def should_plot_transport(column):
    """Check if a column is a transport coefficient (D or V)."""
    return column in ['D', 'V']


def should_plot_power(column):
    """Check if a column is a power profile (PRFe)."""
    return column in ['PRFe']


# =============================================================================
# Plot initialization
# =============================================================================

def initialize_plots(radial_positions):
    """Create and return the initialized plot figures."""
    
    # Get initial range from radial positions
    r_min = min(radial_positions) if radial_positions else 0.5
    r_max = max(radial_positions) if radial_positions else 1.2
    
    # Plot 1: ne & Te vs radial position - use explicit Range1d
    ne_te_plot = figure(
        width=PLOT_WIDTH,
        height=PLOT_HEIGHT,
        title="ne & Te vs Radial Position",
        x_axis_label="R [m]",
        y_axis_label="ne [m⁻³]",
        x_range=Range1d(start=r_min, end=r_max),
        y_range=Range1d(start=0, end=1e18),
    )
    
    # Plot 2: Concentration (density) vs radial position (log scale)
    concentration_plot = figure(
        width=PLOT_WIDTH,
        height=PLOT_HEIGHT,
        title="Concentration vs Radial Position",
        x_axis_label="R [m]",
        y_axis_label="Density [m⁻³]",
        y_axis_type="log",
        x_range=Range1d(start=r_min, end=r_max),
        y_range=Range1d(start=1e10, end=1e20),
    )
    
    # Plot 3: Temperature vs radial position
    temperature_plot = figure(
        width=PLOT_WIDTH,
        height=PLOT_HEIGHT,
        title="Temperature vs Radial Position",
        x_axis_label="R [m]",
        y_axis_label="Temperature [eV]",
        x_range=Range1d(start=r_min, end=r_max),
        y_range=Range1d(start=0, end=100),
    )
    
    # Plot 4: ne & Te vs time - use explicit Range1d to avoid auto-scaling issues
    ne_te_time_plot = figure(
        width=PLOT_WIDTH,
        height=PLOT_HEIGHT,
        title="ne & Te vs Time",
        x_axis_label="Time [s]",
        y_axis_label="ne [m⁻³]",
        x_range=Range1d(start=0, end=1e-3),
        y_range=Range1d(start=0, end=1e18),
    )
    ne_te_time_plot.xaxis.formatter = PrintfTickFormatter(format="%.1e")
    
    # Plot 5: Concentration vs time (log scale) - explicit Range1d
    concentration_time_plot = figure(
        width=PLOT_WIDTH,
        height=PLOT_HEIGHT,
        title="Concentration vs Time",
        x_axis_label="Time [s]",
        y_axis_label="Density [m⁻³]",
        y_axis_type="log",
        x_range=Range1d(start=0, end=1e-3),
        y_range=Range1d(start=1e10, end=1e20),
    )
    concentration_time_plot.xaxis.formatter = PrintfTickFormatter(format="%.1e")
    
    # Plot 6: Temperature vs time - explicit Range1d
    temperature_time_plot = figure(
        width=PLOT_WIDTH,
        height=PLOT_HEIGHT,
        title="Temperature vs Time",
        x_axis_label="Time [s]",
        y_axis_label="Temperature [eV]",
        x_range=Range1d(start=0, end=1e-3),
        y_range=Range1d(start=0, end=100),
    )
    temperature_time_plot.xaxis.formatter = PrintfTickFormatter(format="%.1e")
    
    # Plot 7: D & V (transport coefficients) vs radial position
    transport_plot = figure(
        width=PLOT_WIDTH,
        height=PLOT_HEIGHT,
        title="D & V vs Radial Position",
        x_axis_label="R [m]",
        y_axis_label="D [m²/s]",
        x_range=Range1d(start=r_min, end=r_max),
        y_range=Range1d(start=0, end=10),
    )
    
    # Plot 8: PRFe (coupled RF power) vs radial position
    power_plot = figure(
        width=PLOT_WIDTH,
        height=PLOT_HEIGHT,
        title="PRFe vs Radial Position",
        x_axis_label="R [m]",
        y_axis_label="PRFe [W/m³]",
        x_range=Range1d(start=r_min, end=r_max),
        y_range=Range1d(start=0, end=1e6),
    )
    
    # Plot 9: PRFe (coupled RF power) vs time
    power_time_plot = figure(
        width=PLOT_WIDTH,
        height=PLOT_HEIGHT,
        title="PRFe vs Time",
        x_axis_label="Time [s]",
        y_axis_label="PRFe [W/m³]",
        x_range=Range1d(start=0, end=1e-3),
        y_range=Range1d(start=0, end=1e6),
    )
    power_time_plot.xaxis.formatter = PrintfTickFormatter(format="%.1e")
    
    return (
        ne_te_plot,
        concentration_plot,
        temperature_plot,
        ne_te_time_plot,
        concentration_time_plot,
        temperature_time_plot,
        transport_plot,
        power_plot,
        power_time_plot,
    )


def plot_ne_te_radial(plot, radial_positions, sources_rad):
    """Set up ne and Te on the radial plot with dual y-axes."""
    sources_rad["ne"] = ColumnDataSource(
        data=dict(x=radial_positions, y=[0] * len(radial_positions))
    )
    sources_rad["Te"] = ColumnDataSource(
        data=dict(x=radial_positions, y=[0] * len(radial_positions))
    )
    
    # Add secondary y-axis for Te
    plot.extra_y_ranges = {"Te": Range1d(start=0, end=100)}
    plot.add_layout(LinearAxis(y_range_name="Te", axis_label="Te [eV]"), "right")
    
    # Plot ne on primary axis
    plot.line("x", "y", source=sources_rad["ne"], legend_label="ne", color="blue", line_width=2)
    
    # Plot Te on secondary axis
    plot.line("x", "y", source=sources_rad["Te"], legend_label="Te", color="red",
              line_width=2, y_range_name="Te")
    
    plot.legend.location = "top_left"
    plot.legend.click_policy = "hide"


def plot_ne_te_time(plot, timestamps, sources_time):
    """Set up ne and Te on the time plot with dual y-axes."""
    sources_time["ne"] = ColumnDataSource(
        data=dict(x=timestamps, y=[0] * len(timestamps))
    )
    sources_time["Te"] = ColumnDataSource(
        data=dict(x=timestamps, y=[0] * len(timestamps))
    )
    
    # Add secondary y-axis for Te
    plot.extra_y_ranges = {"Te": Range1d(start=0, end=100)}
    plot.add_layout(LinearAxis(y_range_name="Te", axis_label="Te [eV]"), "right")
    
    # Plot ne on primary axis
    plot.line("x", "y", source=sources_time["ne"], legend_label="ne", color="blue", line_width=2)
    
    # Plot Te on secondary axis
    plot.line("x", "y", source=sources_time["Te"], legend_label="Te", color="red",
              line_width=2, y_range_name="Te")
    
    plot.legend.location = "top_left"
    plot.legend.click_policy = "hide"


def plot_species_radial(plot, column, radial_positions, counter, sources_rad):
    """Add a species line to a radial plot."""
    source = ColumnDataSource(
        data=dict(x=radial_positions, y=[0] * len(radial_positions))
    )
    color_idx = counter % len(COLORS)
    plot.line("x", "y", source=source, legend_label=column, color=COLORS[color_idx], line_width=1.5)
    sources_rad[column] = source


def plot_species_time(plot, column, timestamps, counter, sources_time):
    """Add a species line to a time plot."""
    source = ColumnDataSource(
        data=dict(x=timestamps, y=[0] * len(timestamps))
    )
    color_idx = counter % len(COLORS)
    plot.line("x", "y", source=source, legend_label=column, color=COLORS[color_idx], line_width=1.5)
    sources_time[column] = source


def setup_all_plots(
    ne_te_plot, concentration_plot, temperature_plot,
    ne_te_time_plot, concentration_time_plot, temperature_time_plot,
    transport_plot, power_plot, power_time_plot,
    radial_positions, timestamps, sources_rad, sources_time
):
    """Set up all plot data sources for detected columns."""
    global df
    
    # Set up ne & Te plots
    plot_ne_te_radial(ne_te_plot, radial_positions, sources_rad)
    plot_ne_te_time(ne_te_time_plot, timestamps, sources_time)
    
    counter_conc = 0
    counter_temp = 0
    
    for column in df.columns:
        if column in [TIME_COL, RADIAL_COL, 'ne', 'Te']:
            continue
        
        if should_plot_concentration(column):
            plot_species_radial(concentration_plot, column, radial_positions, counter_conc, sources_rad)
            plot_species_time(concentration_time_plot, column, timestamps, counter_conc, sources_time)
            counter_conc += 1
        
        elif should_plot_temperature(column):
            plot_species_radial(temperature_plot, column, radial_positions, counter_temp, sources_rad)
            plot_species_time(temperature_time_plot, column, timestamps, counter_temp, sources_time)
            counter_temp += 1
    
    # Set up D & V transport plot with dual y-axes
    if 'D' in df.columns:
        sources_rad['D'] = ColumnDataSource(
            data=dict(x=radial_positions, y=[0] * len(radial_positions))
        )
        transport_plot.line("x", "y", source=sources_rad['D'], legend_label="D (ion)", color="blue", line_width=2)
    
    # Add neutral diffusion coefficients (different colors) - scaled by 1/1000 to fit on same axis
    NEUTRAL_D_SCALE = 0.001  # Neutral D is ~1000x larger than ion D
    neutral_colors = {'D_H': 'green', 'D_H2': 'purple', 'D_HeI': 'orange'}
    for col_name, color in neutral_colors.items():
        if col_name in df.columns:
            sources_rad[col_name] = ColumnDataSource(
                data=dict(x=radial_positions, y=[0] * len(radial_positions))
            )
            transport_plot.line("x", "y", source=sources_rad[col_name], 
                               legend_label=f"{col_name}/1000", color=color, line_width=1.5, line_dash="dashed")
    
    if 'V' in df.columns:
        # V can be negative (inward pinch), so use secondary axis
        transport_plot.extra_y_ranges = {"V": Range1d(start=-10, end=10)}
        transport_plot.add_layout(LinearAxis(y_range_name="V", axis_label="V [m/s]"), "right")
        sources_rad['V'] = ColumnDataSource(
            data=dict(x=radial_positions, y=[0] * len(radial_positions))
        )
        transport_plot.line("x", "y", source=sources_rad['V'], legend_label="V", color="red",
                            line_width=2, y_range_name="V")
    
    transport_plot.legend.location = "top_left"
    transport_plot.legend.click_policy = "hide"
    
    # Set up PRFe power plot
    if 'PRFe' in df.columns:
        sources_rad['PRFe'] = ColumnDataSource(
            data=dict(x=radial_positions, y=[0] * len(radial_positions))
        )
        power_plot.line("x", "y", source=sources_rad['PRFe'], legend_label="PRFe", color="orange", line_width=2)
        power_plot.legend.location = "top_left"
        power_plot.legend.click_policy = "hide"
        
        # PRFe vs time
        sources_time['PRFe'] = ColumnDataSource(
            data=dict(x=timestamps, y=[0] * len(timestamps))
        )
        power_time_plot.line("x", "y", source=sources_time['PRFe'], legend_label="PRFe", color="orange", line_width=2)
        power_time_plot.legend.location = "top_left"
        power_time_plot.legend.click_policy = "hide"
    
    # Configure legends
    for plot in [concentration_plot, temperature_plot, concentration_time_plot, temperature_time_plot]:
        plot.legend.location = "top_left"
        plot.legend.click_policy = "hide"


# =============================================================================
# Interactive elements
# =============================================================================

def setup_interactive_elements(ne_te_plot, ne_te_time_plot):
    """Set up interactive selection elements (vertical lines, click handlers)."""
    global sources_rad, sources_time
    
    # Movable point source for radial plot
    movable_point_radial = ColumnDataSource(data=dict(x=[0], y=[0]))
    ne_te_plot.circle("x", "y", source=movable_point_radial, size=8, color="black")
    
    # Vertical line for radial selection
    vertical_line_radial = Span(
        location=0, dimension="height",
        line_color="black", line_dash="dashed", line_width=1
    )
    ne_te_plot.add_layout(vertical_line_radial)
    
    # Label for radial selection
    label_radial = Label(x=0, y=0, text="", text_font_size="10pt", text_color="black")
    ne_te_plot.add_layout(label_radial)
    
    # Movable point source for time plot
    movable_point_time = ColumnDataSource(data=dict(x=[0], y=[0]))
    ne_te_time_plot.circle("x", "y", source=movable_point_time, size=8, color="black")
    
    # Vertical line for time selection
    vertical_line_time = Span(
        location=0, dimension="height",
        line_color="black", line_dash="dashed", line_width=1
    )
    ne_te_time_plot.add_layout(vertical_line_time)
    
    # Label for time selection
    label_time = Label(x=0, y=0, text="", text_font_size="10pt", text_color="black")
    ne_te_time_plot.add_layout(label_time)
    
    # Source for communicating selected values to Python callbacks
    selected_radius_source = ColumnDataSource(data=dict(radius=[0]))
    selected_time_source = ColumnDataSource(data=dict(timestamp=[0]))
    
    return (
        movable_point_radial, vertical_line_radial, label_radial,
        movable_point_time, vertical_line_time, label_time,
        selected_radius_source, selected_time_source,
    )


# =============================================================================
# Update functions
# =============================================================================

def scale_radial_plots(ne_te_plot, concentration_plot, temperature_plot, transport_plot, power_plot, time_val=None):
    """Scale radial plot axes based on current data."""
    global df
    
    if df is None or len(df) == 0:
        return
    
    if time_val is None:
        time_val = df[TIME_COL].max()
    
    df_slice = df[df[TIME_COL] == time_val]
    if len(df_slice) == 0:
        return
    
    # Scale x range (same for all radial plots)
    r_min = df_slice[RADIAL_COL].min()
    r_max = df_slice[RADIAL_COL].max()
    dr = r_max - r_min if r_max > r_min else 0.01
    x_start = r_min - 0.05 * dr
    x_end = r_max + 0.05 * dr
    
    # Scale ne/Te plot - auto-scale based on data range
    min_ne = df_slice['ne'].min() if 'ne' in df_slice else 0
    max_ne = df_slice['ne'].max() if 'ne' in df_slice else 1e18
    min_Te = df_slice['Te'].min() if 'Te' in df_slice else 0
    max_Te = df_slice['Te'].max() if 'Te' in df_slice else 10
    
    # Add 10% padding, but keep minimum at 0 for physical quantities
    ne_range = max_ne - min_ne if max_ne > min_ne else max_ne * 0.1
    Te_range = max_Te - min_Te if max_Te > min_Te else max_Te * 0.1
    
    ne_te_plot.y_range.start = max(0, min_ne - ne_range * 0.05)
    ne_te_plot.y_range.end = max_ne + ne_range * 0.1
    ne_te_plot.extra_y_ranges["Te"].start = max(0, min_Te - Te_range * 0.05)
    ne_te_plot.extra_y_ranges["Te"].end = max_Te + Te_range * 0.1
    ne_te_plot.x_range.start = x_start
    ne_te_plot.x_range.end = x_end
    
    # Scale concentration plot (log scale) - find min/max of all concentration columns
    conc_cols = [c for c in df_slice.columns if should_plot_concentration(c) and c != 'ne']
    if conc_cols:
        all_conc_values = []
        for col in conc_cols:
            vals = df_slice[col].values
            # Filter out zeros and negatives for log scale
            positive_vals = vals[vals > 0]
            if len(positive_vals) > 0:
                all_conc_values.extend(positive_vals)
        
        if all_conc_values:
            min_conc = min(all_conc_values)
            max_conc = max(all_conc_values)
            # Log scale: set range with some padding in log space
            concentration_plot.y_range.start = min_conc * 0.5
            concentration_plot.y_range.end = max_conc * 2.0
    
    concentration_plot.x_range.start = x_start
    concentration_plot.x_range.end = x_end
    
    # Scale temperature plot - find min/max of all temperature columns
    temp_cols = [c for c in df_slice.columns if should_plot_temperature(c)]
    if temp_cols:
        all_temp_values = []
        for col in temp_cols:
            all_temp_values.extend(df_slice[col].values)
        
        if all_temp_values:
            min_temp = min(all_temp_values)
            max_temp = max(all_temp_values)
            temp_range = max_temp - min_temp if max_temp > min_temp else max(max_temp * 0.1, 0.1)
            temperature_plot.y_range.start = min_temp - temp_range * 0.05
            temperature_plot.y_range.end = max_temp + temp_range * 0.1
    
    temperature_plot.x_range.start = x_start
    temperature_plot.x_range.end = x_end
    
    # Scale transport plot (D & V) - include scaled neutral diffusion
    NEUTRAL_D_SCALE = 0.001  # Same scale factor as in data update
    all_D_vals = []
    if 'D' in df_slice.columns:
        all_D_vals.extend(df_slice['D'].values.tolist())
    # Include scaled neutral diffusion values
    for col_name in ['D_H', 'D_H2', 'D_HeI']:
        if col_name in df_slice.columns:
            all_D_vals.extend((df_slice[col_name].values * NEUTRAL_D_SCALE).tolist())
    
    if all_D_vals:
        min_D = min(all_D_vals)
        max_D = max(all_D_vals)
        D_range = max_D - min_D if max_D > min_D else max(max_D * 0.1, 0.1)
        transport_plot.y_range.start = max(0, min_D - D_range * 0.05)
        transport_plot.y_range.end = max_D + D_range * 0.1
    
    if 'V' in df_slice.columns and hasattr(transport_plot, 'extra_y_ranges') and 'V' in transport_plot.extra_y_ranges:
        V_vals = df_slice['V'].values
        min_V = V_vals.min()
        max_V = V_vals.max()
        V_range = max_V - min_V if max_V > min_V else max(abs(max_V) * 0.1, 0.1)
        transport_plot.extra_y_ranges["V"].start = min_V - V_range * 0.1
        transport_plot.extra_y_ranges["V"].end = max_V + V_range * 0.1
    
    transport_plot.x_range.start = x_start
    transport_plot.x_range.end = x_end
    
    # Scale power plot (PRFe)
    if 'PRFe' in df_slice.columns:
        PRFe_vals = df_slice['PRFe'].values
        min_P = PRFe_vals.min()
        max_P = PRFe_vals.max()
        P_range = max_P - min_P if max_P > min_P else max(max_P * 0.1, 1.0)
        power_plot.y_range.start = max(0, min_P - P_range * 0.05)
        power_plot.y_range.end = max_P + P_range * 0.1
    
    power_plot.x_range.start = x_start
    power_plot.x_range.end = x_end


def scale_time_plots(ne_te_time_plot, concentration_time_plot, temperature_time_plot, power_time_plot=None, radius_val=None):
    """Scale time plot axes based on current data."""
    global df
    
    if df is None or len(df) == 0:
        return
    
    radial_positions = df[RADIAL_COL].unique()
    if radius_val is None:
        radius_val = get_closest_value(radial_positions, radial_positions.mean())
    
    df_slice = df[df[RADIAL_COL] == radius_val]
    if len(df_slice) == 0:
        return
    
    # Scale x range (same for all time plots) - use full time range from dataframe
    t_min = df[TIME_COL].min()
    t_max = df[TIME_COL].max()
    
    # Debug output
    print(f"[scale_time_plots] t_min={t_min}, t_max={t_max}, radius={radius_val}")
    
    # Ensure valid range (t_max must be > 0)
    if t_max <= 0:
        t_max = 1e-6
    t_range = t_max - t_min if t_max > t_min else t_max * 0.1
    
    x_start = 0.0  # Time always starts at 0
    x_end = t_max * 1.05  # 5% padding on the right
    
    # Ensure x_end > x_start
    if x_end <= x_start:
        x_end = x_start + 1e-6
    
    # Scale ne/Te plot - auto-scale based on data range
    min_ne = df_slice['ne'].min() if 'ne' in df_slice else 0
    max_ne = df_slice['ne'].max() if 'ne' in df_slice else 1e18
    min_Te = df_slice['Te'].min() if 'Te' in df_slice else 0
    max_Te = df_slice['Te'].max() if 'Te' in df_slice else 10
    
    # Add 10% padding, but keep minimum at 0 for physical quantities
    ne_range = max_ne - min_ne if max_ne > min_ne else max(max_ne * 0.1, 1e10)
    Te_range = max_Te - min_Te if max_Te > min_Te else max(max_Te * 0.1, 1.0)
    
    # Ensure valid y-ranges (end > start)
    y_ne_start = max(0, min_ne - ne_range * 0.05)
    y_ne_end = max(y_ne_start + 1e10, max_ne + ne_range * 0.1)
    y_Te_start = max(0, min_Te - Te_range * 0.05)
    y_Te_end = max(y_Te_start + 0.1, max_Te + Te_range * 0.1)
    
    ne_te_time_plot.y_range.start = y_ne_start
    ne_te_time_plot.y_range.end = y_ne_end
    ne_te_time_plot.extra_y_ranges["Te"].start = y_Te_start
    ne_te_time_plot.extra_y_ranges["Te"].end = y_Te_end
    
    # Debug: show x_range before and after setting
    print(f"[scale_time_plots] Setting x_range: x_start={x_start}, x_end={x_end}")
    print(f"[scale_time_plots] Before: ne_te_time_plot.x_range.start={ne_te_time_plot.x_range.start}, end={ne_te_time_plot.x_range.end}")
    
    ne_te_time_plot.x_range.start = x_start
    ne_te_time_plot.x_range.end = x_end
    
    print(f"[scale_time_plots] After: ne_te_time_plot.x_range.start={ne_te_time_plot.x_range.start}, end={ne_te_time_plot.x_range.end}")
    
    # Scale concentration plot (log scale) - find min/max of all concentration columns
    conc_cols = [c for c in df_slice.columns if should_plot_concentration(c) and c != 'ne']
    if conc_cols:
        all_conc_values = []
        for col in conc_cols:
            vals = df_slice[col].values
            # Filter out zeros and negatives for log scale
            positive_vals = vals[vals > 0]
            if len(positive_vals) > 0:
                all_conc_values.extend(positive_vals)
        
        if all_conc_values:
            min_conc = min(all_conc_values)
            max_conc = max(all_conc_values)
            # Log scale: set range with some padding in log space
            concentration_time_plot.y_range.start = min_conc * 0.5
            concentration_time_plot.y_range.end = max(min_conc * 0.6, max_conc * 2.0)
    
    concentration_time_plot.x_range.start = x_start
    concentration_time_plot.x_range.end = x_end
    
    # Scale temperature plot - find min/max of all temperature columns
    temp_cols = [c for c in df_slice.columns if should_plot_temperature(c)]
    if temp_cols:
        all_temp_values = []
        for col in temp_cols:
            all_temp_values.extend(df_slice[col].values)
        
        if all_temp_values:
            min_temp = min(all_temp_values)
            max_temp = max(all_temp_values)
            temp_range = max_temp - min_temp if max_temp > min_temp else max(max_temp * 0.1, 0.1)
            y_temp_start = min_temp - temp_range * 0.05
            y_temp_end = max(y_temp_start + 0.1, max_temp + temp_range * 0.1)
            temperature_time_plot.y_range.start = y_temp_start
            temperature_time_plot.y_range.end = y_temp_end
    
    temperature_time_plot.x_range.start = x_start
    temperature_time_plot.x_range.end = x_end
    
    # Scale power time plot
    if power_time_plot is not None and 'PRFe' in df_slice.columns:
        PRFe_vals = df_slice['PRFe'].values
        min_P = PRFe_vals.min()
        max_P = PRFe_vals.max()
        P_range = max_P - min_P if max_P > min_P else max(max_P * 0.1, 1.0)
        power_time_plot.y_range.start = max(0, min_P - P_range * 0.05)
        power_time_plot.y_range.end = max(1.0, max_P + P_range * 0.1)
        power_time_plot.x_range.start = x_start
        power_time_plot.x_range.end = x_end


def update_data(
    radial_positions,
    ne_te_plot, concentration_plot, temperature_plot,
    ne_te_time_plot, concentration_time_plot, temperature_time_plot,
    transport_plot, power_plot, power_time_plot,
):
    """Periodic callback to check for file changes and update plots."""
    global DATA_FILE, selected_radius, selected_time, last_modified_time, df
    global sources_rad, sources_time, timestamp_trigger_source, user_selected_time
    
    if DATA_FILE is None:
        return
    
    # Check if file exists and has been modified
    if not os.path.exists(DATA_FILE):
        return
    
    try:
        current_modified_time = os.path.getmtime(DATA_FILE)
    except OSError:
        return
    
    if last_modified_time == current_modified_time:
        return
    
    last_modified_time = current_modified_time
    
    # Read updated data
    try:
        df = pd.read_csv(DATA_FILE)
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        return
    
    if len(df) == 0:
        return
    
    # Update radial positions list from current data
    current_radial_positions = df[RADIAL_COL].unique().tolist()
    current_timestamps = df[TIME_COL].unique().tolist()
    
    # Initialize selected values if not set
    if selected_radius is None:
        selected_radius = get_closest_value(current_radial_positions, np.mean(current_radial_positions))
    else:
        # Ensure selected_radius is still valid in current data
        selected_radius = get_closest_value(current_radial_positions, selected_radius)
    
    # Get latest time
    last_time = df[TIME_COL].max()
    
    # Determine which time to use for radial plots:
    # - If user hasn't selected a time, follow the latest (live mode)
    # - If user has selected a time, keep showing that time
    if not user_selected_time or selected_time is None:
        selected_time = last_time
        display_time = last_time
    else:
        # User selected a specific time - keep it if still valid
        display_time = get_closest_value(current_timestamps, selected_time)
        selected_time = display_time
    
    # Get data slice for radial plots at display_time
    df_radial = df[df[TIME_COL] == display_time]
    
    # Update plot titles
    time_ms = display_time * 1000
    ne_te_plot.title.text = f"ne & Te vs Radial Position at t = {time_ms:.2f} ms"
    concentration_plot.title.text = f"Concentration vs Radial Position at t = {time_ms:.2f} ms"
    temperature_plot.title.text = f"Temperature vs Radial Position at t = {time_ms:.2f} ms"
    
    ne_te_time_plot.title.text = f"ne & Te vs Time at R = {selected_radius:.4f} m"
    concentration_time_plot.title.text = f"Concentration vs Time at R = {selected_radius:.4f} m"
    temperature_time_plot.title.text = f"Temperature vs Time at R = {selected_radius:.4f} m"
    
    # Update radial sources with data at display_time
    NEUTRAL_D_SCALE = 0.001  # Same scale factor as in setup
    for column in sources_rad:
        if column in df_radial.columns:
            y_data = df_radial[column].tolist()
            # Scale neutral diffusion coefficients to fit on same axis as ion D
            if column in ['D_H', 'D_H2', 'D_HeI']:
                y_data = [v * NEUTRAL_D_SCALE for v in y_data]
            sources_rad[column].data = dict(
                x=current_radial_positions,
                y=y_data
            )
    
    # Update time sources with current timestamps for selected radius
    df_radius = df[df[RADIAL_COL] == selected_radius]
    time_values = df_radius[TIME_COL].tolist()
    for column in sources_time:
        if column in df_radius.columns:
            sources_time[column].data = dict(
                x=time_values,
                y=df_radius[column].tolist()
            )
    
    # Trigger timestamp update
    timestamp_trigger_source.data = {'last_time': [last_time]}
    
    # Scale plots based on displayed data
    scale_radial_plots(ne_te_plot, concentration_plot, temperature_plot, transport_plot, power_plot, display_time)
    scale_time_plots(ne_te_time_plot, concentration_time_plot, temperature_time_plot, power_time_plot, selected_radius)


# =============================================================================
# Main document setup
# =============================================================================

def modify_doc(doc):
    """Main function to set up the Bokeh document."""
    global df, DATA_FILE, selected_radius, selected_time
    global sources_rad, sources_time
    
    # Wait for data file to exist and have content
    wait_attempts = 0
    while DATA_FILE is None or not os.path.exists(DATA_FILE) or os.path.getsize(DATA_FILE) == 0:
        if wait_attempts == 0:
            print_blue(f"Waiting for data file: {DATA_FILE}")
        time.sleep(1)
        wait_attempts += 1
        if wait_attempts > 60:
            print("Timeout waiting for data file.")
            return
    
    print_blue(f"Loading data from: {DATA_FILE}")
    
    # Load initial data
    df = pd.read_csv(DATA_FILE)
    radial_positions = df[RADIAL_COL].unique().tolist()
    timestamps = df[TIME_COL].unique().tolist()
    
    # Initialize selection
    selected_radius = get_closest_value(radial_positions, np.mean(radial_positions))
    selected_time = df[TIME_COL].max()
    
    # Create plots
    (
        ne_te_plot, concentration_plot, temperature_plot,
        ne_te_time_plot, concentration_time_plot, temperature_time_plot,
        transport_plot, power_plot, power_time_plot,
    ) = initialize_plots(radial_positions)
    
    # Set up data sources and plot lines
    setup_all_plots(
        ne_te_plot, concentration_plot, temperature_plot,
        ne_te_time_plot, concentration_time_plot, temperature_time_plot,
        transport_plot, power_plot, power_time_plot,
        radial_positions, timestamps, sources_rad, sources_time
    )
    
    # Set up interactive elements
    (
        movable_point_radial, vertical_line_radial, label_radial,
        movable_point_time, vertical_line_time, label_time,
        selected_radius_source, selected_time_source,
    ) = setup_interactive_elements(ne_te_plot, ne_te_time_plot)
    
    # JavaScript callback for radius selection (click on radial plots)
    callback_radial = CustomJS(
        args=dict(
            source=movable_point_radial,
            line=vertical_line_radial,
            label=label_radial,
            ne_source=sources_rad.get("ne"),
            selected_radius_source=selected_radius_source,
        ),
        code="""
            const data = source.data;
            const ne_data = ne_source.data;
            const x = cb_obj.x;

            // Find closest x value
            let closest_x = ne_data['x'][0];
            let closest_y = ne_data['y'][0];
            let min_dist = Math.abs(ne_data['x'][0] - x);
            
            for (let i = 1; i < ne_data['x'].length; i++) {
                const dist = Math.abs(ne_data['x'][i] - x);
                if (dist < min_dist) {
                    min_dist = dist;
                    closest_x = ne_data['x'][i];
                    closest_y = ne_data['y'][i];
                }
            }

            data['x'][0] = closest_x;
            data['y'][0] = closest_y;
            line.location = closest_x;
            label.x = closest_x;
            label.y = closest_y * 1.05;
            label.text = "R = " + closest_x.toFixed(4) + " m";
            
            selected_radius_source.data = {radius: [closest_x]};
            selected_radius_source.change.emit();
            source.change.emit();
        """,
    )
    
    # JavaScript callback for time selection (click on time plots)
    callback_time = CustomJS(
        args=dict(
            source=movable_point_time,
            line=vertical_line_time,
            label=label_time,
            ne_source=sources_time.get("ne"),
            selected_time_source=selected_time_source,
            timestamp_trigger_source=timestamp_trigger_source,
        ),
        code="""
            function updateDisplay(x) {
                const data = source.data;
                const ne_data = ne_source.data;

                // Find closest x value
                let closest_x = ne_data['x'][0];
                let closest_y = ne_data['y'][0];
                let min_dist = Math.abs(ne_data['x'][0] - x);
                
                for (let i = 1; i < ne_data['x'].length; i++) {
                    const dist = Math.abs(ne_data['x'][i] - x);
                    if (dist < min_dist) {
                        min_dist = dist;
                        closest_x = ne_data['x'][i];
                        closest_y = ne_data['y'][i];
                    }
                }

                data['x'][0] = closest_x;
                data['y'][0] = closest_y;
                line.location = closest_x;
                label.x = closest_x;
                label.y = closest_y * 1.05;
                label.text = "t = " + (closest_x * 1000).toFixed(2) + " ms";
                
                selected_time_source.data = {timestamp: [closest_x]};
                selected_time_source.change.emit();
                source.change.emit();
            }
            
            if (cb_obj && cb_obj.x) {
                updateDisplay(cb_obj.x);
            } else {
                updateDisplay(timestamp_trigger_source.data['last_time'][0]);
            }
        """,
    )
    
    # Python callback for radius selection
    def on_radius_selected(attr, old, new):
        global selected_radius, df, sources_time
        
        # Get current radial positions from dataframe
        current_radial_positions = df[RADIAL_COL].unique().tolist()
        
        selected_radius = selected_radius_source.data["radius"][0]
        selected_radius = get_closest_value(current_radial_positions, selected_radius)
        
        # Update time sources for new radius
        df_radius = df[df[RADIAL_COL] == selected_radius]
        time_values = df_radius[TIME_COL].tolist()
        for column in sources_time:
            if column in df_radius.columns:
                sources_time[column].data = dict(
                    x=time_values,
                    y=df_radius[column].tolist()
                )
        
        # Update titles
        ne_te_time_plot.title.text = f"ne & Te vs Time at R = {selected_radius:.4f} m"
        concentration_time_plot.title.text = f"Concentration vs Time at R = {selected_radius:.4f} m"
        temperature_time_plot.title.text = f"Temperature vs Time at R = {selected_radius:.4f} m"
        power_time_plot.title.text = f"PRFe vs Time at R = {selected_radius:.4f} m"
        
        scale_time_plots(ne_te_time_plot, concentration_time_plot, temperature_time_plot, power_time_plot, selected_radius)
    
    # Python callback for time selection
    def on_time_selected(attr, old, new):
        global selected_time, df, sources_rad, user_selected_time
        
        selected_time = selected_time_source.data["timestamp"][0]
        selected_time = get_closest_value(df[TIME_COL].unique(), selected_time)
        
        # Mark that user has explicitly selected a time (disable auto-follow)
        user_selected_time = True
        
        # Update radial sources for new time
        df_time = df[df[TIME_COL] == selected_time]
        radial_pos = df_time[RADIAL_COL].tolist()
        NEUTRAL_D_SCALE = 0.001  # Same scale factor as elsewhere
        for column in sources_rad:
            if column in df_time.columns:
                y_data = df_time[column].tolist()
                # Scale neutral diffusion coefficients to fit on same axis as ion D
                if column in ['D_H', 'D_H2', 'D_HeI']:
                    y_data = [v * NEUTRAL_D_SCALE for v in y_data]
                sources_rad[column].data = dict(
                    x=radial_pos,
                    y=y_data
                )
        
        # Update titles
        time_ms = selected_time * 1000
        ne_te_plot.title.text = f"ne & Te vs Radial Position at t = {time_ms:.2f} ms"
        concentration_plot.title.text = f"Concentration vs Radial Position at t = {time_ms:.2f} ms"
        temperature_plot.title.text = f"Temperature vs Radial Position at t = {time_ms:.2f} ms"
        
        scale_radial_plots(ne_te_plot, concentration_plot, temperature_plot, transport_plot, power_plot, selected_time)
    
    # Attach callbacks
    ne_te_plot.js_on_event("tap", callback_radial)
    selected_radius_source.on_change("data", on_radius_selected)
    
    ne_te_time_plot.js_on_event("tap", callback_time)
    timestamp_trigger_source.js_on_change("data", callback_time)
    selected_time_source.on_change("data", on_time_selected)
    
    # Initial data update
    update_data(
        radial_positions,
        ne_te_plot, concentration_plot, temperature_plot,
        ne_te_time_plot, concentration_time_plot, temperature_time_plot,
        transport_plot, power_plot, power_time_plot,
    )
    
    # Periodic update callback (every 2 seconds)
    doc.add_periodic_callback(
        lambda: update_data(
            radial_positions,
            ne_te_plot, concentration_plot, temperature_plot,
            ne_te_time_plot, concentration_time_plot, temperature_time_plot,
            transport_plot, power_plot, power_time_plot,
        ),
        2000,
    )
    
    # Arrange plots in grid (3x3 with transport and power on bottom row)
    grid = gridplot([
        [ne_te_plot, concentration_plot, temperature_plot],
        [ne_te_time_plot, concentration_time_plot, temperature_time_plot],
        [transport_plot, power_plot, power_time_plot],
    ])
    
    doc.add_root(grid)
    doc.title = "Tomator Results"


# =============================================================================
# Entry point
# =============================================================================

if DATA_FILE:
    print_blue(f"Tomator Plotter - Monitoring: {DATA_FILE}")
    modify_doc(curdoc())
else:
    print("Error: TOMATOR_CSV_FILE environment variable not set.")
