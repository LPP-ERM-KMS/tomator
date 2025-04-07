import os
import time
import pandas as pd
from bokeh.models import (
    ColumnDataSource,
    LinearAxis,
    LogScale,
    Range1d,
    TapTool,
    CustomJS,
    Span,
    Label,
)
from bokeh.plotting import figure, curdoc
from bokeh.layouts import gridplot

# Constants
scale_factor = 1
DATA_FILE = None
PLOT_WIDTH = round((960 * 2 / 3) * scale_factor)
PLOT_HEIGHT = round(445 * scale_factor)
COLORS = [
    "blue",
    "magenta",
    "green",
    "red",
    "cyan",
    "yellow",
    "black",
    "orange",
    "purple",
    "brown",
]
selected_radius = None
selected_time = None
last_modified_time = None
x_value = 0
electron_max_vs_timestamp_line = None
sources_rad = {}
sources_time = {}
timestamp_trigger_source = ColumnDataSource(data={'last_tmain': []})
df = None

def printBlue(text):
    print("\033[94m{}\033[00m" .format(text))


def get_closest_value(array, value):
    closest_distance = float("inf")
    closest_value = None
    for v in array:
        distance = abs(v - value)
        if distance < closest_distance:
            closest_distance = distance
            closest_value = v

    return closest_value


def calculate_temperature(df):
    for column in df.columns:
        if is_energy_column(column):
            element_suffix = column[1:]
            density_col = f"n{element_suffix}"
            if density_col in df.columns:
                temp_col = f"T{element_suffix}"
                df[temp_col] = (2 / 3) * (df[column] / df[density_col])

def make_log_axis(plot):
    plot.yaxis[0].formatter.use_scientific = True
    plot.yaxis[0].formatter.power_limit_high = 0
    plot.yaxis[0].formatter.power_limit_low = 0
    plot.yaxis[0] = LogAxis()

def plot_concentration(plot, column, radial_positions, counter, sources_rad):
    source = ColumnDataSource(
        data=dict(x=radial_positions, y=[0] * len(radial_positions))
    )
    color_cycle_index = counter % len(COLORS)
    plot.line(
        "x", "y", source=source, legend_label=column, color=COLORS[color_cycle_index]
    )
    sources_rad[column] = source


def plot_energy(plot, column, radial_positions, counter, sources_rad):
    source = ColumnDataSource(
        data=dict(x=radial_positions, y=[0] * len(radial_positions))
    )
    color_cycle_index = counter % len(COLORS)
    plot.line(
        "x", "y", source=source, legend_label=column, color=COLORS[color_cycle_index]
    )
    sources_rad[column] = source


def plot_concentration_time(plot, column, counter, sources_time):
    global df
    source = ColumnDataSource(
        data=dict(
            x=df["tmain"].unique().tolist(), y=[0] * len(df["tmain"].unique().tolist())
        )
    )
    color_cycle_index = counter % len(COLORS)
    plot.line(
        "x", "y", source=source, legend_label=column, color=COLORS[color_cycle_index]
    )
    sources_time[column] = source


def plot_energy_time(plot, column, counter, sources_time):
    global df
    source = ColumnDataSource(
        data=dict(
            x=df["tmain"].unique().tolist(), y=[0] * len(df["tmain"].unique().tolist())
        )
    )
    color_cycle_index = counter % len(COLORS)
    # Plot log scale for energy
    plot.line(
        "x", "y", source=source, legend_label=column, color=COLORS[color_cycle_index]
    )

    sources_time[column] = source


def plot_ne_time(electron_max_vs_timestamp, time_stamps, sources_time):
    sources_time["ne"] = ColumnDataSource(
        data=dict(x=time_stamps, y=[0] * len(time_stamps))
    )
    sources_time["Te"] = ColumnDataSource(
        data=dict(x=time_stamps, y=[0] * len(time_stamps))
    )

    electron_max_vs_timestamp.extra_y_ranges = {"Te": Range1d(start=0, end=1)}

    electron_max_vs_timestamp.line(
        "x", "y", source=sources_time["ne"], legend_label="ne", color="blue"
    )

    electron_max_vs_timestamp.add_layout(LinearAxis(y_range_name="Te"), "right")

    electron_max_vs_timestamp.line(
        "x",
        "y",
        source=sources_time["Te"],
        legend_label="Te",
        color="red",
        y_range_name="Te",
    )

    scale_plot_radial(electron_max_vs_timestamp)


def plot_ne_te(ne_ee_plot, RADIAL_POSITIONS, sources_rad):
    sources_rad["ne"] = ColumnDataSource(
        data=dict(x=RADIAL_POSITIONS, y=[0] * len(RADIAL_POSITIONS))
    )
    sources_rad["Te"] = ColumnDataSource(
        data=dict(x=RADIAL_POSITIONS, y=[0] * len(RADIAL_POSITIONS))
    )

    ne_ee_plot.extra_y_ranges = {"Te": Range1d(start=0, end=1)}

    ne_ee_plot.line("x", "y", source=sources_rad["ne"], legend_label="ne", color="blue")

    ne_ee_plot.add_layout(LinearAxis(y_range_name="Te"), "right")
    ne_ee_plot.line(
        "x",
        "y",
        source=sources_rad["Te"],
        legend_label="Te",
        color="red",
        y_range_name="Te",
    )
    scale_plot_time(te_plot=ne_ee_plot)


def plot_columns(
    ne_ee_plot,
    electron_max_vs_timestamp,
    concentration_plot,
    energy_plot,
    concentration_plot_time,
    energy_plot_time,
    RADIAL_POSITIONS,
    sources_rad,
    sources_time,
):
    global df
    counter_concentration = 0
    counter_energy = 0

    plot_ne_time(
        electron_max_vs_timestamp, df["tmain"].unique().tolist(), sources_time
    )
    plot_ne_te(ne_ee_plot, RADIAL_POSITIONS, sources_rad)

    for column in df.columns:
        if column == "tmain" or column == "RadialPositions" or column == "ne":
            continue

        if should_plot_concentration(column):
            plot_concentration(
                concentration_plot,
                column,
                RADIAL_POSITIONS,
                counter_concentration,
                sources_rad,
            )
            plot_concentration_time(
                concentration_plot_time,
                column,
                counter_concentration,
                sources_time,
            )
            counter_concentration += 1

        elif should_plot_temperature(column):
            plot_energy(
                energy_plot, column, RADIAL_POSITIONS, counter_energy, sources_rad
            )
            plot_energy_time(energy_plot_time, column, counter_energy, sources_time)
            counter_energy += 1


def initialize_interactive_elements(ne_ee_plot, electron_max_vs_timestamp):
    # Initialize the 'ne' source first
    global sources_rad, electron_max_vs_timestamp_line
    source_ne = sources_rad["ne"]
    source_ne_timestamp = sources_time["ne"]
    ne_ee_plot.line("x", "y", source=source_ne, legend_label="ne", color="blue")

    # Add a circle glyph for the movable point
    movable_point_source = ColumnDataSource(data=dict(x=[90], y=[0]))
    ne_ee_plot.circle(
        "x",
        "y",
        source=movable_point_source,
        size=5,
        color="black",
        name="movable_point",
    )

    movable_point_source_timestamp = ColumnDataSource(data=dict(x=[0], y=[0]))
    electron_max_vs_timestamp.circle(
        "x",
        "y",
        source=movable_point_source_timestamp,
        size=5,
        color="black",
        name="movable_point_timestamp",
    )

    selected_radius_source = ColumnDataSource(data=dict(radius=[0]))

    # Add a vertical line (span) at the position of the movable point
    vertical_line = Span(
        location=90,
        dimension="height",
        line_color="black",
        line_dash="dashed",
        line_width=1,
    )
    ne_ee_plot.add_layout(vertical_line)

    vertical_line_timestamp = Span(
        location=0,
        dimension="height",
        line_color="black",
        line_dash="dashed",
        line_width=1,
    )
    electron_max_vs_timestamp.add_layout(vertical_line_timestamp)

    # Add a label to display the radius position
    label = Label(x=0, y=0, text="", text_font_size="10pt", text_color="black")
    ne_ee_plot.add_layout(label)

    label_timestamp = Label(
        x=0, y=0, text="", text_font_size="10pt", text_color="black"
    )
    electron_max_vs_timestamp.add_layout(label_timestamp)

    if electron_max_vs_timestamp_line is None:
        electron_max_vs_timestamp_line = electron_max_vs_timestamp.line(
            x=[], y=[], legend_label="ne", color="blue"
        )

    return (
        movable_point_source,
        vertical_line,
        label,
        source_ne,
        source_ne_timestamp,
        selected_radius_source,
        movable_point_source_timestamp,
        vertical_line_timestamp,
        label_timestamp,
    )


def modify_doc(doc):
    global df
    df = pd.read_csv(DATA_FILE)
    RADIAL_POSITIONS = df["RadialPositions"].unique().tolist()

    (
        ne_ee_plot,
        concentration_plot,
        energy_plot,
        electron_max_vs_timestamp,
        concentration_plot_time,
        energy_plot_time,
    ) = initialize_plots(RADIAL_POSITIONS)

    global electron_max_vs_timestamp_line, sources_rad, sources_time, timestamp_trigger_source, selected_radius, selected_time

    calculate_temperature(df)

    plot_columns(
        ne_ee_plot,
        electron_max_vs_timestamp,
        concentration_plot,
        energy_plot,
        concentration_plot_time,
        energy_plot_time,
        RADIAL_POSITIONS,
        sources_rad,
        sources_time,
    )

    (
        movable_point_source,
        vertical_line,
        label,
        source_ne,
        source_ne_timestamp,
        selected_radius_source,
        movable_point_source_timestamp,
        vertical_line_timestamp,
        label_timestamp,
    ) = initialize_interactive_elements(ne_ee_plot, electron_max_vs_timestamp)

    # Define a JavaScript callback to update the movable point position
    callback = CustomJS(
        args=dict(
            source=movable_point_source,
            line=vertical_line,
            label=label,
            ne_source=source_ne,
            selected_radius_source=selected_radius_source,
        ),
        code="""
        const data = source.data;
        const ne_data = ne_source.data;
        const x = cb_obj.x;

        // Find the nearest x-value in the ne_data and get the corresponding y-value
        let closest_distance = Infinity;
        let closest_x = 0; // This will hold the closest discrete x value
        let closest_y = 0;
        for (let i = 0; i < ne_data['x'].length; i++) {
            const distance = Math.abs(ne_data['x'][i] - x);
            if (distance < closest_distance) {
                closest_distance = distance;
                closest_x = ne_data['x'][i]; // Store the exact discrete x value
                closest_y = ne_data['y'][i];
            }
        }

        data['x'][0] = closest_x;
        data['y'][0] = closest_y;  // Set y to the value from the ne curve
        line.location = closest_x;
        label.x = closest_x;
        label.y = closest_y + 0.05 * closest_y;  // Adjust this for label positioning
        label.text = "Radius: " + closest_x.toFixed(2);
        selected_radius_source.data = {radius: [closest_x]};
        selected_radius_source.change.emit();
        source.change.emit();
    """,
    )

    def update_title(attr, old, new):
        global selected_radius, df, selected_time
        global electron_max_vs_timestamp_line
        selected_radius = selected_radius_source.data["radius"][0]
        selected_radius = get_closest_value(RADIAL_POSITIONS, selected_radius)
        

        timestamps = df[df["RadialPositions"] == selected_radius]["tmain"].tolist()

        for column in sources_time:
            sources_time[column].data = {
                "x": timestamps,
                "y": df[df["RadialPositions"] == selected_radius][column].tolist(),
            }

        electron_max_vs_timestamp.title.text = (
            f"ne & Te vs Time Stamp at Radius: {selected_radius:.2f}"
        )

        concentration_plot_time.title.text = (
            f"Concentration vs Time Stamp at Radius: {selected_radius:.2f}"
        )

        energy_plot_time.title.text = (
            f"Temperature vs Time Stamp at Radius: {selected_radius:.2f}"
        )
        scale_plot_radial(electron_max_vs_timestamp, selected_radius)
        scale_plot_time(ne_ee_plot, selected_time)
        
        

    # Set the callback for the entire figure to update the movable point on any click
    ne_ee_plot.js_on_event("tap", callback)
    selected_radius_source.on_change("data", update_title)

    selected_timestamp_source = ColumnDataSource(data=dict(timestamp=[0]))
        

    callback_timestamp = CustomJS(
        args=dict(
            source=movable_point_source_timestamp,
            line=vertical_line_timestamp,
            label=label_timestamp,
            time_source=source_ne_timestamp,
            selected_timestamp_source=selected_timestamp_source,
            timestamp_trigger_source=timestamp_trigger_source,
        ),
        code="""
            function updateDisplay(x) {
                console.log(x)
                const data = source.data;
                const timestamp_data = time_source.data;

                let closest_distance = Infinity;
                let closest_x = 0;
                let closest_y = 0;
                for (let i = 0; i < timestamp_data['x'].length; i++) {
                    const distance = Math.abs(timestamp_data['x'][i] - x);
                    if (distance < closest_distance) {
                        closest_distance = distance;
                        closest_x = timestamp_data['x'][i];
                        closest_y = timestamp_data['y'][i];
                    }
                }

                // Use the closest_x value for positioning
                data['x'][0] = closest_x;
                data['y'][0] = closest_y;
                line.location = closest_x;
                label.x = closest_x;
                label.y = closest_y + 0.01 * closest_y;  // Adjust this for label positioning
                label.text = "Timestamp: " + closest_x.toFixed(2); // Display the exact discrete x value
                selected_timestamp_source.data = {timestamp: [closest_x]};
                selected_timestamp_source.change.emit();
                source.change.emit();
                
            }
            
            
            // This part will be called by the tap event
            if (cb_obj && cb_obj.x) {
                updateDisplay(cb_obj.x);
            }
            else {
                updateDisplay(timestamp_trigger_source.data['last_tmain'][0]);
            }
        """,
    )


    def print_timestamp(attr, old, new):
        global x_value, df, selected_time, selected_radius
        x_value = selected_timestamp_source.data["timestamp"][0]
        selected_time = x_value
        df_last = df[df["tmain"] == selected_time]
        
        

        # Update the title of each plot to include the current timestamp
        ne_ee_plot.title.text = (
            f"ne & Te vs Radial Position at {round(selected_time*1000)}ms"
        )
        concentration_plot.title.text = (
            f"Concentration vs Radial Position at {round(selected_time*1000)}ms"
        )
        energy_plot.title.text = (
            f"Temperature vs Radial Position at {round(selected_time*1000)}ms"
        )

        # Update sources_rad
        for column in sources_rad:
            sources_rad[column].data = dict(
                x=RADIAL_POSITIONS, y=df_last[column].tolist()
            )
        scale_plot_time(ne_ee_plot, selected_time)
        scale_plot_radial(electron_max_vs_timestamp, selected_radius)
        
        
            
    timestamp_trigger_source.js_on_change('data', callback_timestamp)
    electron_max_vs_timestamp.js_on_event("tap", callback_timestamp)
    selected_timestamp_source.on_change("data", print_timestamp)

    global timestamp_max_concentration
    timestamp_max_concentration = {}

    doc.add_periodic_callback(
        lambda: update_data(
            RADIAL_POSITIONS,
            sources_rad,
            ne_ee_plot,
            concentration_plot,
            energy_plot,
            concentration_plot_time,
            energy_plot_time,
            electron_max_vs_timestamp,
        ),
        5000,
    )
    arrange_grid(
        doc,
        ne_ee_plot,
        concentration_plot,
        energy_plot,
        electron_max_vs_timestamp,
        concentration_plot_time,
        energy_plot_time,
    )


def scale_plot_time(te_plot, time=None):
    #print(f"Scaling plot for time: {time}")
    original_df = pd.read_csv(DATA_FILE)
    
    if time is None:
        time = original_df['tmain'].max()
    
    latest_df = original_df[original_df['tmain'] == time].copy()
    latest_df["Te"] = (2 / 3) * (latest_df["Ee"] / latest_df["ne"])
    
    max_ne = latest_df["ne"].max()
    max_Te = latest_df["Te"].max()
    #print(f"Max ne: {max_ne}, Max Te: {max_Te}")
    
    te_plot.y_range.start = 0 - max_ne * 0.1
    te_plot.y_range.end = max_ne * 1.1
    
    te_plot.extra_y_ranges["Te"].start = 0 - max_Te * 0.1
    te_plot.extra_y_ranges["Te"].end = max_Te * 1.1
    
    min_radial = latest_df["RadialPositions"].min()
    max_radial = latest_df["RadialPositions"].max()
    
    te_plot.x_range.start = min_radial - 0.1 * (max_radial - min_radial)
    te_plot.x_range.end = max_radial + 0.1 * (max_radial - min_radial)

def scale_plot_radial(te_plot, radius=90.0):
    #print(f"Scaling plot for radius: {radius}")
    original_df = pd.read_csv(DATA_FILE)
    
    # Find closest radius in the DataFrame to the selected radius
    closest_radius = get_closest_value(original_df["RadialPositions"].unique(), radius)
    
    radius_df = original_df[original_df["RadialPositions"] == closest_radius].copy()
    radius_df["Te"] = (2 / 3) * (radius_df["Ee"] / radius_df["ne"])
    
    max_ne = radius_df["ne"].max()
    max_Te = radius_df["Te"].max()
    #print(f"Max ne: {max_ne}, Max Te: {max_Te}")
    
    te_plot.y_range.start = 0 - max_ne * 0.1
    te_plot.y_range.end = max_ne * 1.1
    
    te_plot.extra_y_ranges["Te"].start = 0
    te_plot.extra_y_ranges["Te"].end = max_Te * 1.1
    
    max_time = original_df['tmain'].max()
    
    # Assuming `time` is available in the context or passed as a parameter
    te_plot.x_range.start = 0 - 0.1 * max_time
    te_plot.x_range.end = max_time + 0.1 * max_time 


    


def initialize_plots(radial_positions):
    # Create and return the initialized plots
    ne_ee_plot = figure(
        width=PLOT_WIDTH,
        height=PLOT_HEIGHT,
        title="ne & Te vs Radial Position",
        y_axis_label="ne",
    )

    # Apply log scale to y-axis for concentration_plot
    concentration_plot = figure(
        width=PLOT_WIDTH, 
        height=PLOT_HEIGHT, 
        title="Concentration vs Radial Position",
        y_scale=LogScale(),  # Set the y-axis to logarithmic scale
    )
    
    energy_plot = figure(
        width=PLOT_WIDTH,
        height=PLOT_HEIGHT,
        title="Temperature vs Radial Position"
    )

    electron_max_vs_timestamp = figure(
        width=PLOT_WIDTH,
        height=PLOT_HEIGHT,
        title="Electron Concentration vs Time Stamp",
        x_axis_label="Time Stamp",
        y_axis_label="ne",
    )

    # Apply log scale to y-axis for concentration_plot_time
    concentration_plot_time = figure(
        width=PLOT_WIDTH,
        height=PLOT_HEIGHT,
        title="Concentration vs Time Stamp",
        x_axis_label="Time Stamp",
        y_scale=LogScale(),  # Set the y-axis to logarithmic scale
    )

    energy_plot_time = figure(
        width=PLOT_WIDTH,
        height=PLOT_HEIGHT,
        title="Temperature vs Time Stamp",
        x_axis_label="Time Stamp",
    )
    
    concentration_plot.y_scale = LogScale()
    concentration_plot_time.y_scale = LogScale()

    return (
        ne_ee_plot,
        concentration_plot,
        energy_plot,
        electron_max_vs_timestamp,
        concentration_plot_time,
        energy_plot_time,
    )



def should_plot_concentration(column):
    undesired_subs = ["d", "u", "C", "tnew", "minstopcrit"]
    return "n" in column and not any(sub in column for sub in undesired_subs)


def is_energy_column(column_name):
    """Check if a column name is for energy."""
    return "E" in column_name and "d" not in column_name


def should_plot_temperature(column):
    return column.startswith("T") and column != "Te"


def update_timestamp_max_concentration():
    global sources_time, selected_radius, df
    timestamps = df[df["RadialPositions"] == selected_radius]["tmain"].tolist()

    for column in sources_time:
        sources_time[column].data = {
            "x": timestamps,
            "y": df[df["RadialPositions"] == selected_radius][column].tolist(),
        }


def update_data(
    radial_positions,
    sources_rad,
    ne_ee_plot,
    concentration_plot,
    energy_plot,
    concentration_plot_time,
    energy_plot_time,
    electron_max_vs_timestamp,
):
    global DATA_FILE, timestamp_max_concentration, selected_radius, electron_max_vs_timestamp_line, last_modified_time, timestamp_trigger_source

    current_modified_time = os.path.getmtime(DATA_FILE)

    if last_modified_time == current_modified_time:
        return
    else:
        last_modified_time = current_modified_time

    if selected_radius is None:
        selected_radius = get_closest_value(radial_positions, 90)

    df = pd.read_csv(DATA_FILE)

    calculate_temperature(df)

    update_timestamp_max_concentration()

    # Get last value of tmain
    last_tmain = df["tmain"].iloc[-1]
    df_last = df[df["tmain"] == last_tmain]
    

    # Update the title of each plot to include the current timestamp
    ne_ee_plot.title.text = f"ne & Te vs Radial Position at {round(last_tmain*1000)}ms"
    concentration_plot.title.text = (
        f"Concentration vs Radial Position at {round(last_tmain*1000)}ms"
    )
    energy_plot.title.text = (
        f"Temperature vs Radial Position at {round(last_tmain*1000)}ms"
    )

    electron_max_vs_timestamp.title.text = (
        f"ne & Te vs Time Stamp at Radius: {selected_radius:.2f}"
    )

    concentration_plot_time.title.text = (
        f"Concentration vs Time Stamp at Radius: {selected_radius:.2f}"
    )

    energy_plot_time.title.text = (
        f"Temperature vs Time Stamp at Radius: {selected_radius:.2f}"
    )

    # Update sources_rad
    for column in sources_rad:
        sources_rad[column].data = dict(x=radial_positions, y=df_last[column].tolist())

    # Update sources_time
    for column in sources_time:
        sources_time[column].data = {
            "x": df["tmain"].unique().tolist(),
            "y": df[df["RadialPositions"] == selected_radius][
                column
            ].tolist(),
        }
        
    timestamp_trigger_source.data = {'last_tmain': [last_tmain]}


def arrange_grid(
    doc,
    ne_ee_plot,
    concentration_plot,
    energy_plot,
    electron_max_vs_timestamp,
    concentration_plot_time,
    energy_plot_time,
):
    grid = [
        [ne_ee_plot, concentration_plot, energy_plot],
        [electron_max_vs_timestamp, concentration_plot_time, energy_plot_time],
    ]
    grid = gridplot(grid)
    doc.add_root(grid)


def start_plot_app():
    curdoc().title = "Tomator"
    wait_for_data_file()
    modify_doc(curdoc())


def wait_for_data_file():
    while os.path.getsize(DATA_FILE) == 0:
        time.sleep(1)


with open("InterfaceFuncs/aux.txt", "r") as file:
    DATA_FILE = file.readline().strip()

printBlue(f"Updating according to: {DATA_FILE}")
start_plot_app()
