import tkinter as tk
from tkinter import ttk, filedialog, messagebox, font, scrolledtext
import json


class SimulationParamsApp(tk.Toplevel):
    def __init__(self, master=None):
        super().__init__(master)

        self.DEFAULT_WIDTH = 800
        self.DEFAULT_HEIGHT = 500
        self.DEFAULT_X = self.winfo_screenwidth() // 2 - self.DEFAULT_WIDTH // 2
        self.DEFAULT_Y = self.winfo_screenheight() // 2 - self.DEFAULT_HEIGHT // 2
        self.data = {}

        self.option_widget_map = {}  # Map to store widgets for each type

        self.type_to_parameters_mapping = { # Mapping of parameters to hide for each type
            "bgray": [],
            "bram": [],
            "bfixpowerfrac": [],
            "bnefix": [],
            "bproptone": ["Rdep1", "PRdep1", "Rdep2", "PRdep2", "Rdep3", "PRdep3", "Rdep4", "PRdep4"],
            "bTOMAS": [],
            "bnopower": [],
            "bkipt": [], #"bscaling": [],
            "btunedv": [],
            "bmanuel": [],
        }

        self.title("Simulation Parameters")
        self.after(1, self.setup_geometry_and_ui)

    def setup_geometry_and_ui(self):
        style = ttk.Style()
        style.configure(
            "Centered.TButton", padding=(0, 0, 0, 0)
        )  # (left, top, right, bottom)

        self.option_add("*Font", "Futura")

        style.configure("TLabel", font=("Futura"))
        style.configure("TEntry", font=("Futura"))
        style.configure("TButton", font=("Futura"))
        style.configure("TCheckbutton", font=("Futura"))

        self.geometry(
            f"{self.DEFAULT_WIDTH}x{self.DEFAULT_HEIGHT}+{self.DEFAULT_X}+{self.DEFAULT_Y}"
        )

        self.title("Simulation Parameters")

        self.saved_file_name = None

        # Treeview for organized frames
        self.tree = ttk.Treeview(self)
        self.tree.column("#0", width=300)  # Set the column width
        self.tree.grid(row=0, column=0, padx=10, pady=10, sticky="ns")

        # Frame list for easy access
        self.frames = {}

        # Main group: MAIN DISCHARGE PARAMETERS
        discharge_group = self.tree.insert("", "end", text="Main discharge", open=False)

        self.add_frame(
            discharge_group,
            "Magnetic field",
            [
                ("Bt (T)", tk.DoubleVar(value=format(1.5, ".5e"))),
                ("Bv (T)", tk.DoubleVar(value=format(0.0075, ".5e"))),
                ("Bh (T)", tk.DoubleVar(value=format(0.0001, ".5e"))),
            ],
        )

        self.add_frame(
            discharge_group,
            "Toroidal machine geometry",
            [
                ("R (cm)", tk.DoubleVar(value=format(88.0, ".5e"))),
                ("a (cm)", tk.DoubleVar(value=format(28.0, ".5e"))),
                ("b (cm)", tk.DoubleVar(value=format(75.0, ".5e"))),
                ("lHFS (cm)", tk.DoubleVar(value=format(80.5808, ".5e"))),
                ("lLFS (cm)", tk.DoubleVar(value=format(95.4192, ".5e"))),
                ("nlimiters", tk.IntVar(value=6)),
                ("Vpl (cm3)", tk.DoubleVar(value=format(56105.104899, ".7e"))),
            ],
        )

        self.add_frame(
            discharge_group,
            "Neutral pressure",
            [
                ("pHe (Pa)", tk.DoubleVar(value=format(0.000168, ".5e"))),
                ("pH2 (Pa)", tk.DoubleVar(value=format(0.0000432, ".5e"))),
            ],
        )

        # Main group: COUPLED (RF) POWER
        rf_power_group = self.tree.insert(
            "", "end", text="Coupled (RF) power", open=False
        )
        
        self.add_frame(
            rf_power_group,
            "RF power",
            [
                ("Prf (kW)", tk.DoubleVar(value=format(400.0, ".5e"))),
                ("bselfcol", tk.BooleanVar(value=False)),
                ("freq (MHz)", tk.DoubleVar(value=format(82700.0, ".5e"))),
                ("dtpramp (ms)", tk.DoubleVar(value=format(0.0, ".5e"))),
            ],
        )
        
        self.add_frame(
            rf_power_group,
            "Type",
            [
                ("bgray", tk.BooleanVar(value=False)),
                ("bram", tk.BooleanVar(value=False)),
                ("bfixpowerfrac", tk.BooleanVar(value=False)),
                ("bnefix", tk.BooleanVar(value=True)),
                ("bTOMAS", tk.BooleanVar(value=False)),
                ("bproptone", tk.BooleanVar(value=False)),
                ("bnopower", tk.BooleanVar(value=False)),
                ("bkipt", tk.BooleanVar(value=False)),
                ("blhr", tk.BooleanVar(value=False)),
                ("bmanuel", tk.BooleanVar(value=False)),
                ("bICWC", tk.BooleanVar(value=False)),
            ],
        )

        # Create a 'Types' subgroup under 'rf_power_group'
        types_group = self.tree.insert(rf_power_group, "end", text="ECH", open=False)
        
        self.add_frame(
            types_group,
            "general ec",
            [
                ("Rdep (cm)", tk.DoubleVar(value=format(89.4, ".5e"))),
                ("pecabs0", tk.DoubleVar(value=format(0.1, ".5e"))),
                ("widthech", tk.DoubleVar(value=format(1.8666666666666667, ".5e"))), # UNITS?
                ("echbackground", tk.DoubleVar(value=format(0.0000001, ".5e"))), # UNITS?
                ("harmonic", tk.DoubleVar(value=format(2.0, ".5e"))), # UNITS?
                ("muw", tk.DoubleVar(value=format(0.01, ".5e"))), # UNITS?
            ],
        )

        self.add_frame(
            types_group,
            "GRAY",
            [],
        )

        self.add_frame(
            types_group,
            "RAM",
            [],
        )

        self.add_frame(
            types_group,
            "necfix",
            [
                ("ic", tk.IntVar(value=86)),
                ("necfix (m-3)", tk.DoubleVar(value=format(10584000000000.0, ".5e"))),
                ("Pini", tk.DoubleVar(value=format(3.0, ".5e"))), # UNITS?
                ("tauP", tk.DoubleVar(value=format(0.00005, ".5e"))), # UNITS?
            ],
        )

        self.add_frame(
            types_group,
            "TOMAS",
            [
                ("Rdep1 (cm)", tk.DoubleVar(value=format(86.1875, ".5e"))),
                ("PRdep1", tk.DoubleVar(value=format(0.0, ".5e"))),
                ("Rdep2 (cm)", tk.DoubleVar(value=format(97.65, ".5e"))),
                ("PRdep2", tk.DoubleVar(value=format(1.0, ".5e"))),
                ("Rdep3 (cm)", tk.DoubleVar(value=format(72.4325, ".5e"))),
                ("PRdep3", tk.DoubleVar(value=format(0.0, ".5e"))),
                ("Rdep4 (cm)", tk.DoubleVar(value=format(72.4325, ".5e"))),
                ("PRdep4", tk.DoubleVar(value=format(0.0, ".5e"))),
            ],
        )

        # Create a 'Types' subgroup under 'rf_power_group'
        types_group = self.tree.insert(rf_power_group, "end", text="ICH", open=False)

        self.add_frame(
            types_group,
            "general ic",
            [
                ("Rant (cm)", tk.DoubleVar(value=format(200.0, ".5e"))),
                ("bantlr", tk.BooleanVar(value=False)),
                ("avlr", tk.DoubleVar(value=format(0.23, ".5e"))), # UNITS?
                ("alphalaw", tk.DoubleVar(value=format(0.002, ".5e"))), # UNITS?
                ("HtoHD", tk.DoubleVar(value=format(1.0, ".5e"))), # UNITS?
            ],
        )

        self.add_frame(
            types_group,
            "KIPT",
            [],
        )

        self.add_frame(
            types_group,
            "LHR",
            [
                ("fracpne", tk.DoubleVar(value=format(1.0, ".5e"))),
                ("fraclhr", tk.DoubleVar(value=format(0.25, ".5e"))),
                ("widthlhr (cm)", tk.DoubleVar(value=format(2.5, ".5e"))),
                ("lhrbackground", tk.DoubleVar(value=format(1e-7, ".5e"))),
            ],
        )

        self.add_frame(
            types_group,
            "ERM",
            [],
        )

        self.add_frame(
            types_group,
            "PPPL",
            [],
        )

        self.add_frame(
            rf_power_group,
            "OTHER",
            [
                ("fixedTe", tk.BooleanVar(value=False)),
            ],
        )

        #Main group: TRANSPORT PARAMETERS
        transport_group = self.tree.insert("", "end", text="Transport", open=False)

        self.add_frame(
            transport_group,
            "Diffusion",
            [
                ("bDfix", tk.BooleanVar(value=False)),  # to be added
                ("Dfix (cm2/s)", tk.DoubleVar(value=format(10000.0, ".5e"))),  # to be added
                ("bDbohm", tk.BooleanVar(value=False)), # to be added
                ("bDscaling", tk.BooleanVar(value=True)), # to be added
                ("Dfact", tk.DoubleVar(value=format(1.0, ".5e"))),
            ]
        )

        self.add_frame(
            transport_group,
            "Convection",
            [
                ("bVfix", tk.BooleanVar(value=False)),      # added ?
                ("Vfix (cm/s)", tk.DoubleVar(value=format(100.0, ".5e"))),        # added ?
                ("bVscaling", tk.BooleanVar(value=True)),  # to be added
                ("veq", tk.IntVar(value=8)),
                ("Vfact", tk.DoubleVar(value=format(5.0, ".5e"))),
            ]
        )

        self.add_frame(
            transport_group,
            "Tune D and V",
            [
                ("btunedv", tk.BooleanVar(value=True)),
                ("btunevleft", tk.BooleanVar(value=True)),
                ("il", tk.IntVar(value=73)),
                ("nelfix (m-3)", tk.DoubleVar(value=format(3189800000000.0, ".5e"))),
                ("ir", tk.IntVar(value=114)),
                ("nerfix (m-3)", tk.DoubleVar(value=format(6430100000000.0, ".5e"))),
                ("Vini", tk.DoubleVar(value=format(-2.0, ".5e"))), # UNITS?
                ("tauV", tk.DoubleVar(value=format(0.01, ".5e"))), # UNITS?    
                ("Dini", tk.DoubleVar(value=format(1.0, ".5e"))), # UNITS?
                ("tauD", tk.DoubleVar(value=format(0.0005, ".5e"))), # UNITS?
            ],
        )

        # Main group: CONDITIONS
        conditions_group = self.tree.insert("", "end", text="Conditions", open=False)



        self.add_frame(
            conditions_group,
            "Physics to include",
            [
                ("bH", tk.BooleanVar(value=True)),
                ("bH2", tk.BooleanVar(value=True)),
                ("bHe", tk.BooleanVar(value=True)),
                ("bion", tk.BooleanVar(value=True)),
                ("bcx", tk.BooleanVar(value=True)),
                ("belas", tk.BooleanVar(value=True)),
                ("bcoulomb", tk.BooleanVar(value=True)),
                ("bimpur", tk.BooleanVar(value=False)),
                ("btranspions", tk.BooleanVar(value=True)),
                ("btranspneut", tk.BooleanVar(value=True)),
                ("bedge", tk.BooleanVar(value=True)),
                ("bpol", tk.BooleanVar(value=True)),
                ("vdrift", tk.BooleanVar(value=True)),
                ("bcoll", tk.BooleanVar(value=True)),
            ],
        )

        self.add_frame(
            conditions_group,
            "Initial conditions",
            [
                ("rmaxini (cm)", tk.DoubleVar(value=format(89.4, ".5e"))),
                ("widthini (cm)", tk.DoubleVar(value=format(2.8, ".5e"))),
                ("nebackgroundl (cm-3)", tk.DoubleVar(value=format(0.00001, ".5e"))), # UNITS?
                ("nebackgroundr (cm-3)", tk.DoubleVar(value=format(0.001, ".5e"))), # UNITS?
                ("Ta0 (eV)", tk.DoubleVar(value=format(0.02587, ".5e"))),
                ("Te0 (eV)", tk.DoubleVar(value=format(3.0, ".5e"))),
                ("nevac", tk.DoubleVar(value=format(1.0, ".5e"))),
                ("nH0 (m-3)", tk.DoubleVar(value=format(2000000.0, ".5e"))),
                ("nHi0 (m-3)", tk.DoubleVar(value=format(1000000.0, ".5e"))),
                ("nH2i0 (m-3)", tk.DoubleVar(value=format(1000000.0, ".5e"))),
                ("nH3i0 (m-3)", tk.DoubleVar(value=format(2000000.0, ".5e"))),
                ("nHeII0 (m-3)", tk.DoubleVar(value=format(10053600000000.0, ".5e"))),
                ("nHeIII0 (m-3)", tk.DoubleVar(value=format(2000000.0, ".5e"))),
                ("nCII0 (m-3)", tk.DoubleVar(value=format(1000000000.0, ".5e"))),
                ("nCIII0 (m-3)", tk.DoubleVar(value=format(10.0, ".5e"))),
                ("nCIV0 (m-3)", tk.DoubleVar(value=format(10.0, ".5e"))),
                ("nCV0 (m-3)", tk.DoubleVar(value=format(10.0, ".5e"))),
            ],
        )

        self.add_frame(
            conditions_group,
            "Edge conditions",
            [
                ("RH", tk.DoubleVar(value=format(0.5, ".5e"))), # UNITS?
                ("REH", tk.DoubleVar(value=format(0.9, ".5e"))),   # UNITS?
                ("gEd", tk.DoubleVar(value=format(2.5, ".5e"))), # UNITS?
                ("gEv", tk.DoubleVar(value=format(1.5, ".5e"))), # UNITS?
                ("gEdn", tk.DoubleVar(value=format(2.0, ".5e"))), # UNITS?
                ("gEe", tk.DoubleVar(value=format(2.5, ".5e"))), # UNITS?
            ],
        )

        # Main group: SIMULATION CONTROL SETTINGS
        control_group = self.tree.insert(
            "", "end", text="Simulation control settings", open=False
        )

        self.add_frame(
            control_group,
            "Simulation grid",
            [
                ("nmeshp", tk.IntVar(value=301)),
            ],
        )
        self.add_frame(
            control_group,
            "Input file",
            [
                ("bfinput", tk.BooleanVar(value=False)),
            ],
        )

        self.add_frame(
            control_group,
            "Time step",
            [
                ("t0 (s)", tk.DoubleVar(value=format(0.0, ".5e"))),
                ("tmainend (s)", tk.DoubleVar(value=format(1.0, ".5e"))),
                ("accur (s)", tk.DoubleVar(value=format(0.05, ".5e"))),
                ("dtmax (s)", tk.DoubleVar(value=format(0.000002, ".5e"))),
                ("dtmin (s)", tk.DoubleVar(value=format(0.000000025, ".5e"))),
                ("dtinit (s)", tk.DoubleVar(value=format(0.00000000001, ".5e"))),
            ],
        )

        self.add_frame(
            control_group,
            "Time step for RF coupling",
            [
                ("bupdateRFstep", tk.BooleanVar(value=True)),
                ("updateRF", tk.IntVar(value=1)),
                ("dtRF (s)", tk.DoubleVar(value=format(0.0000001, ".5e"))),
                ("dtRFvar (s)", tk.BooleanVar(value=False)),
                ("dtRFmax (s)", tk.DoubleVar(value=format(0.000005, ".5e"))),
                ("dtRFmin (s)", tk.DoubleVar(value=format(0.00000001, ".5e"))),
                ("dtRFconv (s)", tk.BooleanVar(value=True)),
            ],
        )

        self.add_frame(
            control_group,
            "Advanced time step settings",
            [
                ("dtsmooth", tk.BooleanVar(value=True)),
                ("shokparam", tk.DoubleVar(value=format(2.0, ".5e"))), # UNITS?
                ("maxtstepincrement", tk.DoubleVar(value=format(1.5, ".5e"))), # UNITS?
                ("btendvar", tk.BooleanVar(value=False)),
                ("convcrit", tk.DoubleVar(value=format(100.0, ".5e"))), # UNITS?
                ("convsavetime (s)", tk.DoubleVar(value=format(0.01, ".5e"))),
                ("baccurvar", tk.BooleanVar(value=True)),
                ("accurcrit", tk.DoubleVar(value=format(500, ".5e"))), # UNITS?
                ("minaccur", tk.DoubleVar(value=format(0.05, ".5e"))), # UNITS?
            ],
        )

        self.add_frame(
            control_group,
            "Output parameters",
            [
                ("Nlog", tk.IntVar(value=500)),
                ("cc", tk.IntVar(value=80)),
                ("Nloopsave", tk.IntVar(value=5000)),
                ("bOutdt", tk.BooleanVar(value=True)),
                ("dtsave (s)", tk.DoubleVar(value=format(0.001, ".5e"))),
            ],
        )

        self.add_frame(
            control_group,
            "Solver parameters",
            [
                ("solvertolerance", tk.DoubleVar(value=format(1e-10, ".5e"))), # UNITS?
            ],
        )

        # Bind a double click event to the treeview
        self.tree.bind("<ButtonRelease-1>", self.on_item_click)
        

        # Calculate the position for the button to be centered at the bottom half of the window
        button_width = 150  # or any desired width
        button_height = 30  # or any desired height
        x_position = (self.DEFAULT_WIDTH - button_width) / 2
        y_position = self.DEFAULT_HEIGHT / 2 + (
            self.DEFAULT_HEIGHT / 4 - button_height / 2
        )
        
        # Load button
        self.load_defaults_btn = ttk.Button(
            self,
            text="Load Defaults",
            command=self.load_defaults,
            style="Centered.TButton",
        )
        self.load_defaults_btn.place(
            x=x_position, y=y_position - button_height - 10, width=button_width, height=button_height
        )
        
        # Save button
        self.save_btn = ttk.Button(
            self,
            text="Save to JSON",
            command=self.save_to_json,
            style="Centered.TButton",
        )
        self.save_btn.place(
            x=x_position, y=y_position, width=button_width, height=button_height
        )

    def add_frame(self, parent, name, parameters):
        canvas = tk.Canvas(self)
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        frame = ttk.Frame(canvas)

        type_variable = tk.StringVar()  # Use a single variable for all radio buttons in a group
        type_variable.trace('w', lambda *args: self.on_type_selected(type_variable))  # Trace changes to this variable

        row_num = 1  # To control the grid placement

        # Add widgets to the frame
        for text, var in parameters:
            label = ttk.Label(frame, text=text)
            label.grid(row=row_num, column=0, sticky="w", padx=5, pady=5)

            if isinstance(var, tk.BooleanVar):
                if name == "Type":
                    rad = ttk.Radiobutton(frame, text="", variable=type_variable, value=text, command=lambda: self.on_type_selected(type_variable, parameters))
                    rad.grid(row=row_num, column=1, padx=5, pady=5, sticky="w")
                elif name == "Transport Type":
                    rad = ttk.Radiobutton(frame, text="", variable=type_variable, value=text, command=lambda: self.on_type_selected(type_variable, parameters))
                    rad.grid(row=row_num, column=1, padx=5, pady=5, sticky="w")
                else:
                    # Use checkboxes for all other categories
                    chk = ttk.Checkbutton(frame, text="", variable=var)
                    chk.grid(row=row_num, column=1, padx=5, pady=5, sticky="w")
                    self.option_widget_map[text] = [label, chk]
                    
                if name == "Input file":  
                    style = ttk.Style()
                    style.configure("Centered.TButton", padding=(0, 0, 0, 0), font=("Futura", 10))
                    action_button = ttk.Button(frame, text="Select Input File", command=self.on_button_pressed, width=15, style="Centered.TButton")
                    action_button.grid(row=row_num + 1, column=0, padx=2, pady=5, sticky="w")
                
                # Check if in OTHER and bmanuel is True
                if name == "OTHER":
                    style = ttk.Style()
                    style.configure("Centered.TButton", padding=(0, 0, 0, 0), font=("Futura", 10))
                    action_button = ttk.Button(frame, text="Select Input File", command=self.on_button_pressed_2, width=15, style="Centered.TButton")
                    action_button.grid(row=row_num + 1, column=0, padx=2, pady=5, sticky="w")


            else:
                label = ttk.Label(frame, text=text)
                entry = ttk.Entry(frame, textvariable=var)
                label.grid(row=row_num, column=0, sticky="w", padx=5, pady=5)
                entry.grid(row=row_num, column=1, padx=5, pady=5, sticky="w")

                if text not in self.option_widget_map:
                    self.option_widget_map[text] = []
                self.option_widget_map[text].extend([label, entry])

            row_num += 1

        # Update canvas settings
        canvas.create_window((0, 0), window=frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        frame.bind(
            "<Configure>",
            lambda event: canvas.configure(scrollregion=canvas.bbox("all")),
        )

        self.frames[name] = (canvas, scrollbar, parameters)
        self.tree.insert(parent, "end", text=name, open=False)

    def on_type_selected(self, type_variable, bool_options=None):
        selected_type = type_variable.get()

        parameters_to_hide = set(self.type_to_parameters_mapping.get(selected_type, []))

        if bool_options is not None:
            for _, (name, var) in enumerate(bool_options):
                if name == selected_type:
                    var.set(True)
                else:
                    var.set(False)

        for param, widgets in self.option_widget_map.items():
            if param in parameters_to_hide:
                if len(widgets) > 1:
                    widgets[1].grid_remove()

            else:
                for widget in widgets:
                    widget.grid()
    
    def on_button_pressed(self):
        default_directory = "../Data"
        temp_directory = "Data"       
        
        file_path = filedialog.askopenfilename(
            initialdir=default_directory,
            title="Select CSV input file",
            filetypes=[("CSV files", "*.csv")],
        )
        if file_path:  # Check if the user didn't cancel the dialog
            section_data = self.data.get("file_input", {})
            data_index = file_path.find(f'{temp_directory}')
            file_path = "../" + file_path[data_index:]
            section_data["input_file_path"] = file_path  # Add the file path to the section data
            self.data["file_input"] = section_data
        else:
            messagebox.showinfo("Operation Cancelled", "Input file selection was cancelled.")
            return
        
    def on_button_pressed_2(self):
        default_directory = "../Data"
        temp_directory = "Data"       
        
        file_path = filedialog.askopenfilename(
            initialdir=default_directory,
            title="Select CSV input file",
            filetypes=[("CSV files", "*.csv")],
        )
        if file_path:  # Check if the user didn't cancel the dialog
            section_data = self.data.get("bmanuel_input", {})
            data_index = file_path.find(f'{temp_directory}')
            file_path = "../" + file_path[data_index:]
            section_data["bmanuel_input_file"] = file_path  # Add the file path to the section data
            self.data["bmanuel_input"] = section_data
        else:
            messagebox.showinfo("Operation Cancelled", "Input file selection was cancelled.")
            return

    def on_item_click(self, event):
        # Delay the execution to allow the interface to update
        self.after(100, lambda: self.on_item_single_click(event))

    def on_item_single_click(self, event):
        item_id = self.tree.identify_row(event.y)
        if not item_id:
            return

        # Calculate the indentation level
        indentation_level = 0
        parent_id = self.tree.parent(item_id)
        while parent_id:
            indentation_level += 1
            parent_id = self.tree.parent(parent_id)

        item_bbox = self.tree.bbox(item_id)

        leftmost_area = 20 + (indentation_level * 20)
        
        if item_bbox and 0 <= event.x - item_bbox[0] <= leftmost_area:
            return

        selected_items = self.tree.selection()
        if selected_items:
            item = selected_items[0]
            children = self.tree.get_children(item)
            
            # Check if the item has children to determine if it's expandable
            if children:
                if self.tree.item(item, "open"):
                    # If the item is already open, close it
                    self.tree.item(item, open=False)
                else:
                    # If the item is closed, open it
                    self.tree.item(item, open=True)
            else:
                # Handle non-expandable items here
                for frame_name, (canvas, scrollbar, _) in self.frames.items():
                    canvas.grid_remove()
                    scrollbar.grid_remove()
                    if self.tree.item(item, "text") == frame_name:
                        canvas.grid(row=0, column=1, sticky="w", padx=10, pady=10)
                        scrollbar.grid(row=0, column=2, sticky="ns")
        else:
            print("No item selected")

    def apply_syntax_highlighting(self, text_widget):
        patterns = {
            r'[\{\}\[\]\:\,]': 'special',  # Special characters (brackets, colon, comma)
            r'\"[^\"]*\"\s*(?=\,)': 'value',  # JSON values (anything after a colon)  
            r'\"[^\"]*\"\s*(?=\:)': 'key',  # JSON keys (strings before a colon)          
        }
        
        colors = {
            'grey': '#696969',
            'brown': '#A52A2A',
            'black': '#000000',
        }

        # Configure tag styles
        text_widget.tag_configure('key', foreground=colors['black']) 
        text_widget.tag_configure('value', foreground=colors["brown"])  
        text_widget.tag_configure('special', foreground=colors["grey"])  
        

        # Apply tags to the text
        for pattern, tag in patterns.items():
            start = '1.0'
            while True:
                match = text_widget.search(pattern, start, stopindex=tk.END, regexp=True)
                if not match:
                    break
                end = f"{match}+{len(text_widget.get(match, f'{match} lineend'))}c"
                text_widget.tag_add(tag, match, end)
                start = end

    def load_defaults(self):
        default_directory = "SimParams"
        file_path = filedialog.askopenfilename(
            title="Select file", 
            initialdir=default_directory,
            filetypes=(("JSON files", "*.json"), ("All files", "*.*"))
        )
        if file_path:
            with open(file_path, "r") as file:
                data = json.load(file)
                # Set the values for each parameter
                for frame_name, (_, _, parameters) in self.frames.items():
                    section_data = data.get(frame_name.lower().replace(" ", "_"), {})
                    for param_name, var in parameters:
                        param_key = param_name.split("(")[0].strip()
                        if param_key in section_data:
                            var.set(section_data[param_key])
                        else:
                            print(f"Key '{param_key}' not found in '{frame_name}'")
        else:
            messagebox.showinfo("Operation Cancelled", "File selection was cancelled.")
            return
        
        

    def save_to_json(self):
        # Check compatibility of inputs
        for frame_name, (_, _, parameters) in self.frames.items():
            for param_name, var in parameters:
                try:
                    value_str = var.get()
                    if isinstance(var, tk.DoubleVar):
                        # Explicitly check if the value is float-compatible
                        _ = float(value_str)
                    elif isinstance(var, tk.StringVar):
                        # Check if the value is a float
                        if "." in value_str:
                            raise ValueError
                except (ValueError, tk.TclError):
                    messagebox.showerror(
                        "Input Error", f"Invalid input for '{param_name}'."
                    )
                    return
        
        # Set default directory path
        default_directory = "SimParams/Private"

        # Prompt user for the file name
        file_name = filedialog.asksaveasfilename(
            initialdir=default_directory,
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        )

        # Check if user didn't cancel the dialog
        if not file_name:
            return

        
        for frame_name, (_, _, parameters) in self.frames.items():
            section_data = {}
            for param_name, var in parameters:
                param_key = param_name.split("(")[
                    0
                ].strip()  # Using the label name before the '(' as the key.

                # Check if it's meant to be an int and convert it
                if isinstance(var, tk.StringVar) and var.get().isnumeric():
                    section_data[param_key] = int(var.get())
                else:
                    section_data[param_key] = var.get()

            self.data[frame_name.lower().replace(" ", "_")] = section_data
            
        # Convert data to JSON string
        json_str = json.dumps(self.data, indent=4)

        # Create new window for JSON preview
        preview_window = tk.Toplevel(self)
        preview_window.title("JSON Preview")
        preview_window.geometry(f"{self.DEFAULT_WIDTH}x{self.DEFAULT_HEIGHT}+{self.DEFAULT_X}+{self.DEFAULT_Y}")

        # Add Text widget with scroll bar
        text_area = scrolledtext.ScrolledText(preview_window, wrap=tk.WORD)
        text_area.insert(tk.INSERT, json_str)
        text_area.pack(expand=True, fill='both')

        # Apply syntax highlighting
        self.apply_syntax_highlighting(text_area)

        # Function to save edited JSON
        def save_edited_json():
            edited_json_str = text_area.get("1.0", tk.END)
            try:
                edited_data = json.loads(edited_json_str)
                with open(file_name, "w") as file:
                    json.dump(edited_data, file, indent=4)
                print(f"Saved to {file_name}")
                self.saved_file_name = file_name
                preview_window.destroy()
                self.destroy()
            except json.JSONDecodeError:
                messagebox.showerror("Error", "Invalid JSON format.")

        # Save button
        style = ttk.Style()
        style.configure("Centered.TButton", padding=(0, 0, 0, 0), font=("Futura", 10))
        save_button = ttk.Button(preview_window, text="Save", command=save_edited_json, width=15, style="Centered.TButton")

        # Adjust the button's placement
        save_button.pack(side='bottom', pady=(0, 2))
