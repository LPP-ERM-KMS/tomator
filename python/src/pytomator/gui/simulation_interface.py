"""
Simulation Interface for Tomator1D Python.

A tkinter-based GUI for editing JSON simulation parameters with proper units.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import json
import os
from pathlib import Path


# ============================================================================
# PARAMETER DEFINITIONS WITH UNITS AND DESCRIPTIONS
# ============================================================================
# Format: "json_key": ("Display Name", "Unit", "Description")

PARAMETER_DEFINITIONS = {
    # Magnetic field
    "Bt": ("Toroidal Field", "T", "Toroidal magnetic field at R0"),
    "Bv": ("Vertical Field", "T", "Vertical magnetic field component"),
    "Bh": ("Horizontal Field", "T", "Horizontal magnetic field component"),
    
    # Geometry
    "R": ("Major Radius", "cm", "Major radius of torus"),
    "a": ("Minor Radius", "cm", "Minor radius (half-width in radial direction)"),
    "b": ("Vertical Extent", "cm", "Vertical extent (half-height)"),
    "lHFS": ("HFS Limiter", "cm", "High-field side limiter position (R - a + offset)"),
    "lLFS": ("LFS Limiter", "cm", "Low-field side limiter position (R + a + offset)"),
    "nlimiters": ("Number of Limiters", "", "Number of poloidal limiters"),
    "Vpl": ("Plasma Volume", "cm³", "Total plasma volume"),
    
    # Neutral pressure
    "pHe": ("Helium Pressure", "Pa", "Neutral helium partial pressure"),
    "pH2": ("Hydrogen Pressure", "Pa", "Neutral H₂ partial pressure"),
    
    # RF power
    "Prf": ("RF Power", "kW", "Total injected RF power"),
    "bselfcol": ("Self-collision", "", "Include RF-induced self-collisions"),
    "freq": ("Frequency", "MHz", "RF frequency"),
    "dtpramp": ("Power Ramp Time", "s", "Duration of power ramp-up"),
    
    # Power coupling type
    "bgray": ("GRAY Code", "", "Use GRAY ray-tracing for EC"),
    "bram": ("RAM Code", "", "Use RAM full-wave code for EC"),
    "bfixpowerfrac": ("Fixed Power Fraction", "", "Fixed absorbed power = Prf × pecabs0"),
    "bnefix": ("Fixed ne Coupling", "", "PID controller maintains ne at control point"),
    "bTOMAS": ("TOMAS Mode", "", "Multi-antenna TOMAS power deposition"),
    "bproptone": ("Proportional to ne", "", "Power deposition proportional to ne"),
    "bnopower": ("No Power", "", "Disable RF power coupling"),
    "bkipt": ("KIPT Mode", "", "KIPT-style power coupling"),
    "blhr": ("LH Resonance", "", "Lower Hybrid Resonance heating"),
    "bmanuel": ("Manual Power", "", "Read power profile from file"),
    "bICWC": ("ICWC Mode", "", "Ion Cyclotron Wall Conditioning mode"),
    
    # EC general
    "Rdep": ("Deposition Radius", "cm", "Radial location of power deposition center"),
    "pecabs0": ("Absorbed Fraction", "", "Initial/fixed absorbed power fraction"),
    "widthech": ("Deposition Width", "cm", "Gaussian width of EC power deposition"),
    "echbackground": ("EC Background", "", "Background absorption level"),
    "harmonic": ("EC Harmonic", "", "Cyclotron harmonic number"),
    "muw": ("μw Parameter", "", "Wave-plasma coupling parameter"),
    
    # necfix (PID controller)
    "ic": ("Control Index", "", "Mesh index of density control point"),
    "necfix": ("Target ne", "cm⁻³", "Target electron density for PID"),
    "P_KI_ini": ("Initial Integral", "", "Initial integral term value"),
    "P_KP": ("PID Kp", "", "Proportional gain"),
    "P_KI": ("PID Ki", "1/s", "Integral gain"),
    "P_KD": ("PID Kd", "s", "Derivative gain"),
    
    # TOMAS
    "Rdep1": ("Antenna 1 Position", "cm", "Radial position of antenna 1"),
    "PRdep1": ("Antenna 1 Power", "", "Power fraction to antenna 1"),
    "Rdep2": ("Antenna 2 Position", "cm", "Radial position of antenna 2"),
    "PRdep2": ("Antenna 2 Power", "", "Power fraction to antenna 2"),
    "Rdep3": ("Antenna 3 Position", "cm", "Radial position of antenna 3"),
    "PRdep3": ("Antenna 3 Power", "", "Power fraction to antenna 3"),
    "Rdep4": ("Antenna 4 Position", "cm", "Radial position of antenna 4"),
    "PRdep4": ("Antenna 4 Power", "", "Power fraction to antenna 4"),
    
    # IC general
    "Rant": ("Antenna Radius", "cm", "IC antenna radial position"),
    "bantlr": ("Antenna LR", "", "Antenna left-right configuration"),
    "avlr": ("LR Asymmetry", "", "Left-right power asymmetry factor"),
    "alphalaw": ("Alpha Law", "", "Power deposition decay parameter"),
    "HtoHD": ("H/(H+D) Ratio", "", "Hydrogen to hydrogen+deuterium ratio"),
    
    # LHR
    "fracpne": ("Power-ne Fraction", "", "Fraction of power proportional to ne"),
    "fraclhr": ("LHR Fraction", "", "Fraction deposited at LH resonance"),
    "widthlhr": ("LHR Width", "cm", "Width of LHR deposition layer"),
    "lhrbackground": ("LHR Background", "", "Background absorption at LHR"),
    
    # Other
    "fixedTe": ("Fixed Te", "", "Keep electron temperature fixed"),
    
    # Diffusion
    "bDfix": ("Fixed D", "", "Use constant diffusion coefficient"),
    "Dfix": ("D Fixed Value", "cm²/s", "Fixed diffusion coefficient value"),
    "bDbohm": ("Bohm D", "", "Use Bohm diffusion: D = Dfact × T/B"),
    "bDscaling": ("Scaling D", "", "Use temperature-dependent scaling"),
    "Dfact": ("D Factor", "", "Diffusion scaling factor (Dfsave)"),
    
    # Advection
    "bVfix": ("Fixed V", "", "Use constant advection velocity"),
    "Vfix": ("V Fixed Value", "cm/s", "Fixed advection velocity value"),
    "bVscaling": ("Scaling V", "", "Use pressure-driven pinch: V = -Vfact×D×∇p/p"),
    "veq": ("V Equation", "", "Advection model equation number"),
    "Vfact": ("V Factor", "", "Advection scaling factor"),
    
    # Tune D and V (PID)
    "btunedv": ("Tune D&V", "", "Enable PID tuning of D and V"),
    "btunevleft": ("Tune V Left", "", "Tune V on left boundary"),
    "il": ("Left Index", "", "Left control point mesh index"),
    "nelfix": ("Left ne Target", "cm⁻³", "Target ne at left control point"),
    "ir": ("Right Index", "", "Right control point mesh index"),
    "nerfix": ("Right ne Target", "cm⁻³", "Target ne at right control point"),
    "Vini": ("Initial V", "cm/s", "Initial advection velocity for tuning"),
    "tauV": ("V Time Constant", "s", "PID time constant for V tuning"),
    "Dini": ("Initial D", "cm²/s", "Initial diffusion for tuning"),
    "tauD": ("D Time Constant", "s", "PID time constant for D tuning"),
    
    # Physics to include
    "bH": ("H Atom Physics", "", "Include atomic hydrogen reactions"),
    "bH2": ("H₂ Physics", "", "Include molecular hydrogen reactions"),
    "bHe": ("He Physics", "", "Include helium reactions"),
    "bion": ("Ionization", "", "Include ionization reactions"),
    "bcx": ("Charge Exchange", "", "Include charge exchange reactions"),
    "belas": ("Elastic Collisions", "", "Include elastic collisions"),
    "bcoulomb": ("Coulomb Collisions", "", "Include Coulomb collisions"),
    "bimpur": ("Impurities", "", "Include impurity species"),
    "btranspions": ("Ion Transport", "", "Solve ion transport equations"),
    "btranspneut": ("Neutral Transport", "", "Solve neutral transport equations"),
    "bedge": ("Edge Effects", "", "Include edge/SOL effects"),
    "bpol": ("Poloidal Losses", "", "Include poloidal magnetic field losses"),
    "vdrift": ("Drift Velocity", "", "Include drift velocity effects"),
    "bcoll": ("Collision Freq", "", "Compute collision frequencies"),
    
    # Initial conditions
    "rmaxini": ("Initial Peak R", "cm", "Radial position of initial density peak"),
    "widthini": ("Initial Width", "cm", "Gaussian width of initial profile"),
    "nebackgroundl": ("Left Background", "cm⁻³", "Background density on HFS"),
    "nebackgroundr": ("Right Background", "cm⁻³", "Background density on LFS"),
    "Ta0": ("Ambient T", "eV", "Ambient/wall temperature"),
    "Te0": ("Initial Te", "eV", "Initial electron temperature"),
    "nevac": ("Vacuum Density", "cm⁻³", "Minimum density floor"),
    "nH0": ("Initial nH", "cm⁻³", "Initial atomic hydrogen density"),
    "nHi0": ("Initial nH⁺", "cm⁻³", "Initial H⁺ ion density"),
    "nH2i0": ("Initial nH₂⁺", "cm⁻³", "Initial H₂⁺ ion density"),
    "nH3i0": ("Initial nH₃⁺", "cm⁻³", "Initial H₃⁺ ion density"),
    "nHeII0": ("Initial nHe⁺", "cm⁻³", "Initial He⁺ ion density"),
    "nHeIII0": ("Initial nHe²⁺", "cm⁻³", "Initial He²⁺ ion density"),
    "nCII0": ("Initial nC⁺", "cm⁻³", "Initial C⁺ impurity density"),
    "nCIII0": ("Initial nC²⁺", "cm⁻³", "Initial C²⁺ impurity density"),
    "nCIV0": ("Initial nC³⁺", "cm⁻³", "Initial C³⁺ impurity density"),
    "nCV0": ("Initial nC⁴⁺", "cm⁻³", "Initial C⁴⁺ impurity density"),
    
    # Edge conditions
    "RH": ("Particle Reflection", "", "Particle reflection coefficient at wall"),
    "REH": ("Energy Reflection", "", "Energy reflection coefficient at wall"),
    "gEd": ("γ_Ed Factor", "", "Energy diffusion enhancement factor"),
    "gEv": ("γ_Ev Factor", "", "Energy advection enhancement factor"),
    "gEdn": ("γ_Edn Factor", "", "Neutral energy transport factor"),
    "gEe": ("γ_Ee Factor", "", "Electron energy transport factor"),
    
    # Simulation grid
    "nmeshp": ("Mesh Points", "", "Number of radial mesh points"),
    
    # Input file
    "bfinput": ("Use Input File", "", "Read initial conditions from file"),
    
    # Time step
    "t0": ("Start Time", "s", "Simulation start time"),
    "tmainend": ("End Time", "s", "Simulation end time"),
    "accur": ("Accuracy", "", "Target relative accuracy per time step"),
    "dtmax": ("Max dt", "s", "Maximum allowed time step"),
    "dtmin": ("Min dt", "s", "Minimum allowed time step"),
    "dtinit": ("Initial dt", "s", "Initial time step size"),
    
    # RF coupling time step
    "bupdateRFstep": ("Update RF Step", "", "Update RF coupling each step"),
    "updateRF": ("RF Update Interval", "", "Steps between RF updates"),
    "dtRF": ("RF dt", "s", "Time step for RF coupling"),
    "dtRFvar": ("Variable RF dt", "", "Allow variable RF time step"),
    "dtRFmax": ("Max RF dt", "s", "Maximum RF coupling time step"),
    "dtRFmin": ("Min RF dt", "s", "Minimum RF coupling time step"),
    "dtRFconv": ("RF Convergence", "", "Check RF convergence"),
    
    # Advanced time step
    "dtsmooth": ("Smooth dt", "", "Smooth time step changes"),
    "shokparam": ("Shock Parameter", "", "Shock detection sensitivity"),
    "maxtstepincrement": ("Max dt Increase", "", "Maximum dt increase factor"),
    "btendvar": ("Variable End", "", "Allow variable end time"),
    "convcrit": ("Convergence Criterion", "", "Convergence threshold"),
    "convsavetime": ("Conv Save Time", "s", "Time between convergence saves"),
    "baccurvar": ("Variable Accuracy", "", "Allow variable accuracy"),
    "accurcrit": ("Accuracy Criterion", "", "Accuracy threshold"),
    "minaccur": ("Min Accuracy", "", "Minimum allowed accuracy"),
    
    # Output
    "Nlog": ("Log Interval", "", "Steps between log outputs"),
    "cc": ("Central Cell", "", "Mesh index for diagnostic output"),
    "Nloopsave": ("Save Interval", "", "Steps between file saves"),
    "bOutdt": ("Output dt", "", "Include dt in output"),
    "dtsave": ("Save dt", "s", "Time interval for saving results"),
    
    # Solver
    "solvertolerance": ("Solver Tolerance", "", "Linear solver tolerance"),
}

# Section groupings for the tree view
SECTION_GROUPS = {
    "Main Discharge Parameters": [
        ("magnetic_field", "Magnetic Field"),
        ("toroidal_machine_geometry", "Geometry"),
        ("neutral_pressure", "Neutral Pressure"),
    ],
    "Coupled (RF) Power": [
        ("rf_power", "RF Power"),
        ("type", "Power Coupling Type"),
        ("general_ec", "EC Parameters"),
        ("necfix", "Density Control (PID)"),
        ("tomas", "TOMAS Antennas"),
        ("general_ic", "IC Parameters"),
        ("lhr", "LH Resonance"),
        ("other", "Other Power Options"),
    ],
    "Transport": [
        ("diffusion", "Diffusion"),
        ("advection", "Advection"),
        ("tune_d_and_v", "D&V Tuning (PID)"),
    ],
    "Physics & Conditions": [
        ("physics_to_include", "Physics to Include"),
        ("initial_conditions", "Initial Conditions"),
        ("edge_conditions", "Edge Conditions"),
    ],
    "Simulation Control": [
        ("simulation_grid", "Simulation Grid"),
        ("input_file", "Input File"),
        ("time_step", "Time Step"),
        ("time_step_for_rf_coupling", "RF Coupling Time Step"),
        ("output_parameters", "Output"),
        ("solver_parameters", "Solver"),
    ],
}


class SimulationInterface(tk.Tk):
    """Main simulation parameter interface."""
    
    def __init__(self):
        super().__init__()
        
        self.title("Tomator1D - Simulation Interface")
        self.option_add("*Font", ("Helvetica", 10))
        
        # Data storage
        self.data = {}
        self.current_file = None
        self.widgets = {}  # Store parameter widgets for updating
        
        # Window setup
        self.after(1, self._setup_ui)
    
    def _setup_ui(self):
        """Set up the main UI layout."""
        # Window size and position
        width, height = 1000, 700
        x = self.winfo_screenwidth() // 2 - width // 2
        y = self.winfo_screenheight() // 2 - height // 2
        self.geometry(f"{width}x{height}+{x}+{y}")
        
        # Configure style
        style = ttk.Style()
        style.configure("TLabel", font=("Helvetica", 10))
        style.configure("TButton", font=("Helvetica", 10), padding=5)
        style.configure("Header.TLabel", font=("Helvetica", 12, "bold"))
        style.configure("Unit.TLabel", font=("Helvetica", 9), foreground="gray")
        style.configure("Treeview", font=("Helvetica", 10))
        style.configure("Treeview.Heading", font=("Helvetica", 10, "bold"))
        
        # Main paned window
        self.paned = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        self.paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Left panel: Tree navigation
        left_frame = ttk.Frame(self.paned)
        self.paned.add(left_frame, weight=1)
        
        # Tree view
        self.tree = ttk.Treeview(left_frame, show="tree")
        self.tree.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)
        
        tree_scroll = ttk.Scrollbar(left_frame, orient=tk.VERTICAL, command=self.tree.yview)
        tree_scroll.pack(fill=tk.Y, side=tk.RIGHT)
        self.tree.configure(yscrollcommand=tree_scroll.set)
        
        # Build tree structure
        self._build_tree()
        
        # Right panel: Parameter editor
        right_frame = ttk.Frame(self.paned)
        self.paned.add(right_frame, weight=3)
        
        # Header with description
        self.header_frame = ttk.Frame(right_frame)
        self.header_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.section_label = ttk.Label(self.header_frame, text="Select a section", 
                                       style="Header.TLabel")
        self.section_label.pack(anchor=tk.W)
        
        self.description_label = ttk.Label(self.header_frame, text="", wraplength=600)
        self.description_label.pack(anchor=tk.W, pady=(5, 0))
        
        # Scrollable parameter frame
        self.param_canvas = tk.Canvas(right_frame)
        self.param_scrollbar = ttk.Scrollbar(right_frame, orient=tk.VERTICAL, 
                                             command=self.param_canvas.yview)
        self.param_frame = ttk.Frame(self.param_canvas)
        
        self.param_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10)
        self.param_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.param_canvas.configure(yscrollcommand=self.param_scrollbar.set)
        self.param_window = self.param_canvas.create_window((0, 0), window=self.param_frame, 
                                                            anchor=tk.NW)
        
        self.param_frame.bind("<Configure>", self._on_frame_configure)
        self.param_canvas.bind("<Configure>", self._on_canvas_configure)
        
        # Bind mouse wheel
        self.param_canvas.bind_all("<MouseWheel>", self._on_mousewheel)
        self.param_canvas.bind_all("<Button-4>", self._on_mousewheel)
        self.param_canvas.bind_all("<Button-5>", self._on_mousewheel)
        
        # Bottom buttons
        button_frame = ttk.Frame(self)
        button_frame.pack(fill=tk.X, padx=10, pady=10)
        
        ttk.Button(button_frame, text="Load JSON", command=self._load_json).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Save JSON", command=self._save_json).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Run Simulation", command=self._run_simulation).pack(side=tk.LEFT, padx=5)
        
        # File label
        self.file_label = ttk.Label(button_frame, text="No file loaded")
        self.file_label.pack(side=tk.RIGHT, padx=10)
        
        # Tree selection binding
        self.tree.bind("<<TreeviewSelect>>", self._on_tree_select)
    
    def _build_tree(self):
        """Build the tree view structure."""
        for group_name, sections in SECTION_GROUPS.items():
            group_id = self.tree.insert("", tk.END, text=group_name, open=False)
            for section_key, section_name in sections:
                self.tree.insert(group_id, tk.END, text=section_name, 
                               values=(section_key,), tags=(section_key,))
    
    def _on_frame_configure(self, event):
        """Update scroll region when frame size changes."""
        self.param_canvas.configure(scrollregion=self.param_canvas.bbox("all"))
    
    def _on_canvas_configure(self, event):
        """Update frame width when canvas size changes."""
        self.param_canvas.itemconfig(self.param_window, width=event.width - 20)
    
    def _on_mousewheel(self, event):
        """Handle mouse wheel scrolling."""
        if event.num == 4 or event.delta > 0:
            self.param_canvas.yview_scroll(-1, "units")
        elif event.num == 5 or event.delta < 0:
            self.param_canvas.yview_scroll(1, "units")
    
    def _on_tree_select(self, event):
        """Handle tree selection - show parameters for selected section."""
        selection = self.tree.selection()
        if not selection:
            return
        
        item = selection[0]
        values = self.tree.item(item, "values")
        
        if values:
            section_key = values[0]
            section_name = self.tree.item(item, "text")
            self._show_section(section_key, section_name)
    
    def _show_section(self, section_key: str, section_name: str):
        """Display parameters for the selected section."""
        # Update header
        self.section_label.config(text=section_name)
        
        # Clear existing widgets
        for widget in self.param_frame.winfo_children():
            widget.destroy()
        self.widgets = {}
        
        # Get section data
        section_data = self.data.get(section_key, {})
        
        if not section_data:
            ttk.Label(self.param_frame, text="No parameters in this section or no file loaded.").grid(
                row=0, column=0, pady=20)
            return
        
        # Create header row
        ttk.Label(self.param_frame, text="Parameter", font=("Helvetica", 10, "bold")).grid(
            row=0, column=0, sticky=tk.W, padx=5, pady=5)
        ttk.Label(self.param_frame, text="Value", font=("Helvetica", 10, "bold")).grid(
            row=0, column=1, sticky=tk.W, padx=5, pady=5)
        ttk.Label(self.param_frame, text="Unit", font=("Helvetica", 10, "bold")).grid(
            row=0, column=2, sticky=tk.W, padx=5, pady=5)
        ttk.Label(self.param_frame, text="Description", font=("Helvetica", 10, "bold")).grid(
            row=0, column=3, sticky=tk.W, padx=5, pady=5)
        
        ttk.Separator(self.param_frame, orient=tk.HORIZONTAL).grid(
            row=1, column=0, columnspan=4, sticky=tk.EW, pady=5)
        
        # Add parameters
        row = 2
        for key, value in section_data.items():
            # Get parameter info
            param_info = PARAMETER_DEFINITIONS.get(key, (key, "", ""))
            display_name, unit, description = param_info
            
            # Parameter name
            ttk.Label(self.param_frame, text=display_name).grid(
                row=row, column=0, sticky=tk.W, padx=5, pady=3)
            
            # Value widget
            if isinstance(value, bool):
                var = tk.BooleanVar(value=value)
                widget = ttk.Checkbutton(self.param_frame, variable=var)
            else:
                var = tk.StringVar(value=self._format_value(value))
                widget = ttk.Entry(self.param_frame, textvariable=var, width=15)
            
            widget.grid(row=row, column=1, sticky=tk.W, padx=5, pady=3)
            self.widgets[key] = (var, type(value))
            
            # Unit
            ttk.Label(self.param_frame, text=unit, style="Unit.TLabel").grid(
                row=row, column=2, sticky=tk.W, padx=5, pady=3)
            
            # Description
            ttk.Label(self.param_frame, text=description, wraplength=300).grid(
                row=row, column=3, sticky=tk.W, padx=5, pady=3)
            
            row += 1
        
        # Configure column weights
        self.param_frame.columnconfigure(3, weight=1)
    
    def _format_value(self, value) -> str:
        """Format a value for display."""
        if isinstance(value, float):
            if abs(value) < 0.001 or abs(value) >= 10000:
                return f"{value:.4e}"
            else:
                return f"{value:.6g}"
        return str(value)
    
    def _parse_value(self, string_val: str, original_type: type):
        """Parse a string value back to its original type."""
        if original_type == bool:
            return string_val.lower() in ('true', '1', 'yes')
        elif original_type == int:
            return int(float(string_val))
        elif original_type == float:
            return float(string_val)
        return string_val
    
    def _update_data_from_widgets(self):
        """Update self.data from current widget values."""
        # Get current section
        selection = self.tree.selection()
        if not selection:
            return
        
        item = selection[0]
        values = self.tree.item(item, "values")
        if not values:
            return
        
        section_key = values[0]
        
        # Update section data
        if section_key not in self.data:
            self.data[section_key] = {}
        
        for key, (var, orig_type) in self.widgets.items():
            try:
                if orig_type == bool:
                    self.data[section_key][key] = var.get()
                else:
                    self.data[section_key][key] = self._parse_value(var.get(), orig_type)
            except (ValueError, tk.TclError) as e:
                messagebox.showerror("Error", f"Invalid value for {key}: {e}")
                return False
        return True
    
    def _load_json(self):
        """Load a JSON file."""
        # Default to examples directory
        default_dir = Path(__file__).parent.parent / "examples"
        if not default_dir.exists():
            default_dir = Path.cwd()
        
        filepath = filedialog.askopenfilename(
            initialdir=default_dir,
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            title="Load Simulation Parameters"
        )
        
        if not filepath:
            return
        
        try:
            with open(filepath, 'r') as f:
                self.data = json.load(f)
            
            self.current_file = filepath
            self.file_label.config(text=f"File: {os.path.basename(filepath)}")
            
            # Show description if present
            if "description" in self.data:
                self.description_label.config(text=self.data["description"])
            else:
                self.description_label.config(text="")
            
            # Refresh current view
            self._on_tree_select(None)
            
            messagebox.showinfo("Success", f"Loaded: {os.path.basename(filepath)}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {e}")
    
    def _save_json(self):
        """Save current parameters to JSON."""
        # Update data from widgets first
        if not self._update_data_from_widgets():
            return
        
        # Default filename
        default_dir = Path(__file__).parent.parent / "examples"
        if not default_dir.exists():
            default_dir = Path.cwd()
        
        initial_file = os.path.basename(self.current_file) if self.current_file else "simulation.json"
        
        filepath = filedialog.asksaveasfilename(
            initialdir=default_dir,
            defaultextension=".json",
            initialfile=initial_file,
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            title="Save Simulation Parameters"
        )
        
        if not filepath:
            return
        
        try:
            # Show preview window
            self._show_json_preview(filepath)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save file: {e}")
    
    def _show_json_preview(self, filepath: str):
        """Show JSON preview before saving."""
        preview = tk.Toplevel(self)
        preview.title("JSON Preview")
        preview.geometry("700x500")
        
        # Text area
        text_area = scrolledtext.ScrolledText(preview, wrap=tk.WORD, font=("Courier", 10))
        text_area.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Insert JSON
        json_str = json.dumps(self.data, indent=4)
        text_area.insert(tk.END, json_str)
        
        # Buttons
        btn_frame = ttk.Frame(preview)
        btn_frame.pack(fill=tk.X, padx=10, pady=5)
        
        def save_and_close():
            try:
                # Parse edited JSON
                edited_json = text_area.get("1.0", tk.END)
                self.data = json.loads(edited_json)
                
                # Save to file
                with open(filepath, 'w') as f:
                    json.dump(self.data, f, indent=4)
                
                self.current_file = filepath
                self.file_label.config(text=f"File: {os.path.basename(filepath)}")
                
                preview.destroy()
                messagebox.showinfo("Success", f"Saved: {os.path.basename(filepath)}")
                
            except json.JSONDecodeError as e:
                messagebox.showerror("Error", f"Invalid JSON: {e}")
        
        ttk.Button(btn_frame, text="Save", command=save_and_close).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Cancel", command=preview.destroy).pack(side=tk.LEFT, padx=5)
    
    def _run_simulation(self):
        """Run the simulation with current parameters."""
        if not self.current_file and not self.data:
            messagebox.showwarning("Warning", "Please load or create parameters first.")
            return
        
        # Update data from widgets
        if not self._update_data_from_widgets():
            return
        
        # Ask for confirmation
        result = messagebox.askyesno(
            "Run Simulation",
            f"Run simulation with parameters from:\n{self.current_file or 'unsaved parameters'}?"
        )
        
        if not result:
            return
        
        # Save to temp file if needed
        if not self.current_file:
            import tempfile
            fd, temp_path = tempfile.mkstemp(suffix=".json")
            with os.fdopen(fd, 'w') as f:
                json.dump(self.data, f, indent=4)
            run_file = temp_path
        else:
            # Save current changes
            with open(self.current_file, 'w') as f:
                json.dump(self.data, f, indent=4)
            run_file = self.current_file
        
        # Run simulation
        try:
            from ..solver import run_simulation
            from .. import TomatorResults
            
            output_dir = Path(TomatorResults.__file__).parent / "output"
            
            self.withdraw()  # Hide main window
            
            run_simulation(run_file, output_dir=str(output_dir), show_plotter=True)
            
        except ImportError as e:
            messagebox.showerror("Error", f"Could not import solver: {e}")
        except Exception as e:
            messagebox.showerror("Error", f"Simulation failed: {e}")
        finally:
            self.deiconify()  # Show main window again


def main():
    """Launch the simulation interface."""
    app = SimulationInterface()
    app.mainloop()


if __name__ == "__main__":
    main()
