import tkinter as tk
from tkinter import ttk, filedialog
import os
import subprocess

# Assuming your original script is saved as `select_program.py`
import InterfaceFuncs.ChooseParameters as select_program


class ParamSelector(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Simulation Interface")
        self.option_add("*Font", "Futura")
        self.after(1, self.setup_geometry_and_ui)

    def setup_geometry_and_ui(self):
        DEFAULT_WIDTH = 800
        DEFAULT_HEIGHT = 500
        DEFAULT_X = self.winfo_screenwidth() // 2 - DEFAULT_WIDTH // 2
        DEFAULT_Y = self.winfo_screenheight() // 2 - DEFAULT_HEIGHT // 2

        self.geometry(f"{DEFAULT_WIDTH}x{DEFAULT_HEIGHT}+{DEFAULT_X}+{DEFAULT_Y}")
        self.configure_main_interface()

    def configure_main_interface(self):
        # Set a frame to hold the widgets
        self.main_frame = ttk.Frame(self, padding="20")
        self.main_frame.pack(padx=40, pady=40, expand=True, fill="both")

        # Label at the top
        label = ttk.Label(
            self.main_frame, text="Simulation Parameters Selection", font=("Arial", 16)
        )
        label.pack(pady=20)

        # Button to allow user to define parameters
        self.define_btn = ttk.Button(
            self.main_frame,
            text="Define Parameters",
            command=self.launch_define_params,
            width=30,
        )
        self.define_btn.pack(pady=10)

        # Button to allow user to select pre-defined parameters
        self.select_btn = ttk.Button(
            self.main_frame,
            text="Select Pre-defined Parameters",
            command=self.launch_select_params,
            width=30,
        )
        self.select_btn.pack(pady=10)

        # Style configurations
        self.style = ttk.Style()
        self.style.configure("TButton", font=("Futura", 12), padding=10)
        self.style.configure("TLabel", font=("Futura", 16))

    def run_simulation(self, filename):
        executable_path = "./../src/build/Tomator1D"
        
        if not os.path.exists(executable_path):
            print("Executable not found, please compile as instructed in the documentation")
            return
        
        self.destroy()
        print("Running simulation...")
        subprocess.run([executable_path, filename])

    def display_run_or_cancel(self, filename):
        # Create a new frame to ask user whether to run or cancel
        self.run_or_cancel_frame = ttk.Frame(self, padding="20")
        self.run_or_cancel_frame.pack(padx=40, pady=40, expand=True, fill="both")

        label = ttk.Label(
            self.run_or_cancel_frame,
            text=f"Do you want to run the simulation according to ",
            font=("Futura", 16),
        )
        label.pack(pady=0)

        label2 = ttk.Label(
            self.run_or_cancel_frame,
            text=f"{os.path.basename(filename)} parameters?",
            font=("Futura", 16),
        )
        label2.pack(pady=10)

        run_btn = ttk.Button(
            self.run_or_cancel_frame,
            text=f"Run simulation",
            command=lambda: self.run_simulation(filename),
            width=30,
        )
        run_btn.pack(pady=0)

        cancel_btn = ttk.Button(
            self.run_or_cancel_frame,
            text="Cancel",
            command=self.return_to_main,
            width=30,
        )
        cancel_btn.pack(pady=20)

    def display_json_content(self, filename):
        # Remove the previous frame
        if hasattr(self, "display_frame"):
            self.display_frame.destroy()

        # Display the run or cancel frame
        self.display_run_or_cancel(filename)

    def launch_define_params(self):
        self.withdraw()
        # Launch the SimulationParamsApp from your original script
        app = select_program.SimulationParamsApp(self)
        self.wait_window(app)

        # Get the saved filename from the app and print its content
        saved_file_name = getattr(app, "saved_file_name", None)
        if saved_file_name:
            self.deiconify()
            self.main_frame.destroy()
            self.display_json_content(saved_file_name)
        else:
            self.deiconify()

    def launch_select_params(self):
        default_directory = "../examples/InputFiles"

        # Prompt user to select a file
        file_name = filedialog.askopenfilename(
            initialdir=default_directory,
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        )

        # Check if user didn't cancel the dialog
        if not file_name:
            return

        # Clear the main interface
        self.main_frame.destroy()

        # Display the content of the selected file
        self.display_json_content(file_name)

    def return_to_main(self):
        if hasattr(self, "display_frame"):
            self.display_frame.destroy()
        if hasattr(self, "run_or_cancel_frame"):
            self.run_or_cancel_frame.destroy()
        self.configure_main_interface()


if __name__ == "__main__":
    selector = ParamSelector()
    selector.mainloop()
