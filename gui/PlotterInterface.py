import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
from bokeh.plotting import figure, show, gridplot, output_file
from bokeh.models import HoverTool
import subprocess
import pandas as pd
import math
import socket


def printToAux(file_path):
    with open("InterfaceFuncs/aux.txt", "w") as f:
        f.write(file_path)


def find_free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("", 0))
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return s.getsockname()[1]


class DataInterface(tk.Tk):
    def __init__(self):
        super().__init__()
        self.processes = []

        self.title("Data Plotting Interface")
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
        # self.process_frame = ttk.Frame(self.main_frame)
        # self.process_frame.pack(pady=20)

        # Label at the top
        label = ttk.Label(self.main_frame, text="Data Plotting", font=("Futura", 16))
        label.pack(pady=20)

        # Button to allow user to run simulation
        self.follow_simulation_btn = ttk.Button(
            self.main_frame,
            text="Plot Simulation",
            command=self.follow_simulation,
            width=30,
        )
        self.follow_simulation_btn.pack(pady=10)

        self.terminate_bokeh_btn = ttk.Button(
            self.main_frame,
            text="Terminate Server",
            command=self.terminate_bokeh_server,
            width=30,
        )
        self.terminate_bokeh_btn.pack(pady=10)
        self.terminate_bokeh_btn["state"] = tk.DISABLED

        # Listbox to display active Bokeh servers
        self.bokeh_server_listbox = tk.Listbox(self.main_frame)
        self.bokeh_scrollbar = ttk.Scrollbar(
            self.main_frame, orient="vertical", command=self.bokeh_server_listbox.yview
        )
        self.bokeh_server_listbox.config(yscrollcommand=self.bokeh_scrollbar.set)

        # Style configurations
        self.style = ttk.Style()
        self.style.configure("TButton", font=("Futura", 12), padding=10)
        self.style.configure("TLabel", font=("Futura", 16))

    def center_window(self, win):
        """Centers a tkinter window"""
        win.update_idletasks()
        width = win.winfo_width()
        height = win.winfo_height()
        x = (win.winfo_screenwidth() // 2) - (width // 2)
        y = (win.winfo_screenheight() // 2) - (height // 2)
        win.geometry("{}x{}+{}+{}".format(width, height, x, y))

    def follow_simulation(self):
        # This will open a file dialog window to select a CSV file
        default_directory = "../Data"
        file_path = filedialog.askopenfilename(
            initialdir=default_directory,
            title="Select a CSV File",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if file_path:
            printToAux(file_path)
            port = find_free_port()
            process = subprocess.Popen(
                [
                    "bokeh",
                    "serve",
                    "--show",
                    "--port",
                    str(port),
                    "InterfaceFuncs/plot_app.py",
                ]
            )
            print("Process ID: ", process.pid, "on port", port)
            self.processes.append((port, process))
            self.terminate_bokeh_btn["state"] = tk.NORMAL

    def terminate_bokeh_server(self):
        self.terminate_window = tk.Toplevel(self)
        self.terminate_window.title("Terminate Bokeh Server")

        ttk.Label(self.terminate_window, text="Active Bokeh Servers:").pack(pady=10)

        self.ports_to_terminate = tk.StringVar(
            value=[proc[0] for proc in self.processes]
        )
        self.terminate_listbox = tk.Listbox(
            self.terminate_window,
            listvariable=self.ports_to_terminate,
            selectmode=tk.MULTIPLE,
        )
        self.terminate_listbox.pack(pady=10, padx=10)

        ttk.Button(
            self.terminate_window,
            text="Terminate Selected",
            command=self.kill_selected_server,
        ).pack(pady=10)

        self.center_window(self.terminate_window)

    def kill_selected_server(self):
        selected_ports = [
            int(self.terminate_listbox.get(i))
            for i in self.terminate_listbox.curselection()
        ]
        processes_to_remove = []
        for port in selected_ports:
            for i, (proc_port, proc) in enumerate(self.processes):
                if port == proc_port:
                    proc.terminate()
                    processes_to_remove.append(i)

        # Removing terminated processes from the list
        for index in sorted(processes_to_remove, reverse=True):
            del self.processes[index]

        # Refresh the listbox
        self.ports_to_terminate.set([proc[0] for proc in self.processes])

        # If no processes are running, disable the "Terminate Server" button
        if not self.processes:
            self.terminate_bokeh_btn["state"] = "disabled"

        # Destroy the termination window to go back to the main interface
        self.terminate_window.destroy()

if __name__ == "__main__":
    app = DataInterface()
    app.mainloop()
