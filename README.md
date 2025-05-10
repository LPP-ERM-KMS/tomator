# **Tomator: A 1D Plasma Simulation inside a Tokamak** 

![logo](logo/logocolor.svg)

Tomator performs a 1D simulation of the plasma inside a tokamak. The project is
built in C++ and uses cmake to build it's binary.  The binary takes as argument
a json file which specifies the simulation. The software comes with two main
interfaces: one for setting up the json file and executing the simulation, and another for
the real-time visualization of simulation results. 

The documentation can be found "here" (to be turned into website link),
for those already familiar with the software, a quick start guide is found
below

## **Table of Contents** 

- [Getting Started](#getting-started) 
- [Setting Up & Running a Simulation](#setting-up--running-a-simulation) 
- [Plotting the Simulation Output](#plotting-the-simulation-output) 
- [Additional Features](#additional-features) 
- [Adding a New Coupled Power Function](#adding-a-new-coupled-power-function) 
- [How to Set Up Remote Execution](#how-to-set-up-remote-execution)
- [Conclusion](#conclusion) 

## **Getting Started** 

### **Installation** 

First, clone the repository to a location you'd like::

```console
LPP@ERM/KMS:~$ git clone https://github.com/LPP-ERM-KMS/tomator.git
```

Then navigate to the src folder in tomator::

```console
LPP@ERM/KMS:~$ cd tomator/src

```

Create a build directory in this location::

```console
LPP@ERM/KMS:~$ mkdir build

```

navigate into this directory::

    cd build

Generate the make file::

    cmake ..

and finally build the software (replace the number after j with the number 
of threads, or if you don't know omit the flag)::

    make -j8

This will have built an executable called 'Tomator1D' which is the primary
binary. Optionally you may install the software by running::

    sudo make install

Now export the TOMATORSOURCE and  TOMATORRESULTS environmental
variable as the absolute location of the tomator source folder and the location
of where you want your results stored, e.g in bashrc add::

    export TOMATORSOURCE=/home/lpp/programs/tomator
    export TOMATORRESULTS=/home/lpp/TomatorResults

adding these to bashrc will make the variables persist across sessions. You are
now done and may move on to :doc:`Usage`.  If modifications were made to the
software or you wish to rebuild for other reasons make sure to clean first
before rebuilding::

    make clean && make -j8

And optionally re-install.

To kickstart the simulation interface, utilize the following command: 

```bash 
python SimulationInterface.py 
``` 

Within this interface, users are presented with two options: 

1. **Run a Predefined Simulation**: Opt for a preset parameter file to begin a simulation with established parameters. 

2. **Customize & Run**: Modify the parameters according to your needs. This action will spawn a new parameter file, post which you can execute the simulation with these bespoke parameters. 

Parameters are neatly categorized in the `SimParams` directory, which is bifurcated into: 

- **Public**: Houses parameter files that are thoroughly vetted, tested, and universally accessible. 

- **Private**: Designated for individual experimentation and personal parameter file storage. 

### **Plotting the Simulation Output** 

To activate the plotting interface, deploy the following command: 

```bash 
python PlotterInterface.py 
``` 

This interface empowers users to: 

- **Observe Real-time Results**: As the simulation is in progress, this platform will dynamically update with the plotted results, eliminating manual intervention. 

- **Interactive Plotting**: A mere click on the first plot allows users to pinpoint a radius. Subsequently, the third plot unravels the electron density values corresponding to the chosen radius, encompassing the entire simulation duration. If the simulation is ongoing, this plot will continually refresh. 

## **Additional Features** 

- **Concurrent Simulations**: Tomator is adept at operating and supervising multiple simulations simultaneously. 

- **Terminate Server**: For users with active simulations, a termination server is at your disposal. Hit it to designate the server port you intend to deactivate. The corresponding port of the simulation can be identified in the browser via the `localhost:` prefix, followed by the port number. 

## **Adding a New Coupled Power Function** 

To extend the simulation with a new coupled power function, follow these steps: 

### Step 1: Create the Coupling Function File 
- Create a C++ source file named `coupledB<name>.cpp` and include the `coupledpower.h` header. 
- Define your coupling functions within this file. 

### Step 2: Update `coupledpower.cpp` 
- Integrate your function by adding an if-statement to call it based on the `<name>` boolean identifier. 

### Step 3: Define the New Boolean Parameter 
- Declare a boolean parameter in `simparam.h` and set its default value to `false` in `simparam.cpp`. 

### Step 4: Modify the Python Interface 
- Add a new entry for your function's boolean in `ChooseParameters.py` - `Type` group, and update `self.type_to_parameters_mapping` with the parameters that aren't important for your function. 

### Step 5: Configure the JSON File Using the Simulation Interface 
- Use `SimulationInterface.py` to generate a JSON configuration file that includes your new parameter. 

### Step 6: Update the Tomator Code 
- In `extractor.cpp`, add an if-condition to the `extract_type` function for your parameter. 

After completing these steps, run the simulation with the new coupled power function integrated. 

## **How to Set Up Remote Execution**

Follow these steps if you choose to run your simulation on a remote cluster:

### Step 1: 
- Create or choose the simulation parameters within the local Tomator interface. A file with these parameters will be generated.

### Step 2: 
- When prompted by the interface to run the simulation, select 'Cancel' if you wish to execute it remotely instead of locally.

### Step 3: 
- Locate the provided Bash script that is intended for remote execution. This script will need to be edited to include your specific cluster configurations.

### Step 4: 
- Open the Bash script and manually replace the placeholder with the actual name of the parameter file you created. For example, if your parameter file is named 'experiment1', change the line in the script to ./Tomator1D experiment1.

### Step 5: 
- Additionally, configure any other necessary settings in the script, such as the destination cluster's address, user credentials, and file paths.

### Step 6: 
- Execute the Bash script. It will establish a connection to the remote cluster, transfer the necessary files, and initiate the simulation using the command ./Tomator1D <nameSimParams>.

## **Conclusion** 

Tomator is a powerful tool that caters to both novices and seasoned users, enabling them to simulate and elucidate the behavior of plasma inside a tokamak in a 1D space. Irrespective of using stock parameters or personalized setups, the software delivers real-time, interactive insights, promising an immersive user journey. 

> **Note**: Ensure the installation of all prerequisite dependencies and adherence to the suggested Python version for an optimal experience. 

## **Requirements** 

### **Python Packages** 

- os 
- time 
- pandas 
- bokeh 
- tkinter 
- subprocess 
- math 
- socket 

Ensure you have these dependencies installed before running the simulation. You can typically install the Python packages using pip: 

```bash 
pip install pandas bokeh tkinter 
``` 

> **Note**: Many of the Python modules like os, time, math, and socket are part of the standard 
