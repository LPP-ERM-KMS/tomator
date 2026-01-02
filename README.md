# **Tomator: A 1D Plasma Simulation inside a Tokamak** 
![logo](logo/logocolortorus.svg)

[![CMake on ubuntu](https://github.com/LPP-ERM-KMS/tomator/actions/workflows/cmake-single-platform.yml/badge.svg?branch=master)](https://github.com/LPP-ERM-KMS/tomator/actions/workflows/cmake-single-platform.yml)

Tomator performs a 1D simulation of the plasma inside a tokamak which
has succesfully reproduced a TCV helium ECRH plasma [[1]](#1) and an ECWC
plasma in ASDEX Upgrade [[2]](#2). 

This repository contains two versions:
- **C++ version** (`src/`): The original implementation using cmake
- **T1D-lite** (`python/`): A Python/FEniCSx implementation that is faster than the C++ version but not yet fully tested

The full documentation can be found [here](https://lpp-erm-kms.github.io/tomator/) and
is the suggested reference, for those already familiar with the software, a
quick start guide is found below.

---

## T1D-lite (Python Version)

T1D-lite is a Python port of Tomator using the FEniCSx/dolfinx finite element framework. It offers faster execution than the C++ version but is still under active development and testing.

### Installation

1. Create a conda environment with the required dependencies:

```console
LPP@ERM/KMS:~$ cd tomator/python
LPP@ERM/KMS:~/tomator/python$ conda env create -f environment.yml
LPP@ERM/KMS:~/tomator/python$ conda activate t1dl-env
```

2. Set the results directory (optional):

```console
LPP@ERM/KMS:~$ export TOMATORRESULTS=~/TomatorResults
```

### Running a Simulation

Run a simulation using the JSON input file:

```console
LPP@ERM/KMS:~/tomator/python$ python -m tomator_dolfinx.solver tomator_dolfinx/examples/TCV5151X_fixPDV.json
```

Example input files are provided in `python/tomator_dolfinx/examples/`:
- `TCV5151X_fixPDV.json` - Fixed power fraction coupling
- `TCV5151X_fixneDV.json` - Fixed ne through power coupled with PID control

### Input File Format

JSON input files use a value-unit format for clarity but note that these cannot be changed (yet):

```json
"magnetic_field": {
    "Bt": {"value": 1.54, "unit": "T"},
    "Bv": {"value": 0.0035, "unit": "T"},
    "Bh": {"value": 0.0001, "unit": "T"}
}
```

### Simulation Interface 

Not yet supported, but you may edit the json input files in a text editor

### Plotting Results

Use the Bokeh-based plotter interface (it should open automatically when launching a simulation):

```console
LPP@ERM/KMS:~/tomator/python$ python -m tomator_dolfinx.gui.plot_app
```

Then open `http://localhost:5006` in your browser.

---

## C++ Version (Original) 

### Building & Installation

The build steps are the same as other cmake software:

```console
LPP@ERM/KMS:~$ git clone https://github.com/LPP-ERM-KMS/tomator.git
LPP@ERM/KMS:~$ cd tomator/src
LPP@ERM/KMS:~/tomator/src$ mkdir build
LPP@ERM/KMS:~/tomator/src$ cd build
LPP@ERM/KMS:~/tomator/src/build$ cmake ..
LPP@ERM/KMS:~/tomator/src/build$ make -jX
```
Where X is the amount of threads your cpu has (use nproc to find out or omit
the j flag if in doubt)

This will have built an executable called 'Tomator1D' which is the primary
binary. Optionally you may install the software:

```console
LPP@ERM/KMS:~/tomator/src/build$ sudo make install
```

Now export the TOMATORSOURCE and  TOMATORRESULTS environmental
variable as the absolute location of the tomator source folder and the location
of where you want your results stored, e.g add them in bashrc (change
the first directory to the one where you installed tomator and the second
to where you want the results stored, make sure you created the directory)::

```console
LPP@ERM/KMS:~$ echo "export TOMATORSOURCE=~/tomator" >> ~/.bashrc
LPP@ERM/KMS:~$ echo "export TOMATORRESULTS=~/TomatorResults" >> ~/.bashrc
```

Adding these to bashrc will make the variables persist across sessions.  If
modifications were made to the software or you wish to rebuild for other
reasons make sure to clean first before rebuilding:

```console
LPP@ERM/KMS:~/tomator/src/build$ make clean && make -jX
```

And optionally re-install.

## Setting up system parameters, running and monitoring the simulation (C++)

Two gui applications were created, SimulationInterface and PlotterInterface,
both located in the gui folder:

```console
LPP@ERM/KMS:~/tomator/gui$ python SimulationInterface.py 
``` 

```console
LPP@ERM/KMS:~/tomator/gui$ python PlotterInterface.py 
``` 

The SimulationInterface program allows a user to define a simulation and/or load a pre-defined simulation
but also run the simulation. The PlotterInterface software allows a user to monitor the simulation result
in real time, by navigating it to the output csv file (in the folder defined in TOMATORRESULTS/inputfilename).

<a id="1">[1]</a> 
Wauters, T., Buermans, J., Haelterman, R., Moiseenko, V., Ricci, D., Verhaeghe, T., … the EUROfusion MST1 team. (2020). RF plasma simulations using the TOMATOR 1D code : a case study for TCV helium ECRH plasmas. PLASMA PHYSICS AND CONTROLLED FUSION, 62(10). https://doi.org/10.1088/1361-6587/aba767

<a id="2">[2]</a> 
Wauters, T., Buermans, J., Cavalier, J., Huett, E., Ragona, R., Svoboda, J., … Eurofusion Mst Team, [missing]. (2023). Characterisation of electron cyclotron wall conditioning plasma in ASDEX Upgrade. NUCLEAR FUSION, 63(6). https://doi.org/10.1088/1741-4326/acc674
