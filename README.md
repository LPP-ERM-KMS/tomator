# **Tomator: A 1D Plasma Simulation inside a Tokamak** 
![logo](logo/logocolor.svg)

[![CMake on ubuntu](https://github.com/LPP-ERM-KMS/tomator/actions/workflows/cmake-single-platform.yml/badge.svg?branch=master)](https://github.com/LPP-ERM-KMS/tomator/actions/workflows/cmake-single-platform.yml)

Tomator performs a 1D simulation of the plasma inside a tokamak which
has succesfully reproduced a TCV helium ECRH plasma [[1]](#1) and an ECWC
plasma in ASDEX Upgrade [[2]](#2). The project is
built in C++ and uses cmake to build it's binary.  The binary takes as argument
a json file which specifies the simulation. The software comes with two main
interfaces: one for setting up the json file and executing the simulation, and another for
the real-time visualization of simulation results. 

The full documentation can be found "here" (to be turned into website link) and
is the suggested reference, for those already familiar with the software, a
quick start guide is found below


## Building & **Installation** 

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

## Setting up system parameters, running and monitoring the simulation

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
