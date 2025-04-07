#! /bin/bash

# Set up environment for compilation
source /usr/share/Modules/init/sh
module use /work/imas/etc/modules/all # What was inside the module file?? 

module load GCCcore/10.2.0
module load Python/3.8.6-GCCcore-10.2.0
module load SciPy-bundle/2020.11-foss-2020b
module load bokeh/2.2.3-foss-2020b
module load Tkinter/3.8.6-GCCcore-10.2.0