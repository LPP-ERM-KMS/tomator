Installation
============

Cloning & Building
------------------

First, clone the repository to a location you'd like::

    git clone https://github.com/LPP-ERM-KMS/tomator.git

Then navigate to the src folder in tomator::

    cd tomator/src

Create a build directory in this location::

    mkdir build

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

Troubleshooting
---------------

If you are running macos on one of the ARM macs, openmp is not directly supported under the latest
clang, to fix this install llvm using homebrew and export the variables::

    export CC=$(brew --prefix llvm)/bin/clang
    export CXX=$(brew --prefix llvm)/bin/clang++
