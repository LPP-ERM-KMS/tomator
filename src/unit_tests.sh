#! /bin/bash

# This script should be run from the root directory of the project (tomator) with the comand "sh src/unit_tests.sh".
# It compiles the project, runs the simulation and checks the results.
# If the last message is "Files are similar within the specified tolerance." then the test is passed.
cd src && make clean && cd ..
sh Scripts/compile.sh
sh Scripts/run.sh
sh Scripts/check.sh