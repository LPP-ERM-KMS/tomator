# #!/bin/bash

# This script compiles the project, runs the simulation and checks the results.
# If the last message is "Files are similar within the specified tolerance." then the test is passed.
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

cd "$parent_path"

cd ../src && make clean && cd ..
sh tests/Scripts/compile.sh
sh tests/Scripts/run.sh
sh tests/Scripts/check.sh
