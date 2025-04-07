#! /bin/bash

source Scripts/environment.sh

make -C src clean
make -C src -j4 all
