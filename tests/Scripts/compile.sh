#! /bin/bash

make -C src clean
make -C src -j20 all
