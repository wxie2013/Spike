#!/bin/bash

mkdir Build
cd Build
cmake ../
make -j8
make install
