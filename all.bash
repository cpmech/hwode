#!/bin/bash

SOURCE=`pwd`
cd /tmp
mkdir -p build-hwode
cd build-hwode/
cmake -S $SOURCE
make

