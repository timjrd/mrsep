#!/usr/bin/env bash
g++ "$@" -std=c++17 src/main.cc -O3 -lgurobi80 -lgurobi_c++ -o mrsep.ilp
