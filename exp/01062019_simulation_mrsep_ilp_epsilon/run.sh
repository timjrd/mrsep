#! /usr/bin/env nix-shell
#! nix-shell --pure -i bash
(
  cd ../19052019_simulation_coverage_statistics
  test -d results || tar -xf results.tar.xz
)

export PYTHONPATH=$PYTHONPATH:../19052019_simulation_coverage_statistics/src
python3 src/main.py
