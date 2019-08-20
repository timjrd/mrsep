#! /usr/bin/env nix-shell
#! nix-shell --pure -i bash
function unpack {
  test -d $1 || tar -xf $1.tar.xz
}
unpack results-000
unpack results-000-no-weights
unpack results-001
unpack results-002
unpack results-003
unpack results-004-null-filter
unpack results-004-null-const-filters

(
  cd ../19052019_simulation_coverage_statistics
  test -d results || tar -xf results.tar.xz
)

jupyter-notebook
