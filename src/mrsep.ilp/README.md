# Maximum Read Subset Explanation Problem ILP implementation with Gurobi

## Setup
Install the [Nix package manager](https://nixos.org/nix/):
```
curl https://nixos.org/nix/install | sh
source ~/.nix-profile/etc/profile.d/nix.sh
```

Obtain a [Gurobi License](http://www.gurobi.com/downloads/licenses/license-center). You can run `grbgetkey` as follows:
```
nix-shell --run "grbgetkey KEY"
```

## Build
```
./build.sh
```

## Run
```
./mrsep.ilp 0 data/strain_types_and_read_mappings.txt
```
