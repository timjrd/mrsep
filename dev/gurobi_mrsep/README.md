# Maximum Read Subset Explanation Problem ILP implementation with Gurobi

## Setup
Install the [Nix package manager](https://nixos.org/nix/):
```
curl https://nixos.org/nix/install | sh
source ~/.nix-profile/etc/profile.d/nix.sh
```

Obtain a [Gurobi License](http://www.gurobi.com/downloads/licenses/license-center). You can run `grbgetkey` as follows:
```
NIXPKGS_ALLOW_UNFREE=1 nix-shell --run "grbgetkey KEY"
```

## Build
```
NIXPKGS_ALLOW_UNFREE=1 nix-build
```

## Run
```
./result/bin/gurobi-mrsep 0 strain_types_and_read_mappings.txt
```
