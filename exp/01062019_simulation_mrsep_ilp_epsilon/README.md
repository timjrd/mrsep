# MRSEP input simulation from real data

## Setup
Install the [Nix package manager](https://nixos.org/nix/):
```
curl https://nixos.org/nix/install | sh
source ~/.nix-profile/etc/profile.d/nix.sh
```

## Looking at the results
```
./notebook.sh
```
Copy/paste the given URL into your browser, then open `doc/epsilon_plots.ipynb`.

## Reproducing the experiment
Edit `src/main.py`, then:
```
./run.sh
```
