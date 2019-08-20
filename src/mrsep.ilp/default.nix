{ pkgs ? import ./pkgs.nix }:
let config = (import pkgs {}).config;
in with import pkgs {config = config // {
  allowUnfree = true;
};};

stdenv.mkDerivation {
  name         = "mrsep.ilp";
  src          = ./.;
  buildInputs  = [gurobi];
  buildPhase   = "bash build.sh";
  installPhase = "mkdir -p $out/bin ; mv mrsep.ilp $out/bin";
}
