{ pkgs ? import ./pkgs.nix }:
with pkgs;

stdenv.mkDerivation {
  name         = "gurobi-mrsep";
  src          = ./.;
  buildInputs  = [gurobi];
  buildPhase   = "g++ main.cc -O3 -lgurobi80 -lgurobi_c++ -o gurobi-mrsep";
  installPhase = "mkdir -p $out/bin ; mv gurobi-mrsep $out/bin";
}
