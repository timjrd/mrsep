args@{ pkgs ? import ./pkgs.nix }:
with import pkgs {};

stdenvNoCC.mkDerivation {
  name = "filters";
  buildInputs = [
    ((python3.withPackages (x: with x; [
      #biopython      
      notebook
      matplotlib      
    ])).override (_:{ ignoreCollisions = true; }))
  ];
}
