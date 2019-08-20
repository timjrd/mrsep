args@{ pkgs ? import ./pkgs.nix }:
with import pkgs {};

stdenvNoCC.mkDerivation {
  name = "simulation-coverage";
  buildInputs = [
    ((python3.withPackages (x: with x; [
      biopython      
    ])).override (_:{ ignoreCollisions = true; }))
  ];
  ART_ILLUMINA = "${import ./art.nix args}/bin/art_illumina";
  EMBOSS_MSBAR = "${emboss}/bin/msbar";
  MUSCLE       = "${muscle}/bin/muscle";
}
