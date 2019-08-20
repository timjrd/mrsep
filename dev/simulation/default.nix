args@{ pkgs ? import ./pkgs.nix }:
with pkgs;

stdenvNoCC.mkDerivation {
  name = "simulation";
  buildInputs = [
    (python.withPackages (x: with x; [
      biopython
      pysam
    ]))
  ];
  ART_ILLUMINA = "${import ./art.nix args}/bin/art_illumina";
}
