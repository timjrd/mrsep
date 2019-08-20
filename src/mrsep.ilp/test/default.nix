args@{ pkgs ? import ./pkgs.nix }:
with import pkgs {};

(haskellPackages.developPackage {
  root = ./.;
}).overrideAttrs (_: {
  SOLVER = "${import ../. args}/bin/mrsep.ilp";
})
