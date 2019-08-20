{ pkgs ? import ./pkgs.nix }:
with import pkgs {};

stdenv.mkDerivation {
  name = "art";
  src = fetchurl {
    url = "https://www.niehs.nih.gov/research/resources/assets/docs/artsrcmountrainier2016.06.05linux.tgz";
    sha256 = "127lxzgs9i8a1wwai13kicsacz2b4ya55bis0kg4if2fi1hdxbk9";
  };
  buildInputs = [
    gsl
  ];
  preBuild = "make clean";
}
