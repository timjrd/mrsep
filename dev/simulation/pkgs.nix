# nixos-18.09 on Apr 3, 2019 
let rev    = "1d36ad6d16dbf1d3937f899a087a4360332eb141";
    sha256 = "0rf1n61xlbvanrknh7g9884qjy6wmwc5x42by3f9vxqmfhz906sq";
in import (fetchTarball {
  url = "https://github.com/NixOS/nixpkgs-channels/archive/${rev}.tar.gz";
  inherit sha256;
}) {}
