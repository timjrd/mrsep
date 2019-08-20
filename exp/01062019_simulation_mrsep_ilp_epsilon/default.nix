args@{ nixpkgs ? import ./pkgs.nix }:
let pkgs = import nixpkgs {}; in
pkgs.stdenvNoCC.mkDerivation rec {
  name = "simulation-epsilon";

  python = (pkgs.python3.withPackages (x: with x; [
    biopython
    notebook
    matplotlib
  ])).override (_:{ ignoreCollisions = true; });
  
  buildInputs = [
    python
    pkgs.tree
  ];
  
  SOLVER = "${import ../../src/mrsep.ilp args}/bin/mrsep.ilp";

  phases = [ "installPhase" ];

  installPhase = ''
    bundle=$out/share/${name}
    mkdir -p $bundle/{sim,src}
    mkdir -p $out/bin
        
    cp ${../19052019_simulation_coverage_statistics/results.tar.xz} $bundle/data.tar.xz
    cp ${../19052019_simulation_coverage_statistics/src}/*.py $bundle/sim
    cp ${./src}/*.py $bundle/src
    
    bin=$out/bin/${name}
    echo '#!${pkgs.bash}/bin/bash' > $bin
    echo "${launcher}" >> $bin
    chmod +x $bin
  '';
  launcher = ''
    PATH=${pkgs.coreutils}/bin:${pkgs.gnutar}/bin:${pkgs.xz}/bin
    if mkdir ${name}; then
      cd ${name}
      tar -xf $bundle/data.tar.xz
      mv results data
      export SOLVER=${SOLVER}
      export PYTHONPATH=$bundle/sim
      ${python}/bin/python3 $bundle/src/main.py
    fi
  '';
}
