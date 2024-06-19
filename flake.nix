{
  description = "motifcluster";
  inputs = {
    nixpkgs.url = github:NixOS/nixpkgs/nixos-24.05;
    flake-utils.url = github:numtide/flake-utils;
  };
  outputs = {
    self,
    nixpkgs,
    flake-utils,
  }:
    with flake-utils.lib;
      eachSystem allSystems (system: let
        pkgs = nixpkgs.legacyPackages.${system};
        python = pkgs.python3.withPackages (ps:
          with ps; [
            build
            matplotlib
            networkx
            numpy
            pandas
            pylint
            pytest
            pytest-cov
            scikit-learn
            scipy
            seaborn
            sphinx
            sphinx-rtd-theme
            twine
          ]);
        R = pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [
            covr
            devtools
            filesstrings
            igraph
            lintr
            mclust
            RSpectra
          ];
        };
        julia = pkgs.julia_19.withPackages [
          "Aqua"
          "Colors"
          "Coverage"
          "CSV"
          "DataFrames"
          "Dates"
          "Distributions"
          "Documenter"
          "JuliaFormatter"
          "Plots"
          "PyPlot"
          "Random"
          "Suppressor"
          "Test"
        ];
        hexsticker = pkgs.python3Packages.buildPythonPackage rec {
          pname = "hexsticker";
          version = "1.2.0";
          src = pkgs.python3Packages.fetchPypi {
            inherit pname version;
            sha256 = "sha256-HjKMjWUqBe2Cnmb1a33IqfIPqRoWxbD6kU9/QjaSaWs=";
          };
          postPatch = ''
            echo "click" > requirements.txt
            echo "pillow" >> requirements.txt
            sed -i "s/ANTIALIAS/LANCZOS/g" hexsticker/create.py
          '';
          doCheck = false;
          propagatedBuildInputs = [
            pkgs.python3Packages.click
            pkgs.python3Packages.pillow
          ];
        };
      in {
        devShells.default = pkgs.mkShell {
          buildInputs = [
            julia
            python
            R
            hexsticker
            pkgs.qpdf
            pkgs.graphviz
            pkgs.gprof2dot
            pkgs.texlive.combined.scheme-full
            pkgs.imagemagick
            pkgs.ghostscript
            pkgs.optipng
          ];
        };
      });
}
