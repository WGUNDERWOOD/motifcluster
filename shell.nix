let
pkgs = import <nixpkgs> { };

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

in pkgs.mkShell {

    buildInputs = with pkgs; [

        # python
        python3
        python3Packages.build
        python3Packages.matplotlib
        python3Packages.networkx
        python3Packages.numpy
        python3Packages.pandas
        python3Packages.pylint
        python3Packages.pytest
        python3Packages.pytest-cov
        python3Packages.scikit-learn
        python3Packages.scipy
        python3Packages.seaborn
        python3Packages.sphinx
        python3Packages.sphinx-rtd-theme
        python3Packages.twine

        # R
        R
        qpdf
        rPackages.covr
        rPackages.devtools
        rPackages.filesstrings
        rPackages.igraph
        rPackages.lintr
        rPackages.mclust
        rPackages.RSpectra

        # julia
        julia

        # other
        graphviz
        gprof2dot
        texlive.combined.scheme-full
        imagemagick
        ghostscript
        hexsticker
        optipng
    ];
}#
