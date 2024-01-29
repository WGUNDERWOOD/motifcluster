let
pkgs = import <nixpkgs> { };

in pkgs.mkShell {

    buildInputs = with pkgs; [

        # python
        python3
        python3Packages.pandas
        python3Packages.pytest
        python3Packages.pytest-cov
        python3Packages.numpy
        python3Packages.networkx
        python3Packages.scipy
        python3Packages.scikit-learn
        python3Packages.pylint
        python3Packages.sphinx
        python3Packages.sphinx-rtd-theme
        python3Packages.build
        python3Packages.twine

        # R
        R
        qpdf
        rPackages.igraph
        rPackages.devtools
        rPackages.lintr
        rPackages.covr
        rPackages.filesstrings
        rPackages.RSpectra
        rPackages.mclust

        # julia
        # TODO

       # other
       graphviz
       gprof2dot
       # TODO hexsticker
    ];
}
