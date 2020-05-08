# motifcluster 0.1.0

## R CMD check results

### Local install

There were no ERRORs, WARNINGs or NOTEs on the local install:

- Arch Linux 5.6.4-arch-1, R-release (R 4.0.0)

### Travis CI

There were no ERRORs, WARNINGs or NOTEs on any of the Travis CI builds:

- Ubuntu 16.04.6 on Travis CI, R-release (R 4.0.0)
- Ubuntu 16.04.6 on Travis CI, R-devel (R under development r78384)
- Ubuntu 16.04.6 on Travis CI, R-oldrel (R 3.6.3)
- macOS High Sierra 10.13.6 on Travis CI, R-release (R 4.0.0)
- macOS High Sierra 10.13.6 on Travis CI, R-oldrel

### Win-Builder

Checks on Win-Builder all return a single NOTE, as expected with a new package:

- Win-Builder, R-release (R 4.0.0)
  - checking CRAN incoming feasibility ... NOTE
    Maintainer: 'William George Underwood <wgu2@princeton.edu>'

- Win-Builder, R-devel (R under development r78369)
  - checking CRAN incoming feasibility ... NOTE
    Maintainer: 'William George Underwood <wgu2@princeton.edu>'

- Win-Builder, R-oldrel (R 3.6.3)
  - checking CRAN incoming feasibility ... NOTE
    Maintainer: 'William George Underwood <wgu2@princeton.edu>'
