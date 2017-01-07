dynnet
======

Description
-----------

`dynnet` provides a set of latent space models for static and dynamic
networks. Currently, the package provides the original Euclidean model of Hoff,
Raftery, & Handcock (HRH; 2002) along with a variation of the model that uses
reference nodes as a method of parameter identification. Extensions to the
Projection Model and dynamic networks are planned for the near future.

*This is currently experimental software. You have been warned.*

Installation
------------

There are current no plans to release this package to CRAN. You will instead
want to install directly from GitHub using Hadley Wickham's `devtools` package.

-   Install dependencies

    ``` r
    install.packages(c("Rcpp", "RcppArmadillo", "igraph", "coda"))
    ```

    On Mac, and issue has been found with the `Rcpp` dependency. A description of
    the problem and a work-around are available here:

    [https://github.com/jason-morgan/dynnet/issues/4](Issue 4: Install issue on Mac)

-   Once the dependencies are installed, you should be able to install using
    `devtools`:

    ``` r
    devtools::install_github("jason-morgan/dynnet")
    ```

Problems
--------

Please report any issues you have using GitHub's issue tracker. This package is
in the development stages, and any input and bug reports would be greatly
appreciated.
