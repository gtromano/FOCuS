



# FOCuS

FOCuS, acronym for Functional Online CuSum, is an algorithm for detecting changes in mean in real-time. This is achieved by a recursive update of a piecewise quadratic, whose maximum is the CUSUM test statistic for a change. FOCuS can be applied to settings where either the pre-change mean is known or unknown. Furthermore, FOCuS can detect changes in presence of point outliers.


## Installation and Requirements

### Installing the package

To install the package from Github:


```r
devtools::install_github("gtromano/FOCuS")
library(FOCuS)
```


Alternatively one could clone this repository, and run:


```r
install.packages("./FOCuS", repos = NULL, type = "source")
library(FOCuS)
```


### Requirements for the installation

The packages requires `Rcpp (>= 1.0.5)` with compiler support for the `std` library with the `g++14` standard.


### Bugs and further queries

If any bug should be spotted, feel free to open an issue on this repo. For any information regarding this package, please email the package mantainer: `g` dot `romano` at `lancaster.ac.uk`.


## Usage and examples

After installing, it should be possible to access the full documentation with:


```r
help(FOCuS)
```


### Running FOCuS

FOCuS can run in offline mode (for testing purposes) via:


```r
set.seed(42)
y <- c(rnorm(3e5, 1), rnorm(1e4, 0))
res <- FOCuS(y, 18)
```


Once happy with the parameters, one can run FOCuS online via a call to a data-generating function. This data generating function is expected to return one observation, pulled from the process of interest.

For example:


```r
set.seed(42)
databuffer <- c(rnorm(3e5, 1), rnorm(1e4, 0))

f <- function() {
  out <- databuffer[i]     # simulating a pull from a buffer
  i <<- i + 1
  out
}

i <- 1; res <- FOCuS(f, 18)
```
