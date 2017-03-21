## Test environments
* GNU/Linux (kernel 4.4.0-53-generic, Ubuntu 14.04.5 LTS), R version 3.3.2 (RStudio Server)
* GNU/Linux (kernel 3.13.0-103-generic, Ubuntu 12.04.5 LTS), R version 3.3.2 (TravisCI)
* win-builder (devel and release) TODO run this
* MS Windows (TODO Carlo: insert windows and R version)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking installed package size ... NOTE
  installed size is  5.3Mb
  sub-directories of 1Mb or more:
    extdata   4.0Mb

  The included raw data files in extdata serve as the core example we discuss in 
  the paper submitted to Methods of Ecology and Evolution. 
  The paper contains usage instructions which unfortunately require the raw data 
  to be included in the package, rather than provided via an auxiliary data package.

## Downstream dependencies
We have checked downstream dependencies of vortexR using devtools::revdep_check().
Results: No errors or warnings found.
