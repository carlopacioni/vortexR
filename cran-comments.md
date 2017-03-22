## Test environments
* GNU/Linux (kernel 4.4.0-59-generic, arch x86_64, Ubuntu 16.04.2 LTS), 
  R version 3.3.3 (RStudio Desktop)
* GNU/Linux (kernel 4.4.0-53-generic, arch x86_64, Ubuntu 14.04.5 LTS), 
  R version 3.3.2 (RStudio Server)
* GNU/Linux (kernel 3.13.0-103-generic, arch x86_64, Ubuntu 12.04.5 LTS), 
  R version 3.3.2 (TravisCI)
* win-builder (devel, release, oldrelease)
* MS Windows (TODO Carlo: insert Windows version of "big laptop"), R version 3.3.3
* MS Windows (TODO Carlo: insert Windows version of "small laptop"), R version 3.3.2

## R CMD check results
R CMD check succeeded in all test environments without errors, warnings, or notes.

## Downstream dependencies
We have checked downstream dependencies of vortexR using devtools::revdep_check().
Results: No errors or warnings found.

## External data
To reduce the installed package's size, we have split the raw data we use in
data handling examples into a separate package `vortexRdata` which we will
submit to CRAN in advance.
