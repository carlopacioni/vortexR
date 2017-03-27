## Test environments

* GNU/Linux Ubuntu 16.04.2 LTS, 64 bit, R 3.3.3 (RStudio Desktop)
* GNU/Linux Ubuntu 14.04.5 LTS, 64 bit, R 3.3.2 (RStudio Server)
* MS Windows Server 2008, 64 bit, R "release", "devel", "oldrelease" (winbuilder)
* MS Windows 7 Professional SP1, 64 bit, R v3.3.3 (RStudio Desktop)
* MS Windows 10 Home, 64 bit, R v3.3.2 (RStudio Desktop)

## R CMD check results
R CMD check succeeded in all test environments without errors, warnings, or notes
except 

* checking CRAN incoming feasibility ... NOTE
    ```
    Maintainer: ‘Carlo Pacioni <C.Pacioni@Murdoch.edu.au>’

    New submission
    ```
    
    We would like to confirm that this is a new submission.
  
* checking package dependencies ... ERROR
    ```
    Package required and available but unsuitable version: 'vortexRdata'
    ```
    
    Occurs only on winbuilder "devel" 3.4.0, but not on winbuilder "devel" 3.3.4 
    or any other test environments.
    
    The R manual section on [package dependencies](
    https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Package-Dependencies)
    advises that data packages should be a `Suggests`, however this causes an error
    `Namespace dependency not required: vortexRdata`.

## Downstream dependencies
We have checked downstream dependencies of vortexR using `devtools::revdep_check()`.
Results: No errors or warnings found.

## External data
To reduce the installed package's size, we have split the raw data we use in
data handling examples into a separate package `vortexRdata` which we will
submit to CRAN in advance.
