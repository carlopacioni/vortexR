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
    
    This ERROR only occurs on winbuilder "devel" 3.4.0, but not on winbuilder 
    "devel" 3.3.4 or any other test environments.
    
## Reviewer notes
### Version 1.0.3
The package was submitted as v1.0.2 to CRAN. The reviewer suggested:
* Description: Facilitate Post Vortex Simulation Analysis by offering tools 
   to collate multiple Vortex (v10) output files into one R 
   object, and analyse the collated output statistically. 

  For non-experts, can you pls elaborate a bit more what 'Vortex' is about?            Perhaps add something like

  See <http...> for more information about ...

The field `Description` in `DESCRIPTION` has been modified accordingly.

## Downstream dependencies
We have checked downstream dependencies of vortexR using `devtools::revdep_check()`.
Results: No errors or warnings found.

## External data
To reduce the installed package's size, we have split the raw data we use in
data handling examples into a separate package `vortexRdata` which is now available from CRAN.
