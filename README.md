[![Build Status](https://travis-ci.org/florianm/vortexR.svg?branch=master)](https://travis-ci.org/florianm/vortexR)
[![Coverage Status](https://coveralls.io/repos/carlopacioni/vortexR/badge.svg?branch=master&service=github)](https://coveralls.io/github/carlopacioni/vortexR?branch=master)
[![Documentation Status](https://readthedocs.org/projects/vortexr/badge/?version=latest)](https://readthedocs.org/projects/vortexr/?badge=latest)

# vortexR
An R package for Post Vortex Simulation Analysis.  

Using this package, data of population viability analysis (PVA) that were generated with the software Vortex (Lacy & Pollak 2013), can be collated, plotted and analysed using basic (e.g. pairwise comparisons of scenarios) or more advanced statistics such as fitting regression models.

## Quickstart
Install the package from version control from within R:
```
library(devtools)
install_github("carlopacioni/vortexR", build_vignette=TRUE)
```
If you are on Windows and have not used `devtools` before, then you have to download the Rtools executable file from CRAN webpage and run it. `devtools` can be installed from within R with 
```
install.packages("devtools")
```
Also, if you have problem with the option ```build_vignette=TRUE```, you can just use ```install_github("carlopacioni/vortexR")```.

If you have problem installing and loading ```glmulti``` package and/or ```rJava```. This is most likely because you do not have Java installed or are not using the same version as R. Make sure that if you are using a 64-bit version of R, you also have installed a 64 bit version of Java (most automatic installation via web browser will install a 32-bit version). 

## Documentation
Use `help(vortexR)` `?vortexR` or `??vortexR` to see a broad description of the package.
Use `help(package = "vortexR")` to see the documentations available. Read the vignette for a more comprehensive description of the package. Be aware, if you did not use ```build_vignette=TRUE``` this documentation may not be available from within R. In these cases, download the PDF of these documents. On some platforms, `??vortexR` may not work even if you used ```build_vignette=TRUE```.

## Citation
If you use `vortexR`, please cite:
Pacioni, C., and Mayer, F. (2017). vortexR: an R package for post Vortex simulation analysis.
