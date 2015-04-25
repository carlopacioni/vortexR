# vortexR
An R package for Post Vortex Simulation Analysis

## Quickstart
Install the package from version control from within R:
```
library(devtools)
install_github("carlopacioni/vortexR")
```
If you are on Windows and have not used `devtools` before, then you have to download the Rtools executable file from CRAN webpage and run it. `devtools` can be installed from within R with 
```
install.packages("devtools")
```
Also, if you have problem installing and loading ```glmulti``` package and/or ```rJava```. This is most likely because you do not have Java installed or are not using the same version as R. Make sure that if you are using a 64-bit version of R, you also have installed a 64 bit version of Java (most automatic installation via web browser will install a 32-bit version). 

## Documentation
Use `help(vortexR)` or `?vortexR` to see a broad description of the package.
Use `help(package = "vortexR")` to see the documentations available.
Read the vignette for a more comprehensive description of the package.

## Citation
If you use `vortexR`, please cite:
Pacioni, C., and Mayer, F. (in prep). vortexR: an R package for post Vortex simulation analysis.
