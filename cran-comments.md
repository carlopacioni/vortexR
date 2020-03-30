## Test environments
* Ubuntu 16.04.6 LTS, 64 bit, R 3.6.2 (Travis CI)
* MS Windows, 64 bit, R "release", "devel" (winbuilder)
* MS Windows 10 Pro, 64 bit, R v3.6.0 (RStudio Desktop)

## R CMD check results
R CMD check succeeded in all test environments without errors, warnings, or notes,
except one note with devel:
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Carlo Pacioni <C.Pacioni@Murdoch.edu.au>'
which we understand is only to confirm the valid email address of the maintener 
and can be ignored if the latter is correct.

## Downstream dependencies
There are currently no downstream dependencies for this package 
