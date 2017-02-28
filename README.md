[![Build Status](https://travis-ci.org/carlopacioni/vortexR.svg?branch=master)](https://travis-ci.org/carlopacioni/vortexR)
[![Coverage Status](https://coveralls.io/repos/github/carlopacioni/vortexR/badge.svg?branch=review_CRAN_CP)](https://coveralls.io/github/carlopacioni/vortexR?branch=review_CRAN_CP)
[![Documentation Status](https://readthedocs.org/projects/vortexr/badge/?version=latest)](https://readthedocs.org/projects/vortexr/?badge=latest)

# vortexR
An R package for Post Vortex Simulation Analysis.  

Using this package, data of population viability analysis (PVA) generated with the software Vortex (Lacy & Pollak 2013) can be collated, plotted and analysed using basic (e.g. pairwise comparisons of scenarios) or more advanced statistics (e.g. fitting regression models).

## Quickstart
Install the package from version control from within R:
```
install.packages("devtools")
library(devtools)
install_github("carlopacioni/vortexR", build_vignette=TRUE)
``` 
Fresh Windows installations of R will require [Rtools for Windows](https://cran.r-project.org/bin/windows/Rtools/) to run `devtools`.

If installation with the ```build_vignette=TRUE``` option fails, you can just use ```install_github("carlopacioni/vortexR")```.

## Java-based packages
The packages ```glmulti``` and ```rJava``` require a Java Development Kit (JDK) installed and registered correctly with R. Make sure to install a 64-bit JDK if you are using a 64-bit version of R. 

A typical installation path for Java-based packages like `glmulti` in a Linux-based operating system (here: Ubuntu 16.04 LTS) along the lines of [DigitalOcean's tutorial](https://www.digitalocean.com/community/tutorials/how-to-install-java-on-ubuntu-with-apt-get):

* On the terminal, install Java: `sudo apt-get install openjdk-8-jdk` and verify with `java -version`
* Set default Java installation: `sudo update-alternatives --config java`
* Set default Java compiler: `sudo update-alternatives --config javac`
* Set the environment variable `JAVA_HOME` to your preferred Java installation 
  (here we use /usr/lib/jvm/java-8-openjdk-amd64/jre) by 
  appending the line `JAVA_PATH=/path/to/your/Java/binary` to `/etc/environment`:
  `sudo echo "JAVA_PATH=\"/usr/lib/jvm/java-8-openjdk-amd64/jre\"" >> /etc/environment`, 
  `source /etc/environment` to instantly export the new environment variable `JAVA_HOME`,
  verify with `echo $JAVA_HOME` which should print the `/path/to/your/Java/binary`.
* Register Java with R: `sudo R CMD javareconf`
* In R, install `rJava` with `install.packages("rJava")`
* Install `glmulti` with `install.packages("glmulti")`

Independently of vortexR, a `sudo R CMD javareconf` and possibly installation of `rJava` and Java-usiing packages like `glmulti` will be required after each update of R and / or your Java installation.

A typical installation under Windows could follow:

* Download and install [Java](https://java.com/en/)
* Find "Environment Variables", add variable `JAVA_HOME` 
  ([Windows 10 tutorial](https://javatutorial.net/set-java-home-windows-10), 
  [Windows 7 tutorial](http://www.robertsindall.co.uk/blog/setting-java-home-variable-in-windows/)),
  verify on Command Prompt (Win + r, "cmd", Enter): `echo %JAVA_HOME%`
* Follow [this tutorial](https://support.microsoft.com/en-au/help/3103813/qa-when-i-try-to-load-the-rjava-package-using-the-library-command,-i-get-an-error) to set `CLASSPATH`, `PATH`, `RPATH` and `RHOME` correctly.
* On the Command Prompt, run `R CMD javareconf`
* In R, install `rJava` and `glmulti`

If you have still problem installing and loading ```glmulti``` package and/or ```rJava```. This is most likely because you do not have Java installed or are not using the same version as R. Make sure that if you are using a 64-bit version of R, you also have installed a 64 bit version of Java (most automatic installation via web browser will install a 32-bit version). 

Your mileage may vary depending on your operating system and your versions of Java and R. 
Stack Overflow's [R community](http://stackoverflow.com/questions/tagged/r) is a great source for troubleshooting.

## Documentation
Use `help(vortexR)` `?vortexR` or `??vortexR` to see a broad description of the package.
Use `help(package = "vortexR")` to see the documentations available. Read the vignette for a more comprehensive description of the package. Be aware, if you did not use ```build_vignette=TRUE``` this documentation may not be available from within R. In these cases, download the PDF of these documents. On some platforms, `??vortexR` may not work even if you used ```build_vignette=TRUE```.

## Citation
If you use `vortexR`, please use the citation generated from `citation('vortexR')`.
