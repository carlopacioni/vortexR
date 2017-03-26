[![Build Status](https://travis-ci.org/carlopacioni/vortexR.svg?branch=master)](https://travis-ci.org/carlopacioni/vortexR)
[![Coverage Status](https://coveralls.io/repos/github/carlopacioni/vortexR/badge.svg?branch=master)](https://coveralls.io/github/carlopacioni/vortexR?branch=master)

# vortexR
An R package for Post Vortex Simulation Analysis.  

Using this package, data of population viability analysis (PVA) generated with 
the software Vortex (Lacy & Pollak 2013) can be collated, plotted and analysed 
using basic (e.g. pairwise comparisons of scenarios) or more advanced statistics 
(e.g. fitting regression models).

## Install
To install the latest stable release from CRAN:

```
install.packages("vortexR")
```

To install the latest development version from source:

```
install.packages("devtools")
library(devtools)
install_github("carlopacioni/vortexR", build_vignette=TRUE)
``` 
Fresh Windows installations of R will require 
[Rtools for Windows](https://cran.r-project.org/bin/windows/Rtools/).

If installation with ```build_vignette=TRUE``` fails, 
you can run ```install_github("carlopacioni/vortexR")```.

### Java-based packages
The packages ```glmulti``` and ```rJava``` require a Java Development Kit (JDK) 
installed and registered correctly with R. 
Make sure to install a 64-bit JDK if you are using a 64-bit version of R. 

A typical installation path for Java-based packages like `glmulti` in a 
GNU/Linux-based operating system (here: Ubuntu 16.04 LTS) along the lines of 
[DigitalOcean's tutorial](https://www.digitalocean.com/community/tutorials/how-to-install-java-on-ubuntu-with-apt-get):

* On the terminal, install Java: `sudo apt-get install openjdk-8-jdk` and verify with `java -version`
* Set default Java installation: `sudo update-alternatives --config java`
* Set default Java compiler: `sudo update-alternatives --config javac`
* Set the environment variable `JAVA_HOME` to your preferred Java installation: 
  (here we use /usr/lib/jvm/java-8-openjdk-amd64/jre) by appending the line 
  `JAVA_PATH=/path/to/your/Java/binary` to `/etc/environment`:
  `sudo echo "JAVA_PATH=\"/usr/lib/jvm/java-8-openjdk-amd64/jre\"" >> /etc/environment`
* Run `source /etc/environment` to instantly export the new environment variable `JAVA_HOME`
* Verify that `echo $JAVA_HOME` prints the `/path/to/your/Java/binary`
* Register Java with R: `sudo R CMD javareconf`
* In R, install `rJava` with `install.packages("rJava")`
* Install `glmulti` with `install.packages("glmulti")`

Independently of vortexR, a `sudo R CMD javareconf` (and possibly the installation 
of `rJava` and Java-using packages like `glmulti`) will be required after each 
update of R and / or your Java installation.

A typical installation under Windows could follow:

* Download and install [Java](https://java.com/en/)
* Find "Environment Variables", add variable `JAVA_HOME` 
  ([Windows 10 tutorial](https://javatutorial.net/set-java-home-windows-10), 
  [Windows 7 tutorial](http://www.robertsindall.co.uk/blog/setting-java-home-variable-in-windows/)),
  verify on Command Prompt (Win + r, "cmd", Enter): `echo %JAVA_HOME%`
* Make sure that the environment variable 'JAVA_HOME' is set properly and points to a 64-bit version jof Java(JDK or SE
  Runtime). For 64-bit Revolution R you need to make sure that you have a 64-bit version of Java installed.
* Make sure that the Java 'CLASSPATH' variable is set properly.
  For the 64-bit version of the JDK, this could be e.g. `C:\Program Files\Java\jdk1.7.0_02\jre\lib\ext`
* Check your system PATH settings. On Windows, goto `Control Panel ---> System --> Advanced Settings ---> Environment Variables`.
  Your path should include the directories:   `C:\Program Files\Java\jre6\bin\server` 
  and the R installation directory, e.g. `C:\Revolution\R-Enterprise-7.3\R-3.1.1\bin\x64`.
* Make sure that the environment variables `RPATH` and `RHOME` are BOTH set and point to the R installation directory, e.g.
  `C:\Revolution\R-Enterprise-7.3\R-3.1.1\bin\x64`
* In R, install `rJava` and `glmulti`

If you have still problem installing and loading `glmulti` package and/or 
`rJava`, you may not have Java installed or are not using the same Java version as R. 
Make sure that if you are using a 64-bit version of R, you also have installed 
a 64 bit version of Java (most automatic installation via web browser will 
install a 32-bit version). 

Your mileage may vary depending on your operating system and your versions of Java and R. 
Stack Overflow's [R community](http://stackoverflow.com/questions/tagged/r) 
is a great source for troubleshooting.

## Learn
Use `help(vortexR)` `?vortexR` or `??vortexR` to see a broad description of the package.
Use `help(package = "vortexR")` to see the documentations available. 
Read the vignette for a more comprehensive description of the package. 
Be aware, if you did not use `build_vignette=TRUE` this documentation may 
not be available from within R. In these cases, download the PDF of these documents. 
On some platforms, `??vortexR` may not work even if you used `build_vignette=TRUE`.

## Use
If you use `vortexR`, please use the citation generated from `citation('vortexR')`.

## Contribute
We are happy to receive feedback and contributions through bug reports and pull requests.

Notes on vortexR's code style:

* Most of the code is formatted following Google's R Style Guide. 
* Some parts have been then re-formatted following Hadley's style suggestions from his book:
  r-packages. The latter re-formatting mainly concerns function and parameter names.
* Exported functions and arguments should use only snake_case, avoid CamelCase and dot.case.
* CamelCase and dot.case are OK in the code.
* Keep a space around operators except with `=` for function parameters.
* "if" statements only go in one line IF there is one command only (and with 
  curly brackets), otherwise split them over multiple lines.
* No full stop in the title and arguments' descriptions within roxygen's comments.
