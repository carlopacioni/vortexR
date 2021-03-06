% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vortexR_pkg_data.R
\name{pac.run.lhs}
\alias{pac.run.lhs}
\title{Collated results from Vortex scenarios - Pacioni et al. (2017)}
\format{
A named list of two \code{data.frame}s:
  run (153 obs, 7 var), lrun (153 obs, 8 var).
}
\source{
Pacioni, C., Williams, M., Lacy RC, Spencer, P.B.S. and Wayne, A.F. (2017)
 Predators and genetic fitness: key threatening factors for the conservation
 of bettong species. Pacific Conservation Biology.
 \href{https://doi.org/10.1071/PC17002}{DOI:10.1071/PC17002}
}
\description{
Data from Pacioni et al. (2017) used to conduct a sensitivity
 analysis on demographic parameters. Vortex outputs, from the project named
 'Pacioni_et_al' and (Latin Hypercube Sampling) sensitivity test scenario
 'ST_LHS' (.run files), were collated with \code{collate_run}.
}
\examples{
data("pac.run.lhs")
head(pac.run.lhs)
}
