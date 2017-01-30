#------------------------------------------------------------------------------#
# Package info
#------------------------------------------------------------------------------#
#' vortexR: an R package for Post Vortex Simulation Analysis
#'
#' \code{vortexR} facilitates Post Vortex Simulation Analysis (PVSA) by offering
#' tools to collate multiple Vortex (v10) output files into one R object, generate
#' plots and conduct basic analysis (e.g. pairwise comparisons of scenarios) and
#' more advanced statistics such as fitting of a Generalised Linear Model (GLM)
#' to investigate the main and the interaction effects of the variables of
#' interest.
#'
#' \code{vortexR} has a number of functions that are useful during the
#' development of a Vortex project and for its analysis after completion.
#' \code{vortexR} makes it easy to automatise the creation of plots and
#' computation of basic statistics to inspect the effect of changes carried out
#' in the Vortex project. Once the project development is completed, the same
#' framework used in \code{vortexR} during the development of the project can be
#' refined and extended to include more advanced statistical analyses or can be
#' easily included in Markdown documents for the creation of reports (by
#' converting them into pdf) or update web-pages.
#'
#' The use of \code{vortexR} ensures reproducibility and standardises analytical
#' approaches in population viability analysis.
#'
#' @section Documentations:
#' Use \code{help(package = "vortexR")} for a list of \code{vortexR} functions
#' and their specific documentations.
#'
#' A more detailed description of the package and functions can be opened with:
#' \code{vignette(package="vortexR", topic="User-guide")}.
#'
#' More vignettes may be come available in the future. Use
#' \code{vignette(package="vortexR")} to see all the available vignettes.
#'
#' @section Citation:
#' If you use \code{vortexR}, please cite:
#' Pacioni, C., and Mayer, F. (in prep). vortexR: an R package for post Vortex
#' simulation analysis.
#'
#' @section Get in touch:
#' Please, use \url{https://github.org/carlopacioni/vortexR/issues} to report
#' any issues with vortexR. If unsure, or for feedback, contact me at:
#' carlo.pacioni 'at' gmail.com.
#'
#' @section Publications:
#' \emph{Below there are listed a few publications that used \code{vortexR}.}
#'
#' Campbell et al. (2016). Assessing the economic benefits of starling
#'  detection and control to Western Australia. Australasian Journal of
#'  Environmental Management, 23, 81-99.
#'
#' Pacioni, C., Williams, M., Lacy RC, Spencer, P.B.S. and Wayne, A.F. (in prep)
#'  Predators and genetic fitness: key threatening factors for the conservation
#'  of bettong species. Pacific Conservation Biology.
#'
#' @docType package
#'
#' @name vortexR
NULL



#------------------------------------------------------------------------------#
# Data
#------------------------------------------------------------------------------#

#' @name sta.main
#' @title Collated results from Vortex scenarios - Campbell et al (2016)
#' @description A dataset with the results from the main Vortex scenarios used
#'  in Campbell et al (2016). Vortex outputs, from the project named
#'  "Starlingv3PopBased" (.dat files), were collated with \code{collate_dat}.
#' @usage data(sta.main)
#' @format a \code{data.frame} with 1632 observations of 44 variables.
#' @source Campbell et al. (2016). Assessing the economic benefits of starling
#'  detection and control to Western Australia. Australasian Journal of
#'  Environmental Management, 23, 81-99
NULL

#' @name sta.evy5
#' @title Collated results from Vortex scenarios - Campbell et al (2016)
#' @description A dataset with the results from Vortex scenarios used in
#'  Campbell et al (2016) to simulate major application of control measures
#'  in every 5 year cycle. Vortex outputs, from the project named
#'  "Starlingv3PopBased" and the sensitivity test scenario "MReductEvy5"
#'  (.stdat files), were collated with \code{collate_dat}.
#' @usage data(sta.evy5)
#' @format a \code{data.frame} with 1020 observations of 47 variables.
#' @source Campbell et al. (2016). Assessing the economic benefits of starling
#'  detection and control to Western Australia. Australasian Journal of
#'  Environmental Management, 23, 81-99
NULL

#' @name sta.evy5.b11
#' @title Collated results from Vortex scenarios - Campbell et al (2016)
#' @description A dataset with the results from Vortex scenarios used in Campbell
#'  et al (in press) to simulate major application of control measures in every
#'  5 year cycle, maintaining 2011 levels of investment. Vortex outputs, from
#'  the project named "Starlingv3PopBased" and the sensitivity test scenario
#'  "MReduction_B11_09Evy5" (.stdat files), were collated with \code{collate_dat}.
#' @usage data(sta.evy5.b11)
#' @format a \code{data.frame} with 1020 observations of 47 variables.
#' @source Campbell et al. (2016). Assessing the economic benefits of starling
#'  detection and control to Western Australia. Australasian Journal of
#'  Environmental Management, 23, 81-99
NULL

#' @name pac.clas
#' @title Collated results from Vortex scenarios - Pacioni et al. (in prep)
#' @description Subset (only 3 runs) of data from Pacioni et al. (in prep) used
#'   to conduct a sensitivity analysis on demographic parameters. Vortex
#'   outputs, from the project named "Pacioni_et_al" and (Single-Factor)
#'   sensitivity test scenario "ST_Classic" (.stdat files), were collated with
#'   \code{collate_dat}.
#' @usage data(pac.clas)
#' @format a \code{data.frame} of 2904 observations of 68 variables.
#' @source Pacioni, C., Williams, M., Lacy RC, Spencer, P.B.S. and Wayne, A.F. (in prep)
#'   Predators and genetic fitness: key threatening factors for the conservation
#'   of bettong species. Pacific Conservation Biology.
NULL

#' @name pac.lhs
#' @title Collated results from Vortex scenarios - Pacioni et al. (in prep)
#' @description Data from Pacioni et al. (in prep) used to conduct a sensitivity
#'  analysis on demographic parameters. Vortex outputs, from the project named
#'  "Pacioni_et_al" and (Latin Hypercube Sampling) sensitivity test scenario
#'  "ST_LHS" (.stdat files), were collated with \code{collate_dat}.
#' @usage data(pac.lhs)
#' @format A \code{data.frame} of 6171 observations of 68 variables.
#' @source Pacioni, C., Williams, M., Lacy RC, Spencer, P.B.S. and Wayne, A.F. (in prep)
#'  Predators and genetic fitness: key threatening factors for the conservation
#'  of bettong species. Pacific Conservation Biology.
NULL

#' @name pac.run.lhs
#' @title Collated results from Vortex scenarios - Pacioni et al. (in prep)
#' @description Data from Pacioni et al. (in prep) used to conduct a sensitivity
#'  analysis on demographic parameters. Vortex outputs, from the project named
#'  "Pacioni_et_al" and (Latin Hypercube Sampling) sensitivity test scenario
#'  "ST_LHS" (.run files), were collated with \code{collate_run}.
#' @usage data(pac.run.lhs)
#' @format A named list of two \code{data.frame}s:
#'   run (153 obs, 7 var), lrun (153 obs, 8 var).
#' @source Pacioni, C., Williams, M., Lacy RC, Spencer, P.B.S. and Wayne, A.F. (in prep)
#'  Predators and genetic fitness: key threatening factors for the conservation
#'  of bettong species. Pacific Conservation Biology.
NULL

#' @name pac.yr
#' @title Collated results from Vortex scenarios - Pacioni et al. (in prep)
#' @description Data from Pacioni et al. (in prep) used to conduct a sensitivity
#'  analysis on demographic parameters. Vortex outputs, from the project named
#'  "Pacioni_et_al" and (Single-Factor) sensitivity test scenario "ST_Classic"
#'  (.yr files), were collated with \code{collate_yr}.
#' @usage data(pac.yr)
#' @format A named list of two \code{element}s:
#'   all (8712 obs, 26 var), means (2904 obs, 25 var).
#' @source Pacioni, C., Williams, M., Lacy RC, Spencer, P.B.S. and Wayne, A.F. (in prep)
#'  Predators and genetic fitness: key threatening factors for the conservation
#'  of bettong species. Pacific Conservation Biology.
NULL

#' @name pac.clas.Nadults
#' @title Harmonic mean of adults and population sizes
#' @description Data from Pacioni et al. (in prep) - sensitivity test scenario
#' "ST_Classic" - were used to calculate the harmonic mean of adults and population
#' sizes using \code{Nadults}.
#' @usage data(pac.clas.Nadults)
#' @format A \code{data.frame} with 24 observations of 4 variables.
#' @source Pacioni, C., Williams, M., Lacy RC, Spencer, P.B.S. and Wayne, A.F. (in prep)
#'  Predators and genetic fitness: key threatening factors for the conservation
#'  of bettong species. Pacific Conservation Biology.
NULL

#' @name pac.clas.Ne
#' @title Effective population size
#' @description Data from Pacioni et al. (in prep) - sensitivity test scenario
#' "ST_Classic" - were used to calculate the effective population size
#' sizes using \code{Ne}.
#' @usage data(pac.clas.Ne)
#' @format A \code{data.frame} with 24 observations of 2 variables.
#' @source Pacioni, C., Williams, M., Lacy RC, Spencer, P.B.S. and Wayne, A.F. (in prep)
#'  Predators and genetic fitness: key threatening factors for the conservation
#'  of bettong species. Pacific Conservation Biology.
NULL

#' @name pac.clas.lookup
#' @title Look-up table
#' @description Data from Pacioni et al. (in prep) - sensitivity test scenario
#' "ST_Classic" - were used to generate a look-up table
#' sizes using \code{lookup_table}.
#' @usage data(pac.clas.lookup)
#' @format A \code{data.frame} with 24 observations of 8 variables.
#' @source Pacioni, C., Williams, M., Lacy RC, Spencer, P.B.S. and Wayne, A.F. (in prep)
#'  Predators and genetic fitness: key threatening factors for the conservation
#'  of bettong species. Pacific Conservation Biology.
NULL

#' @name pac.clas.pairw
#' @title Results of pairwise comparisons of simulation scenarios
#' @description Results of pairwise comparisons of simulation scenarios included
#' in the sensitivity test scenario "ST_Classic" using \code{pairwise}.
#' @usage data(pac.clas.pairw)
#' @format A named list of 12 \code{element}s. See documentation for details.
#' @source Pacioni, C., Williams, M., Lacy RC, Spencer, P.B.S. and Wayne, A.F. (in prep)
#'  Predators and genetic fitness: key threatening factors for the conservation
#'  of bettong species. Pacific Conservation Biology.
NULL

#------------------------------------------------------------------------------#
# Data - Raw Vortex outputs
#------------------------------------------------------------------------------#
#' @name pacioni
#' @title Raw Vortex outputs from Pacioni et al. (in prep)
#' @description A folder with Vortex outputs from Pacioni et al. (in prep) used
#'  to run examples and Vortex project file. NOTE: these simulations are shorter
#'  than those presented in the paper (only 3 runs for 120 'Vortex-years').
#' @usage
#' # To retrieve the path to the files use:
#' # system.file("extdata", "pacioni", package="vortexR")
#' @format One .xml file and several .run and .stdat files.
#' @source Pacioni, C., Williams, M., Lacy RC, Spencer, P.B.S. and Wayne, A.F. (in prep)
#'  Predators and genetic fitness: key threatening factors for the conservation
#'  of bettong species. Pacific Conservation Biology.
NULL

#' @name campbell
#' @title Raw Vortex outputs from Campbell et al (2016)
#' @description A folder with Vortex outputs from Campbell et al (2016).
#' @usage
#' # To retrieve the path to the files use:
#' # system.file("extdata", "campbell", package="vortexR")
#' @format Several .dat and .stdat files.
#' @source Campbell et al. (2016). Assessing the economic benefits of starling
#'  detection and control to Western Australia. Australasian Journal of
#'  Environmental Management, 23, 81-99
NULL
