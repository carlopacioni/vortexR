#' Calculate the mean recovery rate (Pacioni et al 2017) and compare
#' scenarios
#'
#' \code{rRec} calculates the mean and standard deviation growth rate between
#' the time \code{yr0} and \code{yrt}, which was defined as 'recovery rate' by
#' Pacioni et al (in press). The function then calculates the strictly
#' standardised mean difference (SSMD, Zhang 2007) for each scenario, and each population
#' contained in the data. \code{rRec} uses this statistic to compare each scenario
#' (providing associated p-values) with a baseline scenario.
#'
#' The means and standard deviations are calculated as:
#' **check how to insert formula in Roxy tags**
#' rRec = sigma(Ni*Mi) / sigma(Ni)
#' (N1*M1+N2*M2+N3*M3)/(N1+N2+N3)
#' SD ={N1*S1+N2*S2+N3*S3}/(N1+N2+N3)
#'
#' Where M is the mean growth rate in each year, N is the sample size (number
#' of simulation runs) and S is the standard deviation.
#'
#' The baseline scenario is selected with the argument \code{scenario}. However,
#' if the simulations were part of a sensitivity testing (as indicated by
#' \code{ST}) then the baseline scenario is selected using the scenario with the
#' suffix '(Base)'.
#'
#' @inheritParams pairwise
#' @inheritParams Ne
#' @inheritParams collate_dat
#' @return A table (\code{data.table}) with the mean rRec and its SD, the SSMD
#' and its associated p-value for each scenario and population
#' @references Zhang, X. D. 2007. A pair of new statistical parameters for quality control
#' in RNA interference high-throughput screening assays. Genomics 89:552-561.
#'
#' Pacioni, C., and Mayer, F. (2017). vortexR: an R package for post Vortex
#' simulation analysis.
#' @import data.table
#' @export
#' @examples
#' # Using Pacioni et al. example data. See ?pac.clas for more details.
#' data(pac.clas)
#' recov <- rRec(pac.clas, project="Pacioni_et_al", scenario="ST_Classic",
#'               ST=TRUE, runs=3, yr0=1, yrt=120, save2disk=FALSE,
#'               dir_out="DataAnalysis/rRec")


rRec <- function(data, project, scenario, ST=FALSE, runs, yr0=1, yrt,
                 save2disk=TRUE, dir_out="DataAnalysis/rRec") {
    ############################################################################
    # Dealing with no visible global variables
    ############################################################################
    scen.name <-NULL
    Year <- NULL
    r.stoch <- NULL
    SD.r. <- NULL
    rruns <- NULL
    SDruns <- NULL
    Scenario <- NULL
    SD <- NULL
    Population <- NULL
    SSMD <- NULL
    base <- NULL
    SDbase <- NULL
    pvalues <- NULL
    J <- NULL
    . <- NULL
    ###########################################################################

  fname <- if (ST) {paste(project, "_", scenario, sep="")} else {project}
  data <- data.table(data)
  if (ST ) scenario <- grep("(Base)", data[, unique(scen.name)], value=TRUE)
  setkey(data, Year)
  rTable <- data[J(yr0:yrt), .(rruns=r.stoch * runs, SDruns=SD.r. * runs),
                by=c("scen.name", "pop.name")]
  rTable <- rTable[, .(rRec=sum(rruns) / (runs * length(yr0:yrt)),
                     SD=sum(SDruns) / (runs * length(yr0:yrt))),
                     by=c("scen.name", "pop.name")]
  setnames(rTable, c("scen.name", "pop.name"), c("Scenario", "Population"))
  setkey(rTable, Scenario)
  Base <- rTable[scenario, .(rRec, SD), by=Population]
  setnames(Base, c("rRec", "SD"), c("base", "SDbase"))
  rTable <- merge(rTable, Base, by="Population")
  rTable[, SSMD := (rRec - base) / sqrt(SD^2 + SDbase^2)]
  rTable[, pvalues := vortexR::pval(SSMD)]
  if (save2disk) df2disk(rTable, dir_out, fname, ".rTable")
  return(rTable)
}
