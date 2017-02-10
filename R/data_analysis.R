#------------------------------------------------------------------------------#
# Data analysis
#------------------------------------------------------------------------------#

#' Calculate the effective population size (Ne)
#'
#' \code{Ne} calculates the effective population size (Ne) between \code{yr0}
#' and \code{yrt} for several scenarios based on the loss of genetic diversity
#' (expected heterozygosity) using the temporal approach.
#'
#' \code{yr0} is adjusted by adding the number of years of the generation time
#' (rounded to the nearest integer).  In this way the user can provide the same
#' \code{yr0,yrt} and \code{gen} to \code{Nadults} and \code{Ne} and these values
#' are adjusted internally to correctly calculate the Ne/N ratios where relevant.
#' If this behaviour is not desired, use \code{gen=0}.
#'
#' \strong{NOTE:} When a population goes extinct, the results of the calculations
#' are spurious (they are 0.5). This may change in future versions.
#'
#' @param data The output from \code{collate_dat}
#' @param scenarios A vector of scenario names for which Ne needs to be calculated,
#' default: "all"
#' @param gen The generation time express in years
#' @param yr0,yrt The time window to be considered (first and last year respectively)
#' @param save2disk Whether to save the output to disk, default: TRUE
#' @param fname The name of the files where to save the output, defult: "Ne"
#' @param dir_out The local path to store the output. Default: DataAnalysis
#' @return A \code{data.table} (\code{data.frame} if \code{\link[data.table]{data.table}} is not
#'  loaded) with Ne values
#' @import data.table
#' @export
#' @examples
#' # Using Pacioni et al. example data. See ?pac.clas for more details.
#' data(pac.clas)
#' # Calculate Ne for all scenarios in the data. Note the odd value for scenario
#' # 12, consequent to the population going extinct.
#' NeAll <- Ne(data=pac.clas, scenarios="all", gen=2.54, yr0=50, yrt=120,
#'          save2disk=FALSE)
Ne <-  function(data=NULL,
                scenarios="all",
                gen=1,
                yr0=1,
                yrt=2,
                save2disk=TRUE,
                fname="Ne",
                dir_out="DataAnalysis") {

  # Function definitions
  NeOne <- function(scenario) {
          stdatDT <- data.table(data)
          setkeyv(stdatDT, c("scen.name", "Year"))

          # From yr0+gen
          t0 <- round(yr0 + gen)
          tf <- yrt
          Gt0 <- stdatDT[J(scenario, t0), "Het", with=FALSE]
          Gtf <- stdatDT[J(scenario, tf), "Het", with=FALSE]
          t <- (tf - t0) / gen
          sq <- (Gtf[ , Het] / Gt0[ , Het] )^(1 / t)
          den <- sq - 1
          effPopSize <- -(1/den)*(1/2)
          effPopSize <- t(as.data.frame(effPopSize))
          effPopSize <- as.data.frame(effPopSize)
          names(effPopSize) <- stdatDT[J(scenario, t0), pop.name]
          return(effPopSize)
        }
  if (scenarios == "all")
    scenarios <- unique(data$scen.name)

  r.Ne <- lapply(scenarios, NeOne)
  r.Ne <- rbindlist(r.Ne)
  r.Ne <- r.Ne[ , Scenario := scenarios]

  if (save2disk) {
    df2disk(r.Ne, dir_out, fname)
  }

  message(paste("Effective population size based on loss of gene diversity from year",
                round(yr0 + gen), "to year", yrt))
  message("NOTE: The first year used in the calculations is adjusted using the generation time provided (yr0 + gen).
See documentation for more information")
  return(r.Ne)
}

#' Calculate the harmonic mean of the total number of adults
#'
#' \code{Nadults} calculates, for several scenarios, the harmonic mean of the total
#' number of adults between \code{yr0} and \code{yrt}. These can be use to
#' calculate Ne/N ratios where relevant.
#'
#' \code{yrt} is adjusted by subtracting the number of years of the generation
#' time (rounded to the nearest integer). In this way the user can provide the
#' same \code{yr0,yrt} and \code{gen} to \code{Nadults} and \code{Ne} and these
#' values are adjusted internally to correctly calculate the Ne/N ratios where
#' relevant. If this behaviour is not desired, use \code{gen=0}.
#'
#' @param data The second element (census_means) of the output from \code{collate_yr}
#' @param npops_noMeta The total number of populations excluding the metapopulation,
#' default: 1
#' @param appendMeta Whether to calculate data for the metapopulation,
#' default: FALSE
#' @param fname The name of the files where to save the output, defult: "Nadults"
#' @inheritParams Ne
#' @return A \code{data.table} (\code{data.frame} if \code{\link[data.table]{data.table}} is not
#'  loaded) with Nb values
#' @import data.table
#' @export
#' @examples
#' # Using Pacioni et al. example data. See ?pac.yr for more details.
#' data(pac.yr)
#' NadultAll <- Nadults(data=pac.yr[[2]], scenarios="all", gen=2.54, yr0=50, yrt=120,
#'              save2disk=FALSE)
Nadults <- function (data,
                scenarios="all",
                npops_noMeta=1,
                appendMeta=FALSE,
                gen=1,
                yr0=1,
                yrt=2,
                save2disk=TRUE,
                fname="Nadults",
                dir_out="DataAnalysis") {

  # Function definitions
  HarmMean <- function(x) 1/mean(1/(x))

  LongFormat <- function (numPop) {
    k <- numPop - 1
    cstart <- 3 + (k * ncolpop)
    cend <- 3 + ncolpop * numPop - 1
    pop <- rep(paste0("pop", numPop), nrow(data))
    one.pop <- data[ , list(Scenario,Year)]
    one.pop <- one.pop[ , "Population" := pop]
    one.pop[ , (h1) := data[ , cstart:cend, with=FALSE]]
    return(one.pop)
  }

  # Convert in long format

  h <- names(data)

  # Number of cols for each pop except the first 2, which are fixed cols
  ncolpop <- (length(h) - 2) / npops_noMeta
  # Headings without pop suffix
  h1 <- gsub(pattern="pop1", "", h[3:(3 + ncolpop - 1)])

  # Reshape df in long format
  tldata <- lapply(1:npops_noMeta, LongFormat)

  lcensusMeans <- rbindlist(tldata)

  # If more than one pop, do Metapopulation calculations and rbind
  if (npops_noMeta > 1 & appendMeta) {
    message("Doing calculations for Metapopulation...")
    message("Please wait...")

    census <- c("N", "AM", "AF", "Subadults", "Juv",
                "nDams", "nBroods", "nProgeny")
    meta <-lcensusMeans[, lapply(.SD, sum), by=list(Scenario,Year),
                        .SDcols=census]
    message("Done...")
    message("Appending Metapopulation data to CensusMeans data frame...")
    gs <- names(lcensusMeans)[grep(pattern="^GS", names(lcensusMeans))]
    meta[ , Population := rep(paste0("pop", npops_noMeta + 1), nrow(data))]
    setcolorder(meta, c("Scenario", "Year", "Population", census))
    meta[ , (gs) := tldata[[1]][ , gs, with=FALSE]]

    lcensusMeans <- rbindlist(list(lcensusMeans, meta),
                              use.names=TRUE, fill=TRUE)
    message("Done!")
  }

  # Calculate harmonic means
  message(paste("Calculating the harmonic means of total number of individuals",
                "and number of adults from year", yr0, "to year", round(yrt-gen)))
  message("NOTE: The last year used in the calculations is adjusted using the generation time provided (yrt - gen).
          See documentation for more information")

  if (scenarios == "all")
    scenarios <- data[ , unique(Scenario)]

  setkey(lcensusMeans, Scenario)
  slcensusMeans <- lcensusMeans[J(scenarios), ]
  setkey(slcensusMeans, Year)
  slcensusMeans <- slcensusMeans[.(yr0:round(yrt - gen)), ]

  slcensusMeans[ , Nad := sum(.SD), .SDcols=c("AM", "AF"),
                by=list(Scenario,Population,Year)]

  harm.means <- slcensusMeans[ , lapply(.SD, HarmMean), .SDcols=c("Nad", "N"),
                              by=list(Scenario,Population)]
  message("Done!")

  if (save2disk) {df2disk(harm.means, dir_out, fname)}
  return(harm.means)
}

#' Pairwise comparisons and ranks of scenarios
#'
#' \code{pairwise} conducts pairwise comparisons against a baseline scenario
#' using sensitivity coefficients and strictly standardised mean difference. It also
#' ranks scenarios (and/or parameters when relevant) using these statistics.
#'
#' When \code{yrs="max"} (default), VortexR automatically sets  \code{yrs} to
#' the last year of the simulation.
#'
#' Pairwise comparisons against a baseline scenario are conducted using
#' sensitivity coefficients (SC, Drechsler et al. 1998) and strictly
#' standardised mean difference (SSDM, Zhang 2007).
#'
#' \code{pairwise} ranks, for each population, the scenarios (and SVs if
#' relevant, see below) based on the absolute value of the statistics (either SC
#' or SSMD) regardless of the sign. That is, the scenario with the absolute SC
#' or SSMD value most different from zero will have a rank equal to '1'. The
#' actual statistics need to be inspected to evaluate the direction of the change.
#'
#' The Kendall's coefficient of concordance is calculated to test whether the
#' order of ranked scenarios (or SVs if relevant) is statistically consistent
#' across the chosen points in time and parameters (or SVs). For example, if 100
#' years were simulated, \code{yrs=c(50, 100)} and \code{params=c("Nall", "Het")},
#' the consistency of ranking will be tested across the four raters (i.e. Nall
#' at year 50, and at year 100, Het at year 50 and at year 100). Kendall's test
#' operates a listwise deletion of missing data. However, when data in a whole
#' column (i.e. ranks for a parameter) are missing, the column is removed before
#' the statistic is calculated (See vignette for more information).
#'
#' It is possible to evaluate the mean effect of a range of values for certain
#' parameters on their outcome variables of interest (i.e. ranking the parameters,
#' rather than scenarios). This is automatically done when the analysis is
#' conducted on with \code{ST=TRUE,type="Single-Factor"} and there is more than
#' one SV passed with the argument \code{SVs}. Alternatively, it is achievable
#' with a combined use of \code{group.mean=TRUE,SVs}. The first argument result
#' in the calculations, following Conroy and Brook (2003), of the mean SC and
#' SSMD for each group of scenarios that have different parameter values. \code{SVs}
#' provides the names of the parameters to be considered. Parameters are then
#' ranked accordingly (See vignette for more information).
#'
#' The parameter values passed with \code{SVs} are evaluated at year=0. This is
#' done because these parameters may take value 'zero' if the relevant
#' populations goes extinct. There are cases where Vortex may not evaluate these
#' parameters even at year 0. This may happen, for example, when a population is
#' empty at initialization (i.e. the initial population size is zero), or when K
#' is set to zero at the beginning of the simulation. The user has to make sure
#' that the values for the parameters passed in are correct.
#'
#' Note that it only makes sense to rank parameters in a ST run when the
#' Single-Factor option is used in Vortex. This is because with Single-Factor,
#' the parameters are modified one at the time (See vignette for more information).
#'
#' @param data A data.frame generated by \code{collate_dat}
#' @param project The Vortex project name
#' @param scenario The ST Vortex scenario name or the scenario that should be
#' used as baseline if simulations were not conducted with the ST module
#' @param params A character vector with the parameters to be compared,
#'  default: c("PExtinct", "Nextant", "Het", "Nalleles")
#' @param yrs The year(s) to be analysed, default: "max"
#' @param ST Whether files are from sensitivity analysis (TRUE),
#'  or not (FALSE, default)
#' @param type Type of ST. Possible options are: "Sampled",
#' "Latin Hypercube Sampling", "Factorial" or "Single-Factor"
#' @param group.mean Whether calculate the mean of the statistics
#'  (SSMD and Sensitivity Coefficient) by group. See details
#' @param SVs A character vector with the parameters to be used to group
#' scenarios, default: NA
#' @param save2disk Whether to save the output to disk, default: TRUE
#' @param dir_out The local path to store the output. Default: DataAnalysis/Pairwise
#' @return A list of six elements:
#' \itemize{
#'   \item A data.frame with SC values for all scenarios
#'   \item A data.frame with SSMD values
#'   \item A data.frame with p-values for SSMD values
#'   \item A data.frame with the scenario ranks based on SC and one based on SSMD
#'   \item The output of the Kendall's test
#' }
#' If \code{group_mean=TRUE} there will be six additional elements:
#' \itemize{
#'   \item A data.frame with the mean SC values for each parameter
#'   \item A data.frame with the mean SSMD values
#'   \item A data.frame with p-values calculated for the mean SSMD values
#'   \item A data.frame with the parameter ranks based on the mean SC and one
#'         based on the mean SSMD
#'   \item The output of the Kendall's test performed on the ranking of the
#'         parameters
#' }
#' @references Conroy, S. D. S., and B. W. Brook. 2003. Demographic sensitivity and
#' persistence of the threatened white- and orange-bellied frogs of Western
#' Australia. Population Ecology 45:105-114.
#'
#' Drechsler, M., M. A. Burgman, and P. W. Menkhorst. 1998. Uncertainty in
#' population dynamics and its consequences for the management of the
#' orange-bellied parrot \emph{Neophema chrysogaster}. Biological Conservation
#' 84:269-281.
#'
#' Zhang, X. D. 2007. A pair of new statistical parameters for quality control
#' in RNA interference high-throughput screening assays. Genomics 89:552-561.
#'
#' @import data.table
#' @importFrom irr kendall
#' @export
#' @examples
#' # Using Pacioni et al. example data. See ?pac.clas for more details.
#' data(pac.clas)
#' pairw<-pairwise(data=pac.clas, project="Pacioni_et_al", scenario="ST_Classic",
#'                params=c("Nall", "Het"), yrs=c(60,120), ST=TRUE,
#'                type="Single-Factor",
#'                SVs=c("SV1", "SV2", "SV3", "SV4", "SV5", "SV6", "SV7"),
#'                save2disk=FALSE)
pairwise <-  function(data,
                      project,
                      scenario,
                      params=c("PExtinct", "Nextant", "Het", "Nalleles"),
                      yrs="max",
                      ST=FALSE,
                      type=NA,
                      group.mean=FALSE,
                      SVs=NA,
                      save2disk=TRUE,
                      dir_out="DataAnalysis/Pairwise") {

  # Function definitions
  SEname <- function(par) paste("SE.", par, ".", sep="")
  SDname <- function(parSD) paste("SD.", parSD, ".", sep="")
  naming.coef <- function(naming) paste("SC", "_", naming, yr, sep="")
  naming.ssmd <- function(naming.ssmd)
                      paste("SSMD", "_", naming.ssmd, yr, sep="")
  kend <- function(tab.ranks) {
      k <- irr::kendall(tab.ranks[, - c(1:2), with=FALSE], TRUE)
      return(k)
  }

  # Flag (false) columns where all entries are NA
  FColsAllNA <- function(lranks)
                        apply(lranks, 2, function(chk) !all(is.na(chk)))

  # Error handling
  suppressWarnings(if (!yrs == "max" & !is.numeric(yrs))
    stop("invalid value(s) for 'yrs' "))

  fname <- if (ST) {
    paste(project, "_", scenario, sep="")
  } else {
    project
  }

  # set yrs to max
  suppressWarnings(if (yrs == "max") {yrs <- max(data$Year)})

  # set group.mean if needed
  if (ST & type == "Single-Factor" & length(SVs) > 1) {group.mean <- TRUE}

  # Set up headings for params and SE and SD
  params <- make.names(params)
  SE <- sapply(params, SEname)
  if ("r.stoch" %in% params) {SE["r.stoch"] <- "SE.r."}

  SD <- sapply(params, SDname)
  if ("r.stoch" %in% params) {SD["r.stoch"] <- "SD.r."}

  # Create a dataframe for the base scenario.
  if (ST) {
    select.base <- data$scen.name == (paste(scenario, "(Base)", sep=""))
    stbase <- subset(data, select.base)
    scen.name.base <- paste(scenario, "(Base)", sep="")
  } else {
    select.base <- data$scen.name == scenario
    stbase <- data[select.base, ]
    scen.name.base <- scenario
  }

  # calculate sensitivity coefficient & SSMD
  coef <- numeric(0)
  scenario.name <- unique(data$scen.name)
  scenario.name <- scenario.name[!scenario.name == scen.name.base]
  pops.name <- unique(data$pop.name)
  for (yr in yrs) {
    base <- numeric(0)
    for (popName in pops.name) {
      substdat <- subset(data, Year == yr & pop.name == popName &
                           !scen.name == scen.name.base)
      substbase <- subset(stbase, Year == yr & pop.name == popName)
      for (param in params) {
        base <- substbase[, param]
        SDbase <- substbase[, SD[match(param, params)]]
        coef.calc <- if (param == "PExtinct" | param == "PExtant") {
                     function(v) ( - (log(v/(1 - v)) - log(base / (1 - base))))
                     } else {
                     function(v) (v - base) / base
                     }
        coef <- sapply(substdat[, param], coef.calc )
        if (exists("tttable.coef")) {
          tttable.coef <- cbind.data.frame(tttable.coef, coef)
        } else {
          tttable.coef <- cbind.data.frame(scenario.name,
                                           substdat$pop.name, coef)
        }
        ssmd.numerator <- if (param=="PExtinct" | param=="PExtant") {
          (base - substdat[, param])
        } else {
          (substdat[, param] - base)
        }
        ssmd.denominator <- sqrt((substdat[, SD[match(param, params)]])^2 +
                                   (SDbase^2))
        ssmd <- ssmd.numerator / ssmd.denominator
        if (exists("ttssmd.table")) {
          ttssmd.table <- cbind.data.frame(ttssmd.table, ssmd)
        } else {
          ttssmd.table <- cbind.data.frame(scenario.name,
                                           substdat$pop.name, ssmd)
        }
      }
      if (exists("ttable.coef")) {
        ttable.coef <- rbind(ttable.coef, tttable.coef)
      } else {
        ttable.coef <- tttable.coef
      }
      if (exists("tssmd.table")) {
        tssmd.table <- rbind(tssmd.table, ttssmd.table)
      } else {
        tssmd.table <- ttssmd.table
      }
      rm(tttable.coef)
      rm(ttssmd.table)
    }
    if (exists("h.table.coef")) {
      h.table.coef <- c(h.table.coef, sapply(params, naming.coef))
    } else {
      h.table.coef <- c("Scenario", "Population", sapply(params, naming.coef))
    }
    if (exists("h.ssmd.table")) {
      h.ssmd.table <- c(h.ssmd.table, sapply(params, naming.ssmd))
    } else {
      h.ssmd.table <- c("Scenario", "Population", sapply(params, naming.ssmd))
    }
    table.coef.col.n <- dim(ttable.coef)[2]
    ssmd.table.col.n <- dim(tssmd.table)[2]

    if (exists("table.coef")) {
      table.coef <- cbind(table.coef, ttable.coef[3:table.coef.col.n])
    } else {
      table.coef <- ttable.coef
    }
    if (exists("ssmd.table")) {
      ssmd.table <- cbind(ssmd.table, tssmd.table[3:ssmd.table.col.n])
    } else {
      ssmd.table <- tssmd.table
    }
    rm(ttable.coef)
    rm(tssmd.table)
  }

  colnames(table.coef) <- h.table.coef  #replace headings in table.coef
  colnames(ssmd.table) <- h.ssmd.table  #replace headings in table.ssmd

  # vector with colnames for params with prefix for coef
  h.coef <- h.table.coef[3:length(h.table.coef)]
  h.ssmd <- h.ssmd.table[3:length(h.ssmd.table)]

  ssmd.table.pvalues <- if (is.null(dim(ssmd.table[, h.ssmd]))) {
                            sapply(ssmd.table[, h.ssmd], pval)
                            } else {
                              apply(ssmd.table[, h.ssmd], 2, pval)
                            }
  ssmd.table.pvalues <- cbind(ssmd.table[, 1:2], ssmd.table.pvalues)

  # Rank scenarios for each params by population based on sensitivity
  # coefficients and SSMD sensitivity coefficients
  DT.table.coef <- data.table(table.coef)
  setkey(DT.table.coef, Population, Scenario)
  ranks.sc <- DT.table.coef[, lapply(- abs(.SD), rank, na.last="keep"),
                            by=Population, .SDcols=h.coef] # Rank
  ranks.sc[, Scenario := DT.table.coef[,Scenario]]
  setcolorder(ranks.sc, c("Population", "Scenario", h.coef))

  # split by pop in a list where each element is a pop
  lranks.sc.pops <- split(ranks.sc, ranks.sc[, Population])
  # list of vectors with logical of cols to retain
  lsel <- lapply(lranks.sc.pops, FColsAllNA)
  ranks.sc.fin <- list()
  for (i in 1:length(lsel)) {
    # Remove flagged cols of NA from each element
    # TODO : Generate popName with Lapply and remove if (exists())
    ranks.sc.fin[[i]] <- lranks.sc.pops[[i]][, lsel[[i]], with=FALSE]
    popName <- as.character(ranks.sc.fin[[i]][1, Population])
    if (exists("popNames")) {
      popNames <- c(popNames, popName)
    } else {
      popNames <- popName
    }
  }
  names(ranks.sc.fin) <- popNames

  # ssmd
  DT.ssmd <- data.table(ssmd.table)
  setkey(DT.ssmd, Population, Scenario)
  ranks.ssmd <- DT.ssmd[, lapply(- abs(.SD), rank, na.last="keep"),
                        by=Population, .SDcols=h.ssmd] # Rank
  ranks.ssmd[, Scenario := DT.ssmd[, Scenario]]
  setcolorder(ranks.ssmd, c("Population", "Scenario", h.ssmd))

  # split by pop in a list where each element is a pop
  lranks.ssmd.pops <- split(ranks.ssmd, ranks.ssmd[, Population])
  # list of vectors with logical of cols to retain
  lsel <- lapply(lranks.ssmd.pops, FColsAllNA)
  ranks.ssmd.fin <- list()
  rm(popNames)
  for (i in 1:length (lsel)) {
    # Remove flagged cols of NA from each element
    ranks.ssmd.fin[[i]] <- lranks.ssmd.pops[[i]][, lsel[[i]], with=FALSE]
    # TODO : Generate popName with Lapply and remove if (exists())
    popName <- as.character(ranks.ssmd.fin[[i]][1, Population])
    if (exists("popNames")) {
      popNames <- c(popNames, popName)
    } else {
      popNames <- popName
    }
  }
  names(ranks.ssmd.fin) <- popNames
  kendall.out <- list(SC=NULL, SSMD=NULL)
  # NOTE : kendall function handles na listwise
  kendall.out$SC <- capture.output(
                      cat("Rank comparison of sensitivity coefficients", "\n"),
                      lapply(ranks.sc.fin, kend))
  kendall.out$SSMD <- capture.output(
                      cat("Rank comparison of SSMD", "\n"),
                      lapply(ranks.ssmd.fin, kend))
  if (save2disk) {
    # write results
    df2disk(table.coef, dir_out, fname, ".coef.table")
    df2disk(ssmd.table, dir_out, fname, ".SSMD.table")
    df2disk(ssmd.table.pvalues, dir_out, fname, ".SSMD.table.pvalues")
    df2disk(ranks.sc, dir_out, fname, ".ranks.sc")
    df2disk(ranks.ssmd, dir_out, fname, ".ranks.SSMD")
    capture.output(print(kendall.out, quote=FALSE),
                   file=paste0(dir_out, "/", fname,".kendall.txt"))
  }
  # Collate results in a list
  r.OneWay<-list(coef.table=table.coef,
                 SSMD.table=ssmd.table,
                 SSMD.table.pvalues=ssmd.table.pvalues,
                 ranks.SC=ranks.sc,
                 ranks.SSMD=ranks.ssmd,
                 Kendall=kendall.out)
  # if group.mean == TRUE calculate the mean for each SV and rank SVs
  if (group.mean) {

    # Calculate mean sensitivity coefficients and mean ssdm
    subpopsstdat <- subset(data, (Year == 0 & !scen.name == scen.name.base))
    # vector with colnames for params with prefix for coef
    h.coef <- h.table.coef[3:length(h.table.coef)]
    h.ssmd <- h.ssmd.table[3:length(h.ssmd.table)]
    for (SV in SVs) {
      SVbase <- stbase[stbase$Year == 0, SV]
      select.4.SV <- !subpopsstdat[, SV] == SVbase

      # mean coef calculations
      DT.coef <- data.table(table.coef)
      setkey(DT.coef, Scenario) # sort to match order in scen.4.SV
      DT.coef[,scen.4.SV := select.4.SV] # add scen.4.SV
      setkeyv(DT.coef, c("scen.4.SV", "Population"))
      # calculate mean by pop for that SV
      DT.coef <- DT.coef[J(TRUE), lapply(.SD, mean),
                         by=Population, .SDcols=h.coef]
      DT.coef[, SV := SV]

      # mean ssmd calculations
      DT.ssmd <- data.table(ssmd.table)
      setkey(DT.ssmd, Scenario)
      DT.ssmd[, scen.4.SV := select.4.SV]
      setkey(DT.ssmd, scen.4.SV)
      DT.ssmd <- DT.ssmd[J(TRUE), lapply(.SD, mean),
                         by=Population, .SDcols=h.ssmd]
      DT.ssmd[, SV:=SV]

      # appends mean calculations for both statistics & for each SV to table
      if (exists("mean.coef.table")) {
        mean.coef.table <- rbind(mean.coef.table, DT.coef)
      } else {
        mean.coef.table <- DT.coef
      }
      if (exists("mean.ssmd.table")) {
        mean.ssmd.table <- rbind(mean.ssmd.table, DT.ssmd)
      } else {
        mean.ssmd.table <- DT.ssmd
      }
    }

    setcolorder(mean.coef.table, c("Population", "SV", h.coef))
    setkeyv(mean.coef.table, c("Population", "SV"))
    setcolorder(mean.ssmd.table, c("Population", "SV", h.ssmd))
    setkeyv(mean.ssmd.table, c("Population", "SV"))

    mean.ssmd.table.pvalues <- mean.ssmd.table[,
                lapply(.SD, pval), .SDcols=h.ssmd]
    mean.ssmd.table.pvalues <- cbind(mean.ssmd.table[, list(Population, SV)],
                mean.ssmd.table.pvalues)

    # Rank SVs for each params by population based on mean sensitivity
    # coefficients and mean SSMD
    ranks.msc <- mean.coef.table[, lapply(-abs(.SD), rank, na.last="keep"),
                                 by=Population, .SDcols=h.coef] # Rank
    ranks.msc[, SV := SVs] #Add SV col
    setcolorder(ranks.msc, c("Population", "SV", h.coef))
    lranks.msc.pops <- split(ranks.msc,ranks.msc[, Population])
    lsel <- lapply(lranks.msc.pops, FColsAllNA)
    rm(popNames)
    ranks.msc.fin <- list()
    for (i in 1:length (lsel)) {
      ranks.msc.fin[[i]] <- lranks.msc.pops[[i]][,lsel[[i]], with=FALSE]
      popName <- as.character(ranks.msc.fin[[i]][1, Population])
      if (exists("popNames")) {
        popNames <- c(popNames, popName)
      } else {
        popNames <- popName
      }
    }
    names(ranks.msc.fin) <- popNames

    ranks.mssmd <- mean.ssmd.table[, lapply(-abs(.SD), rank, na.last="keep"),
                                   by=Population, .SDcols=h.ssmd] # Rank
    ranks.mssmd[, SV := SVs] #Add SV col
    setcolorder(ranks.mssmd, c("Population", "SV", h.ssmd))
    lranks.mssmd.pops <- split(ranks.mssmd,ranks.mssmd[, Population])
    lsel <- lapply(lranks.mssmd.pops, FColsAllNA)
    ranks.mssmd.fin <- list()
    rm(popNames)
    for (i in 1:length (lsel)) {
      ranks.mssmd.fin[[i]] <- lranks.mssmd.pops[[i]][, lsel[[i]], with=FALSE]
      popName <- as.character(ranks.mssmd.fin[[i]][1, Population])
      if (exists("popNames")) {
        popNames <- c(popNames, popName)
      } else {
        popNames <- popName
      }
    }
    names(ranks.mssmd.fin) <- popNames
    kendall.mean.out <- list(SC=NULL, SSMD=NULL)
    kendall.mean.out$SC <- capture.output(
                      print("Rank comparison of mean sensitivity coefficients"),
                      lapply(ranks.msc.fin, kend))
    kendall.mean.out$SSMD <- capture.output(
                      print("Rank comparison of mean SSMD"),
                      lapply(ranks.mssmd.fin, kend))

    if (save2disk) {
      df2disk(mean.coef.table, dir_out, fname, ".mean.coef.table")
      df2disk(mean.ssmd.table, dir_out, fname, ".mean.SSMD.table")
      df2disk(mean.ssmd.table.pvalues, dir_out, fname, ".mean.SSMD.table.pvalues")
      df2disk(ranks.msc, dir_out, fname, ".ranks.mSC")
      df2disk(ranks.mssmd, dir_out, fname, ".ranks.mSSMD")
      capture.output(print(kendall.mean.out, quote=FALSE),
                     file=paste0(dir_out, "/", fname, ".Kendall.means.txt"))
    }

    # Collate results for means
    r.OneWay$mean.coef.table <- mean.coef.table
    r.OneWay$mean.SSMD.table <- mean.ssmd.table
    r.OneWay$mean.SSMD.table.pvalues <- mean.ssmd.table.pvalues
    r.OneWay$ranks.mean.SC <- ranks.msc
    r.OneWay$ranks.mean.SSMD <- ranks.mssmd
    r.OneWay$Kendall.means <- kendall.mean.out
  }
  return(r.OneWay)
}


#' Search for the best regression model(s)
#'
#'
#' \code{fit_regression} fits either a Generalized Linear Model or a betareg model
#' to the data and search for the best model(s) given a list of predictors using
#' the R package glmulti.
#'
#' \code{fit_regression} fits a different type of regression model depending on the
#' dependent variable. When this is a count (e.g. N or the number of alleles),
#' the function will fit a Generalized Linear Model. The first fit is attempted
#' with a Poisson error distribution and if the dispersion parameter (calculated
#' as \eqn{residual deviance / df} is larger than (the somewhat
#' arbitrary cut off of) 1.5, the model will be refitted with a quasipoisson
#' error distribution (a message is displayed if this happens).
#'
#' \code{fit_regression} establishes whether the dependent variable is a count
#' by searching for it in \code{count_data}. If the users generated their own
#' dependent variable (e.g. through a PS), this has to be included in \code{count_data}
#' to indicate \code{fit_regression} that it is analysing count data.
#'
#' When the number of alleles is the dependent variable (from .run files), this
#' is rounded to integer (to meet R requirement that count data are integers)
#' before a GLM is fitted to the data.
#'
#' If \code{param} is a proportion (e.g. Gene Diversity and Inbreeding), then
#' the function uses a Beta regression from the R package betareg (Cribari-Neto
#' & Zeileis 2010). Different link functions are tested and the one with the
#' lowest AIC value is selected. The selected link function is displayed on the
#' R console and the difference in the AIC scores relative to the best link
#' function is also displayed.
#'
#' In the initial fit of the model the main and interactions effects are included.
#'
#' Successively, a search for the best model is carried out. This is performed
#' with the R package \code{\link[glmulti]{glmulti}} (Calcagno & de Mazancourt 2010).
#' \code{fit_regression} will conduct an exhaustive search if ncand is less or
#' equal to the number of candidate models, otherwise it will use a genetic
#' search method (see glmulti documentations for more details about the search
#' methods). When \code{\link[glmulti]{glmulti}} uses the genetic search method,
#' two small files (with extension \code{.modgen.back} and \code{.mods.back.} are
#' written in the working directory even if \code{save2disk=FALSE}.
#'
#' \code{fit_regression} explicitly ignores NA.
#'
#' Depending on the data, fitting several Beta regression models to complete the
#' search may be a long (and memory hungry) process. Also, the package betareg
#' has the limitation (at least at the moment of writing) that cannot handle
#' analysis of data when the dependent variable takes value of either exactly 0
#' or 1.
#'
#' See vignette for a more detailed explanation of \code{fit_regression}.
#'
#' @param data The long format of census (from \code{conv_l_yr}) or run (lrun,
#' the second element) of the output from \code{collate_run}
#' @param lookup (Optional) A look-up table where the scenario names are listed
#' together with the (missing) variables needed to fit the regression models
#' @param census Whether the input is census data
#' @param yr The year that has to be used in the analysis if census=TRUE
#' @param scenario Vortex scenario name
#' @param popn The sequential number of the population (in integer)
#' @param param The dependent variable
#' @param vs Character vector with independent variable(s)
#' @param count_data Character vector with param(s) that are counts and would use
#' a Poisson error distribution
#' @param ic Information criterion
#' @param l Level for glmulti search: 1 main effects, 2 main effects + interactions
#' @param ncand The threshold of candidate models after which switch to the
#' genetic search method, default: 30
#' @param set_size Value to be used in confsetsize (from \code{\link[glmulti]{glmulti}}
#'  The number of models to be looked for, i.e. the size of the returned confidence
#' set.)
#' @param dir_out The local path to store the output.
#'  Default: DataAnalysis/FitRegression
#' @inheritParams pairwise
#' @return A \code{glmulti} object with the best models found.
#' @references Calcagno, V., and C. de Mazancourt. 2010. glmulti: an R package
#' for easy automated model selection with (generalized) linear models. Journal
#' of Statistical Software 34:1-29.
#'
#' Cribari-Neto, F., and Zeileis, A. (2010) Beta regression in R. Journal of
#' Statistical Software 34(2).
#' @import data.table
#' @importFrom betareg betareg.fit
#' @importFrom glmulti glmulti
#' @export
#' @examples
#' # Using Pacioni et al. example data. See ?pac.run.lhs and ?pac.lhs for more
#' # details.
#' data(pac.run.lhs, pac.lhs)
#'
#' # Remove base scenario from .stdat data
#' pac.lhs.no.base <- pac.lhs[!pac.lhs$scen.name == "ST_LHS(Base)", ]
#'
#' # Use function lookup_table to obtain correct parameter values at year 0
#' lkup.ST_LHS <- lookup_table(data=pac.lhs.no.base, project="Pacioni_et_al",
#'                             scenario="ST_LHS",
#'                             pop="Population 1",
#'                             SVs=c("SV1", "SV2", "SV3", "SV4", "SV5", "SV6", "SV7"),
#'                             save2disk=FALSE)
#'
#' # Remove base scenario from .run output in long format
#' lrun.ST_LHS.no.base <- pac.run.lhs[[2]][!pac.run.lhs[[2]]$Scenario == "ST_LHS(Base)", ]
#'
#' reg <- fit_regression(data=lrun.ST_LHS.no.base, lookup=lkup.ST_LHS,
#'                       census=FALSE,
#'                       project="Pacioni_et_al", scenario="ST_LHS", popn=1,
#'                       param="N", vs=c("SV1", "SV2", "SV3"), l=2,  ncand=30,
#'                       save2disk=FALSE)
#'
#'  # Clean up of residual files written by glmulti
#'  # Note, in some OS (W) these files may be locked because in use by R and have
#'  # to be manually after the R session has been either terminated or restarted
#'  file.remove(c("Pacioni_et_al_ST_LHS_N.modgen.back",
#'                "Pacioni_et_al_ST_LHS_N.mods.back"))
#'
#' # Example of information you can obtained once you have run fit_regression
#'
#' # The formula for the best model
#' bestmodel <- reg@@formulas[1]
#'
#' # The formulae for the best 30 model
#' bestmodels <- reg@@formulas
#'
#' # List of IC values
#' qaicvalues <- reg@@crits
#'
#' # QAIC differences between the first 5 best models (stored in 'delta')
#' delta <- as.vector(NULL)
#' for (i in 1:5) {
#'   del <- qaicvalues[i+1] - qaicvalues[i]
#'   delta <- c(delta, del)
#' }
#'
#' # The best model's coefficients
#' coef.best <- coef(reg@@objects[[1]])
#'
#' # The model averaged coefficients
#' coef.all <- glmulti::coef.glmulti(reg)
#' coefs <- data.frame(Estimate=coef.all[,1],
#'                     Lower=coef.all[,1] - coef.all[,5],
#'                     Upper=coef.all[,1] + coef.all[,5])
#'
#' # Plot IC profile
#' plot(reg, type="p")
#'
#' # Plot of model averaged importance of terms
#' plot(reg, type="s")
fit_regression <-  function(data,
                            lookup=NA,
                            census=TRUE,
                            yr,
                            project,
                            scenario,
                            popn,
                            param="N" ,
                            vs=c("GS1"),
                            count_data=c("Nextant", "Nall", "Nalleles", "N", "AM", "AF", "Subadults",
                                         "Juv", "nDams", "nBroods", "nProgeny", "nImmigrants",
                                         "nEmigrants", "nHarvested", "nSupplemented", "YrExt", "Alleles"),
                            ic="aic",
                            l=1,
                            ncand=30,
                            set_size=NA,
                            save2disk=TRUE,
                            dir_out="DataAnalysis/FitRegression") {
  # Function definitions

  # Selecet method for glmulti search
  select_method <- function(cand, ncand) {
    if (cand > ncand) {
      m <- "g"
    } else {
      m <- "h"
    }
    return(m)
  }

  #  Updated the betareg model with the link x and extract AIC
  LinkTest <- function(x) logLik(update(breg, link=x))

  # vector with available link functions to be used with betareg
  links <- c("logit", "probit", "cloglog", "cauchit", "loglog")

  suppressWarnings(if (!is.na(lookup)) {
    data <- plyr::join(data, lookup, by='Scenario', type="left")
  } )

  # convert data.table
  data <- data.table(data)

  # select yr if census and create pop vector
  if (census) {
    # Select the yr
    setkey(data, Year)
    data <- data[J(yr), ]
    # Select the pop (with its number)
    pop <- paste0("pop", popn)
  } else {
    pop <- levels(data$Population)[popn]
  }

  # select pops
  setkey(data, Population)
  data <- data[pop]

  xs <- paste(vs, collapse="*") # get ready the vs for the formula

  if (param %in% count_data) {
    data[ , (param) := round(.SD), .SDcols=param]
    fam <- "poisson"
  }

  # Plot param
  # Tried DT[, plot(GeneDiv)] to call the plot directily within data.table
  # but when I use param rather than the column names DT return a datatable
  # rather than a vector (as it does with GeneDiv)
  paramvalues <- data[[param]]
  name <- paste(project, scenario, param, sep="_")
  if (save2disk) {
    dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)
    pdf(paste0(dir_out, "/", paste(name, "histogram.pdf", sep="_")))
  }

  hist(paramvalues,  main=paste("Histogram of", param), xlab=param)
  if (save2disk) {dev.off()}

  message(paste("summary of", param))
  print(summary(paramvalues))
  # set up formula to be used in regression models
  formula <- as.formula(paste0(param, "~", xs))

  if (param %in% count_data) {
    # Preliminary fit the GLM
    message("Fitting a GLM...")
    glm1 <- glm(data=data, formula, family=fam, na.action=na.omit)

    # Check number of candidate models.
    cand <- do.call("glmulti", list(glm1, method="d", family=fam,
                                    level=l, na.action=na.omit))

    # Calculates ratio Res deviance to df
    chat <- deviance(glm1) / df.residual(glm1)

    # if the ratio (chat) is >1.5 then set the overdispersion parameter
    if (chat > 1.5) {
      message("Overdispersion was detected in the data.")
      # Fit the GLM with quasi- error distribution to get the dispersion parameter
      glm2 <- glm(data=data, formula, family="quasipoisson", na.action=na.omit)

      # Changed the c value for qaic search with glmulti
      R.utils::setOption("glmulti-cvalue", summary(glm2)$dispersion)
      message(paste("Setting overdispersion parameter to:",
                    getOption("glmulti-cvalue")))
      message("NOTE: the Information Criterion for model search was changed to QAIC")
      ic <- "qaic"
    }
    m <- select_method (cand, ncand)
    if (is.na(set_size)) {
      set_size <- min(cand, ncand)
      message(paste("confsetsize set to", set_size))
    }
    message(paste("Search method set to", m))
    message(paste("Search for best candidate models using level =",
                  l, "started..."))

    # Search for best model(s) with glmulti
    tnp <- system.time(best.mod <- do.call("glmulti",
                                           list(glm1, family=fam, crit=ic,
                                               method=m, confsetsize=set_size,
                                               plotty=FALSE, report=FALSE, level=l,
                                               name=name, na.action=na.omit)))
  } else {
    message("Fitting a beta regression...")
    message("Searching for the best link function...")

    # Preliminary fit of betw regression
    breg <- betareg(formula, data=data, na.action=na.omit)

    # Generate vector with AIC values for the available link functions
    linkAIC <- sapply(links, LinkTest)
    linkpos <- which.min(linkAIC) # Select lower AIC (position)
    message(paste("Best link function:", links[linkpos]))
    message(paste("Link function AIC differences relative to", links[linkpos]))
    print(deltalinkAIC <- linkAIC - linkAIC[linkpos])

    # Check number of candidate models.
    cand <- do.call("glmulti", list(formula, data=data, method="d", level=l,
                                    fitfunc=betareg, link=links[linkpos],
                                    na.action=na.omit))

    # This is repeated but it can't go later because I want this to be displayed
    # before the search starts.
    # The search may be very long and if something is wrong the user have a
    # chance to stop it.
    # It can't go earlier because the cand call is different depending
    # on the model used
    m <- select_method (cand, ncand)
    if (is.na(set_size)) {
      set_size <- min(cand, ncand)
      message(paste("confsetsize set to", set_size))
    }
    message(paste("Using search method:", m))
    message(paste("Search for best candidate models using level =", l,
                  "started..."))
    tnp <- system.time(best.mod <- do.call("glmulti",
                                           list(formula, data=data,  crit=ic,
                                                method=m, confsetsize=set_size,
                                                plotty=FALSE, report=FALSE, level=l,
                                                name=name, fitfunc=betareg,
                                                link=links[linkpos],
                                                na.action=na.omit)))
  }
  message("Done! Elapsed time:")
  print (tnp)

  if (save2disk) {
    message("Best models saved to disk in the file ...best.mod.rda")
    save(best.mod, file=paste0(dir_out, "/", name, "_best.mod.rda"))
    pdf(paste0(dir_out, "/", name, "_IC_plot.pdf"))
    plot(best.mod, type="p")
    dev.off()
  }
  return(best.mod)
}
