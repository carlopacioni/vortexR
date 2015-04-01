#------------------------------------------------------------------------------#
# Data
#------------------------------------------------------------------------------#

#' @name sta
#' @title Campbell example, project "Starlingv3PopBased" (.dat)
#' @description Data from Campbell example, project "Starlingv3PopBased" (.dat)
#' @format a \code{data.frame} with 1632 obervations of 44 variables
#' @source Campbell study
NULL

#' @name sta.evy5
#' @title Campbell example, project "Starlingv3PopBased" (.dat),
#'   sensitivity test scenario "MReductEvy5" (.stdat)
#' @description Data from Campbell example, project "Starlingv3PopBased" (.dat),
#'   sensitivity test scenario "MReductEvy5" (.stdat)
#' @format a \code{data.frame} with 1020 obervations of 47 variables
#' @source Campbell study
NULL

#' @name sta.b11
#' @title Campbell example, project "Starlingv3PopBased" (.dat),
#'   sensitivity test scenario "MReduction_B11_09Evy5" (.stdat)
#' @description Data from Campbell example,
#'   project "Starlingv3PopBased" (.dat),
#'   sensitivity test scenario "MReduction_B11_09Evy5" (.stdat)
#' @format a \code{data.frame} with 1020 obervations of 47 variables
#' @source Campbell study
NULL

#' @name pac.clas
#' @title Pacioni et al. example, project "Pacioni_et_al",
#'   sensitivity test scenario "ST_Classic" (.stdat)
#' @description Data from Pacioni et al. example,
#'   project "Pacioni_et_al",
#'   sensitivity test scenario "ST_Classic" (.stdat)
#' @format a \code{data.frame} of 2904 observations of 68 variables
#' @source Pacioni et al.
NULL

#' @name pac.lhs
#' @title Pacioni et al. example, project "Pacioni_et_al",
#' sensitivity test scenario "ST_LHS" (.stdat)
#' @description Data from Pacioni et al. example, project "Pacioni_et_al",
#'   sensitivity test scenario "ST_LHS" (.stdat)
#' @format A \code{data.frame} of 6171 observations of 68 variables
#' @source Pacioni et al.
NULL

#' @name pac.run.lhs
#' @title Pacioni et al. example, project "Pacioni_et_al",
#' sensitivity test scenario "ST_LHS" (.stdat), "run"
#' @description Data from Pacioni et al. example, project "Pacioni_et_al",
#'   sensitivity test scenario "ST_LHS" (.stdat), "run"
#' @format A named list of two \code{data.frame}s:
#'   run (153 obs, 7 var), lrun (153 obs, 8 var)
#' @source Pacioni et al.
NULL

#' @name pac.yr
#' @title Pacioni et al. example, project "Pacioni_et_al",
#' sensitivity test scenario "ST_Classic" (.yr), "yr"
#' @description Data from Pacioni et al. example, project "Pacioni_et_al",
#'   sensitivity test scenario "ST_Classic" (.yr), "yr"
#' @format A named list of two \code{element}s:
#'   all (8712 obs, 26 var), means (2904 obs, 25 var)
#' @source Pacioni et al.
NULL

#------------------------------------------------------------------------------#
# Exported helper functions
#------------------------------------------------------------------------------#

#' Save a data.frame as both Rdata and CSV
#'
#' \code{df2disk} saves a given data.frame as both Rdata and CSV with a given
#' name and optional name postfix to a given location.
#'
#' \code{df2disk} is used by the \code{collate_} functions when the operator
#' chooses to save2disk.
#'
#' @param df A data.frame
#' @param dirpath The destination path for written files, will be created if
#'   necessary
#' @param fname The file name
#' @param postfix An optional name postfix
#' @usage
#' my.df <- data.frame(1, 1:10, sample(LETTERS[1:3], 10, replace = TRUE))
#' my.folder <- paste0(getwd(), "/test")
#' df2disk(df=my.df, dirpath=getwd(), fname="testname")
#' df2disk(df=my.df, dirpath=my.folder, fname="testname", postfix="_testpostfix")
#' @export
df2disk <- function(df, dirpath, fname, postfix=""){

  dir.create(dirpath, showWarnings=FALSE)

  save(df, file=paste0(dirpath, "/", fname, postfix, ".rda"), compress="xz")

  write.table(df, file=paste0(dirpath, "/", fname, postfix, ".txt"),
              sep=";", row.names=FALSE)
}


#------------------------------------------------------------------------------#
# Internal helper functions
#------------------------------------------------------------------------------#

#' Return file paths of files matching a pattern in a directory
#'
#' @param path The directory to search in
#' @param pattern A pattern to match file names
#' @param fn_name The vortexR function name for verbose messages
#' @param verbose Progress messages, default: FALSE
#' @return A character vector of fully qualified file paths
#' @import gtools
get_file_paths <- function(path, pattern, fn_name, verbose=FALSE){
  files <- gtools::mixedsort(list.files(path=path,
                                        pattern=pattern,
                                        full.names=TRUE))

  if (length(files) == 0) {
    stop(paste0("ERROR vortexr::", fn_name," found no files",
                  " containing '", fname, "' in ", path))
  } else {
    if (verbose){
      msg <- paste0("INFO vortexR::", fn_name," found ", length(files),
                   " matching files in ", path, ":")
      message(msg)
    }
    return(files)
  }
}


#' Standard Error from a vector
#'
#' @param se A standard error of a set of values
#' @param no The number of values
#' @return The standard deviation of the values
se2sd <- function(se, no) {se * sqrt(no)}


#' Return a prefixed and repeated string of characters.
#'
#' @param chars A string of characters (popvalue)
#' @param times The number of repetitions (ncolpop), default: 1
#' @param prefix A text prefix, default: ""
#' @usage
#' PrefixAndRepeat(c("a","b","c"), 3, "pop_")
PrefixAndRepeat <- function(chars, times=1, prefix="") {
  rep(paste0(prefix, chars), times)
}


#' Remove subheadings and downfill scenario name and subheadings as columns.
#'
#' @param iter The iteration (run) number
#' @param filename The fully qualified filename to read from
#' @param n_rows The number of rows to read from the file
#' @param iter_ln The number of rows to skip from the file
#' @param lines An object returned from readLines()
#' @param header A character vector of column names
#' @return A data.frame
CompileIter <- function (iter, filename, n_rows, iter_ln, lines, header) {
  temp.df <- read.table(filename, header=FALSE, sep=";",
                        nrows=n_rows, skip=iter_ln[iter],
                        colClasses="numeric", comment.char="")
  colnames(temp.df) <- header
  Iteration <- rep(iter, length=length(temp.df$Year))
  ScenNameStarts <- attr(regexpr(pattern="\\$Scenario: ", lines[1]),
                         "match.length") + 1
  ScenName <- substr(lines[1], ScenNameStarts, nchar(lines[1]))
  Scenario <- rep(ScenName, length=length(temp.df$Year))
  temp.df <- cbind(Scenario, Iteration, temp.df)
  return(temp.df)
}


#------------------------------------------------------------------------------#
# Data loading
#------------------------------------------------------------------------------#

#' Collate one local Vortex output file into a data.frame.
#'
#' \code{collate_one_dat} parses one Vortex .dat or .stdat file, and returns the
#' data within as one data.frame.
#'
#' @param filename The fully qualified filename of a Vortex .dat or .stdat file
#' @param runs The number of simulation runs
#' @param verbose Progress messages, default: FALSE
#' @return a data.frame with data from one .dat or .stdat file and
#'  population/scenario names as factors
#' @export
collate_one_dat <- function(filename, runs, verbose=FALSE){

  if (verbose) {message(cat("INFO vortexR::collate_one_dat parsing", filename))}
  lines <- readLines(filename)

  # Blocks of population data start with "Population"
  # If there are more than one population, a last block "Metapopulation" exists
  popLn <- grep(pattern="^Population", lines)
  popN <- length(popLn)

  if (popN > 1) {
    popLn <- c(popLn, grep(pattern="^Metapopulation", lines))
    popN <- popN + 1
    readFor <- popLn[2] - popLn[1] - 1
  } else {
    readFor <- length(lines) - popLn[1]
  }

  if (verbose) {
    message(cat("INFO vortexR::collate_one_dat found", popN, "populations"))
  }

  # Column names
  h <- make.names(read.table(filename, header=FALSE, sep=";", nrows=1,
                              skip=(popLn[1] - 3), colClasses="character"))

  # Scenario name
  scen <- read.table(filename, header=FALSE, sep=";", nrows=1,
                     skip=(popLn[1] - 2), colClasses="character")

  # Loop over populations
  for (pop in 1:popN) {

    readAfter <- popLn[pop]

    # Read population data block
    tmp <- read.table(filename, header=FALSE, sep=";",
      nrows=readFor, skip=readAfter,
      colClasses="numeric", comment.char="")
    colnames(tmp) <- h

    # Create vector with population name
    popNameStarts <- attr(regexpr(pattern="^Population ", lines[popLn]),
                          "match.length") + 4
    popName <- substr(lines[popLn], popNameStarts, nchar(lines[popLn]))
    if (pop > 1 & pop == max(popN)) {
      pop.name <- rep("Metapopulation", length=length(tmp$Year))
    } else {
      pop.name <- rep(popName[pop], length=length(tmp$Year))
    }

    # Create vector with scenario name
    scen.name <- unlist(scen.name <- rep(scen, length=length(tmp$Year)))

    # Create vectors of Standard Deviations of PExtant and PExtinct
    SD.PExtant. <- sapply(tmp$SE.PExtant., se2sd, no=runs)
    SD.PExtinct. <- sapply(tmp$SE.PExtinct., se2sd, no=runs)

    # Add scenario name, popname, and probability SDs as columns to data.frame
    tmp <- cbind(scen.name, pop.name, tmp, SD.PExtant., SD.PExtinct.)
    if (pop == 1) {x <- tmp} else {x <- rbind(x, tmp)}
  }
  return(x)
}

#' Collate Vortex .dat or .stdat output files into one data.frame.
#'
#' \code{collate_dat} collates all Vortex output files in a given directory
#' matching a given project name (and scenario name when relevant)into one
#' data.frame using \code{collate_one_dat}.
#'
#' The number of Vortex simulation runs has to be specified, as it cannot be
#' inferred by vortexR.
#'
#' To read Vortex output files from Sensitivity Testing with the extension
#' ".stdat", specify a scenario name. It is assumed that only file names
#' from Sensitivity Testing contain a scenario name.
#'
#' The given dir.in directory may contain other files; only files matching the
#' project (and, optionally, the scenario) name will be read.
#' If no matching files are found in the given directory, NULL will be returned.
#'
#' @param project Vortex project name to be imported
#' @param runs The number of Vortex simulation runs
#' @param scenario The scenario name if ST, default: NULL
#' @param dir.in The local folder containing Vortex files, default: NULL. If
#'   not specified, will fall back to use current working directory.
#' @param save2disk Whether to save the data as rda and csv, default: TRUE
#' @param dir.out The local path to store the output. Default: ProcessedData
#' @param verbose Progress messages, default: TRUE
#' @return a data.frame with data from all matching Vortex files or NULL
#' @import gtools
#' @export
#' @examples
#' # Using Campbell and Pacioni et al. example files
#' camp.dir <- system.file("extdata", "campbell", package="vortexr")
#' pac.dir <- system.file("extdata", "pacioni", package="vortexr")
#'
#' # Campbell example, project "Starlingv3PopBased" (.dat)
#' sta <- collate_dat("Starlingv3PopBased", 1000, dir.in = camp.dir)
#'
#' # Campbell example, project "Starlingv3PopBased",
#' #   sensitivity test scenario "MReductEvy5" (.stdat)
#' sta.evy5 <- collate_dat("Starlingv3PopBased", 1000,
#'   scenario = "MReductEvy5", dir.in = camp.dir)
#'
#' # Campbell example, project "Starlingv3PopBased",
#' #   sensitivity test scenario "MReduction_B11_09Evy5" (.stdat)
#' sta.b11 <- collate_dat("Starlingv3PopBased", 1000,
#'   scenario = "MReduction_B11_09Evy5", dir.in = camp.dir)
#'
#' # Pacioni et al. example, project "Pacioni_et_al",
#' #   sensitivity test scenario "ST_Classic" (.stdat)
#' pac.clas <- collate_dat("Pacioni_et_al", 1000,
#'   scenario = "ST_Classic", dir.in = pac.dir)
#'
#' # Pacioni et al. example, project "Pacioni_et_al",
#' #   sensitivity test scenario "ST_LHS" (.stdat)
#' pac.lhs <- collate_dat("Pacioni_et_al", 1000,
#'   scenario = "ST_LHS", dir.in = pac.dir)
#'
#' # Save collated data as Rda and CSV to ./out
#' \dontrun{
#' pac.lhs <- collate_dat("Pacioni_et_al", 1000,
#'   scenario = "ST_LHS", dir.in = pac.dir,
#'   save2disk = TRUE, dir.out = system.file(getwd(), "out"))
#' }
#'
#' # Dirpath falls back to current working directory
#' setwd(camp.dir)
#' sta <- collate_dat("Starlingv3PopBased", 1000)
collate_dat <- function(project=NULL, runs,
                        scenario=NULL,
                        dir.in=NULL,
                        save2disk=TRUE,
                        dir.out="ProcessedData",
                        verbose=TRUE){

  if (is.null(scenario)) {
    fname <- project
    pat <- paste0("^", fname, ".*\\.dat$")
  } else {
    fname <- paste(project, scenario, sep="_")
    pat <- paste0("^", fname, ".*\\.stdat$")
  }

  if (is.null(dir.in)) {dir.in <- getwd()}

  files <- get_file_paths(path=dir.in,
                          pattern=pat,
                          fn_name="collate_dat",
                          verbose=verbose)

  d <- data.frame()
  if (verbose){message("vortexR::collate_dat is reading:")}
    for (filename in files) {
      if (verbose){message(filename, "\r")}
      d <- rbind(d, collate_one_dat(filename, runs))
    }
    if (save2disk == T) {df2disk(d, dir.out, fname, "_data")}
    return(d)
}


#' Collate Vortex .run output files into one named list.
#'
#' @param npops The total number of simulated populations including the
#'   metapopulation
#' @inheritParams collate_dat
#' @return a list with two elements: run, a data.frame with data from all
#'   Vortex files and lrun, where the same data are re-arranged in long format
#' @export
#' @usage
#' pac.dir <- system.file("extdata", "pacioni", package="vortexr")
#' pac.run.lhs <- collate_run("Pacioni_et_al", "ST_LHS", 1, dir.in = pac.dir)
collate_run <- function(project=NULL,
                        scenario=NULL,
                        npops=1,
                        dir.in=NULL,
                        save2disk=TRUE,
                        dir.out="ProcessedData",
                        verbose=TRUE) {

  run <- data.frame()
  lrun <- data.frame()

  if (is.null(dir.in)) {dir.in <- getwd()}
  fname <- paste(project, scenario, sep="_")

  files <- get_file_paths(path=dir.in,
                          pattern=paste0("^", fname,".*\\.run$"),
                          fn_name="collate_run",
                          verbose=verbose)

  if (verbose){message("vortexR::collate_run is reading:")}
  for (filename in files) {
    if (verbose) {message(filename, "\r")}
    h <- gsub(" ", "", read.table(filename, header=FALSE, sep=";",
            nrows=1, skip=2, colClasses="character", comment.char=""))

    trun <- read.table(filename, header=FALSE, sep=";", skip=3,
                       colClasses="numeric", comment.char="")
    colnames(trun) <- h

    Scenario <- read.table(filename, header=FALSE, sep=":", nrows=1, skip=0,
                           colClasses="character", , comment.char="")
    Scenario <- gsub(" ", "", Scenario)[2]

    Scenario <- rep(Scenario, length=length(trun$Iteration))
    trun <- cbind(Scenario, trun)
    run <- rbind(run, trun)
  }

  # Number of cols for each pop except the intial, fixed cols
  ncolpop <- (length(h) - 1) / npops
  popnames <- read.table(filename, header=FALSE, sep=";", nrows=1, skip=1,
                         colClasses="character", comment.char="")
  popnames <- gsub(" ", "", popnames)
  pop <- unique(popnames[2:length(popnames)])

  # Reshape df in long format
  for (numPop in 1:npops) {
    k <- numPop - 1
    Cstart <- 3 + (k*ncolpop)
    Cend <- 3 + ncolpop * numPop - 1
    Population <- rep(pop[numPop], nrow(run))
    tlrun <- cbind(Population, run[ , 1:2], run[ , Cstart:Cend])
    lrun <- rbind(lrun, tlrun)
  }

  # paste pop names to headings
  h[2:length(h)] <- paste(h[2:length(h)], popnames[2:length(h)], sep="_")
  h <- c("Scenario", h)
  names(run) <- h # replace heading in run

  if (save2disk == T) {
    df2disk(run, dir.out, fname, "_run")
    df2disk(lrun, dir.out, fname, "_lrun")
  }

  r.RunST <- list(run=run, lrun=lrun)
  return(r.RunST)
}


#' Collate Vortex .yr output files into one named list.
#'
#' \code{collate_yr} collates all the .yr output from Vortex into one R object
#' and calculates the mean for each simulated year across all iterations
#'
#' @param npops_noMeta The total number of populations excluding the metapopulation,
#' default: 1
#' @inheritParams collate_dat
#' @return a list with two elements: "census", a data.frame with data from all
#' Vortex files and "census_means", a data.table with the mean of each parameter
#' across all iterations for each simulated year
#' @import data.table
#' @export
#' @usage
#' pac.dir <- system.file("extdata", "pacioni", package="vortexr")
#' pac.yr <- collate_yr("Pacioni_et_al", "ST_Classic", 1, dir.in = pac.dir)
collate_yr <- function(project=NULL,
                        scenario=NULL,
                        npops_noMeta=1,
                        dir.in=NULL,
                        save2disk=TRUE,
                        dir.out="ProcessedData",
                        verbose=TRUE) {

  if (is.null(dir.in)) {dir.in <- getwd()}
  fname <- paste0(project, "_", scenario)

  files <- get_file_paths(path=dir.in,
                          pattern=paste0("^", fname, ".*\\.yr$"),
                          fn_name="collate_yr",
                          verbose=verbose)

  censusData <- vector(mode="list", length=length(files))

  if (verbose){message("vortexR::collate_yr is reading:")}
  for (i in 1:length(files)) {
    lines <- readLines(files[i])
    if (verbose) {message(files[i], "\r")}

    # Header
    header <- as.vector(sapply(strsplit(lines[3], ";"), stringr::str_trim))
    ncolpop <- (length(header) - 1) / npops_noMeta
    hsuff <- as.vector(sapply(1:npops_noMeta, PrefixAndRepeat,
                              times=ncolpop, prefix="pop"))
    header <- c("Year", paste0(header[2:length(header)], hsuff))

    # Line numbers of "Iteration" subheadings
    iter_ln <- grep(pattern="^Iteration ", lines)

    # Number of rows underneath subheading
    n_rows <- if (length(iter_ln) > 1) (iter_ln[2] - iter_ln[1] - 1) else - 1

    # Extract data from each Iteration subheading in file
    one_yr <- lapply(1:length(iter_ln), CompileIter,
                     files[i], n_rows, iter_ln, lines, header)

    censusData[[i]] <- data.table::rbindlist(one_yr)
  }
  censusAll <- data.table::rbindlist(censusData)

#   # Functional implementation without for-loop ------------------------------#
#   # ll is a list of lists (one per file) of strings (one per line)
#   ll = lapply(files, function(x) readLines(x))
#
#   # Assuming the third line contains all columns names
#   header <- as.vector(sapply(strsplit(ll[[1]][[3]], ";"), stringr::str_trim))
#   no_col <- length(header)
#
#   # Wafflestomp the list of lists of character vectors into one data.frame
#   d <- plyr::ldply(ll, function(x){stringr::str_split_fixed(x, ";", no_col)})
#
#   setnames(d, header) # set variable names
#   # TODO downfill scenario name and iteration number, see
#   # http://stackoverflow.com/questions/10554741/r-fill-in-data-frame-with-values-from-rows-above
#   # TODO remove non-data lines (scen, pop, header, iter)
#   # hand over as censusAll
#   # end functional ----------------------------------------------------------#


  data.table::setkey(censusAll, Scenario, Year)
  censusMeansDT <- censusAll[ , lapply(.SD, mean), by="Scenario,Year"]
  censusMeansDT <- censusMeansDT[ , Iteration := NULL]

  if (save2disk == T) {
    df2disk(censusAll, dir.out, fname, "_census")
    df2disk(censusMeansDT, dir.out, fname, "_census_means")
  }

  return(list(census=censusAll, census_means=censusMeansDT))
}

#------------------------------------------------------------------------------#

#' Collate processed data generated by any of the 'collate' functions
#'
#' \code{collate_proc_data} collates multiple dfs generated by any of the
#' 'collate' functions. Only dfs generated by the same function can be comined
#' together, Missing data are filled with NA.
#'
#' @param data A list with a df from \code{collate_[dat, yr, run]} as each element
#' @inheritParams collate_dat
#' @return A data.frame with the collated data. Missing data are filled with NA.
#' @import plyr
#' @export
collate_proc_data <- function(data=NULL,
                              save2disk=TRUE,
                              dir.out="ProcessedData"){

  # TODO use data.table::rbindlist
  dfs <- plyr::rbind.fill(data)

  if (save2disk == T) {df2disk(dfs, dir.out, "CombinedDB")}

  return(dfs)
}

#' Convert 'census' (generated by collate_yr) into long format
#'
#' \code{conv_l_yr} converts the first element of the output from \code{collate_yr}
#' (census) in long format. This can then be fed into downstream analysis
#' (e.g. fit_regression)
#'
#' @param data The df 'census' from \code{collate_yr}
#' @param appendMeta Whether to calculate data for the metapopulation,
#'   default: FALSE
#' @param project Vortex project name (used to name the output)
#' @param scenario Vortex scenario name (used to name the output)
#' @param yrs The year(s) that need to be retained in the output
#' @inheritParams collate_yr
#' @return the census data.frame in long format
#' @import data.table plyr
#' @export
conv_l_yr <- function(data=NULL,
                      npops_noMeta=1,
                      appendMeta=FALSE,
                      project=NULL,
                      scenario=NULL,
                      yrs=c(1, 2),
                      save2disk=TRUE,
                      dir.out="ProcessedData") {
  requireNamespace("plyr", quietly=TRUE)
  requireNamespace("data.table", quietly=TRUE)

  LongFormat <- function(numPop) {
    k <- numPop - 1
    cstart <- 4 + (k * ncolpop)
    cend <- 4 + ncolpop * numPop - 1
    pop <- rep(paste0("pop", numPop), nrow(data))
    one.pop <- data[ , list(Scenario,Iteration,Year)]
    one.pop <- one.pop[ , "Population" := pop]
    one.pop[ , (h1) := data[ , cstart:cend, with=FALSE]]
    return(one.pop)
  }

  h <- names(data)
  # Number of cols for each pop except the first 3, which are fixed cols
  ncolpop <- (length(h) - 3) / npops_noMeta
  # Headings without pop suffix
  h1 <- gsub(pattern="pop1", "", h[4:(4 + ncolpop - 1)])
  setkey(data, Year)
  data<- data[.(yrs), ] # Select the yrs

  # Reshape df in long format
  tldata <- lapply(1:npops_noMeta, LongFormat)

  lcensus <- rbindlist(tldata)

  # If more than one pop, do Metapopulation calculations and rbind
  if (npops_noMeta > 1 & appendMeta == T) {
    message("Doing calculations for Metapopulation...")
    message("Please wait...")

    census <- c("N", "AM", "AF", "Subadults", "Juv",
                "nDams", "nBroods", "nProgeny")
    meta <-lcensus[, lapply(.SD, sum), by=list(Scenario,Iteration,Year),
                 .SDcols=census]
    message("Done...")
    message("Appending Metapopulation data to Census data frame...")
    gs <- names(lcensus)[grep(pattern="^GS", names(lcensus))]
    meta[ , Population := rep(paste0("pop", npops_noMeta + 1), nrow(data))]
    setcolorder(meta, c("Scenario", "Iteration", "Year", "Population", census))
    meta[ , (gs) := tldata[[1]][ , gs, with=FALSE]]

    # with data.table v1.9.3 can use
    # rbindlist(list(lcensus, meta), use.names=TRUE, fill=TRUE)
    lcensus <- rbind.fill(lcensus, meta)
    message("Done!")
  }

  if (save2disk == T) {
    df2disk(lcensus, dir.out, paste0(project, "_", scenario), "_lcensus")
  }
  return(lcensus)
}

#' Line plots of Vortex parameters vs years
#'
#' \code {line_plot_year} generates line plots of the selected Vortex parameters
#' for the selected populations, for all simulated years
#'
#' @param data A df from \code{collate_dat}
#' @param project Vortex project name (used to name the output)
#' @param scenario Vortex scenario name (used to name the output)
#' @param save2disk Whether to save the output to disk, default: TRUE
#' @param params Vortex parameters to be plotted,
#' default: c("PExtinct", "Nextant", "Het", "Nalleles")
#' @param plotpops The populations to be included in the plot, default: 'all'
#' @param save2disk Whether to save the output to disk, default: TRUE
#' @param dir.out The local path to store the output. Default: Plots
#' @return line plot(s)
#' @import ggplot2
#' @import grid
#' @export
line_plot_year <- function(data=NULL,
                           project=NULL,
                           scenario=NULL,
                           params=c("PExtinct", "Nextant", "Het", "Nalleles"),
                           plotpops=c("all"),
                           save2disk=TRUE,
                           dir.out="Plots") {

  require(ggplot2)
  require(grid)

  if (plotpops == "all")
    plotpops <- unique(data$pop.name)

  # Set up headings for params
  params <- make.names(params)

  # set legend size
  nscen <- length(unique(data$scen.name))
  if (nscen < 32){
    legKeySize <- 7
    legTxtSize <- 7
  } else {
    if (nscen > 60) {
      legKeySize <- 2.5
      legTxtSize <- 5
    } else {
      legKeySize <- round(160 / nscen, digits=2)
      legTxtSize <- 5
    }
  }

  fname_root <- if (is.null(scenario)) {
    project
  } else {
    paste0(project, "_", scenario)
  }

  popstdat <- subset(data, pop.name %in% plotpops)

  r.line_plot_year <- list()
  i <- 0
  if (save2disk == T) {
    dir.create(dir.out, showWarnings=FALSE)
    pdf(paste(dir.out, "/", fname_root, "_", "YearVsParams.pdf", sep=""))
  }

  for (param in params) {
    i <- i + 1
    root <- paste0(fname_root, "_", param)
    g <- ggplot(popstdat, aes_string(x="Year", y=param)) +
      geom_line(aes(color=scen.name)) +
      facet_grid(pop.name~.)
    plot <- g + theme(panel.margin=unit(0.2, "inches"),
                      legend.text=element_text(size=legTxtSize),
                      legend.key.size=unit(legKeySize, "mm")) +
      scale_colour_discrete(name="Scenarios")
    print(plot)
    r.line_plot_year[[i]] <- plot
    assign(paste(root, "_", "plot", sep=""), g)
  }

  if (save2disk == T) {
    dev.off()
    save(list=(ls(pattern=paste(fname_root, "_", ".*", "_", "plot", sep=""))),
         file=paste(dir.out, "/", fname_root, "_", "YearVsParams.rda", sep=""))
  }

  names(r.line_plot_year) <- ls(pattern=paste(fname_root, "_", ".*", "_",
                                            "plot", sep=""))
  return(r.line_plot_year)
}

#' Line plots of Vortex parameters vs years
#'
#' \code {line_plot_year_mid} generates line plots of the selected Vortex parameters
#' for the selected populations, from year zero to yrmid
#'
#' TODO working @examples
#'
#' @param yrmid The last year to plot
#' @inheritParams line_plot_year
#' @return line plot(s)
#' @import ggplot2
#' @import grid
#' @export
line_plot_year_mid <-  function(data=NULL,
                                project=NULL,
                                scenario=NULL,
                                yrmid=1,
                                params=c("PExtinct", "Nextant", "Het", "Nalleles"),
                                plotpops=c("all"),
                                save2disk=TRUE,
                                dir.out="Plots") {

  require(ggplot2)
  require(grid)

  if (plotpops == "all")
    plotpops <- unique(data$pop.name)

  # Set up headings for params
  params <- make.names(params)

  # set legend size
  nscen <- length(unique(data$scen.name))
  if (nscen < 32){
    legKeySize <- 5
    legTxtSize <- 7
  } else {
    if (nscen > 60) {
      legKeySize <- 2.5
      legTxtSize <- 5
    } else {
      legKeySize <- round(160 / nscen, digits=2)
      legTxtSize <- 5
    }
  }

  fname_root <- if (is.null(scenario)) {
    project
  } else {
    paste0(project, "_", scenario)
  }

  yrmidstdat <- subset(data, Year <= yrmid & pop.name %in% plotpops)

  r.line_plot_year_mid <- list()
  i <- 0
  if (save2disk == T) {
    dir.create(dir.out, showWarnings=FALSE)
    pdf(paste0(dir.out, "/", fname_root, "_", "YearMidVsParams.pdf"))
  }

  for (param in params) {
    i <- i + 1
    root <- paste(fname_root, param, sep="_")
    g <- ggplot(yrmidstdat, aes_string(x="Year", y=param)) +
      geom_line(aes(color=scen.name)) +
      facet_grid(pop.name~.)
    plot <- g +
      theme(panel.margin=unit(0.2, "inches"),
            legend.text=element_text(size=legTxtSize),
            legend.key.size=unit(legKeySize, "mm")) +
      scale_colour_discrete(name="Scenarios")
    print(plot)
    r.line_plot_year_mid[[i]] <- plot
    assign(paste(root, "_", "Mid", "plot", sep=""), g)
  }
  pat <- paste0(fname_root, "_", ".*", "_", "Mid", "plot")
  if (save2disk == T) {
    dev.off()
    save(list=(ls(pattern=pat)),
         file=paste0(dir.out, "/", fname_root, "_", "YearMidVsParams.rda"))
  }

  names(r.line_plot_year_mid) <- ls(pattern=pat)
  return(r.line_plot_year_mid)
}


#' Dot plots of mean Vortex parameters.
#'
#' \code {dot_plot} generates dot plots of mean parameter values for each population
#' (row) at each year value requested with 'yrs' (columns).
#' Bars represent standard deviation
#'
#'TODO examples
#' @param yrs The years to be included in the plot
#' @param setcolour Variable to be used to set colours of data, default: scen.name
#' @inheritParams line_plot_year
#' @return dot plots of mean parameter values with standard deviation
#' @import ggplot2
#' @import grid
#' @export
dot_plot <- function(data=NULL,
                     project=NA,
                     scenario=NA,
                     yrs=c(1,2),
                     params=c("PExtinct", "Nextant", "Het", "Nalleles"),
                     setcolour="scen.name",
                     plotpops=c("all"),
                     save2disk=TRUE,
                     dir.out="Plots") {

  require(ggplot2)
  require(grid)

  if (plotpops == "all")
    plotpops <- unique(data$pop.name)

  # Set up headings for params
  params <- make.names(params)

  fname_root <- if (is.null(scenario)) {
    project
  } else {
    paste0(project, "_", scenario)
  }

  # set legend size
  nLegItems <- length(unique(data[ , setcolour]))
  if (nLegItems < 32){
    legKeySize <- 5
    legTxtSize <- 7
  } else {
    if (nLegItems > 60) {
      legKeySize <- 2.5
      legTxtSize <- 5
    } else {
      legKeySize <- round(160 / nLegItems, digits=2)
      legTxtSize <- 5
    }
  }

  # Vector of SD names for params
  SDname<-function(parSD) paste("SD.", parSD, ".", sep="")
  SD<-sapply(params, SDname)
  if ("r.stoch" %in% params) SD["r.stoch"] <- "SD.r."

  # dot plots by pops & yrs of mean params with (SD) bars
  popstdat <- subset(data, pop.name %in% plotpops)

  r.dot_plot <- list()
  if (save2disk == T) {
    dir.create(dir.out, showWarnings=FALSE)
    pdf(paste(dir.out, "/", fname_root, "_", "dot_plots.pdf", sep=""))
  }

  for (i in 1:length(params)) {
    yrstdat <- subset(popstdat, Year %in% yrs)
    root <- paste(fname_root, "_", params[i], sep="")
    min <- yrstdat[params[i]] - yrstdat[SD[i]]
    names(min) <- "min"
    max <- yrstdat[params[i]] + yrstdat[SD[i]]
    names(max) <- "max"
    yrstdat <- cbind(yrstdat,min, max)
    limits <- aes(ymax=max, ymin=min)
    d <- ggplot(yrstdat,
                aes_string(color=setcolour, x="scen.name", y=params[i])) +
      geom_point() +
      theme(axis.text.x=element_text(angle=-90, size=5, vjust=1)) +
      xlab("Scenario") +
      facet_grid(pop.name~Year) +
      geom_errorbar(limits, width=0.15) +
      theme(panel.margin=unit(0.2, "inches")) +
      if (setcolour == "scen.name") {
        theme(legend.position="none") } else {
          theme(legend.position="right",
                legend.text=element_text(size=legTxtSize),
                legend.key.size=unit(legKeySize, "mm"))
        }
    print(d)
    assign(paste(root, "_", "dot_plot", sep=""), d)
    r.dot_plot[[i]] <- d
  }
  pat <- paste0(fname_root, "_", ".*", "_", "dot_plot")
  if (save2disk == T) {
    dev.off()
    save(list=(ls(pattern=pat)),
         file=paste0(dir.out, "/", fname_root, "_", "dot_plots.rda"))
  }

  names(r.dot_plot) <- ls(pattern=pat)
  return(r.dot_plot)
}

#' Generates a matrix of scatter plots
#'
#' \code{m_scatter} generates a matrix of pairwise scatter plots
#'
#'TODO examples
#'
#' @param data  The output from \code {collate_dat}, the long format of the
#' output from \code {collate_run} or the output from \code{con_l_yr}
#' @param data_type The type of input data. Possible options are "dat", "yr" or "run"
#' @param lookup A table to add relevant variable matched using the scenarios
#' names
#' @param yr The year to be plotted
#' @param popn The sequential number of the population (in integer)
#' @param param The parameter to be plotted in the last raw
#' @param vs The parameters to be plotted
#' @param fname The name of the files where to save the output
#' @inheritParams line_plot_year
#' @return a matrix of scatter plots
#' @import data.table
#' @export
m_scatter <- function (data=NULL,
                      data_type="dat", # possible options are "dat", "yr" or "run"
                      lookup=NA,
                      yr=1,
                      popn=1,
                      param="N",
                      vs=NA,
                      save2disk=TRUE,
                      fname=NULL,
                      dir.out="Plots") {
  require(GGally)
  require(data.table)

  # Set up headings for param
  param <- make.names(param)

  if (data_type == "dat") {
    data <- data.table(data)

    # Replace col names if dat.
    # This will be changed directly in collate_dat in future versions
    setnames(data, c("pop.name", "scen.name"), c("Population", "Scenario"))

    # Thin data to use less memory
    setkey(data, Year)
    data <- data[J(yr), ]
    pop <- data[ , levels(Population)][popn]
  } else {
    if (data_type == "yr") {
      setkey(data, Year)
      data <- data[J(yr), ]
      pop <- paste0("pop", popn) # The pop is selected with its number.
    } else {
      data <- data.table(data)
      pop <- data[ , levels(Population)][popn]
    }
  }

  setkey(data, Population)
  data <- data[.(pop), ] # Select pop
  suppressWarnings(if (!is.na(lookup)) {
    data <- join(data, lookup, by='Scenario', type="left")
  } )
  corrMtrx <- ggpairs(data[ , c(vs, param), with=FALSE], alpha=0.2,
                      lower = list(continuous="smooth",
                                   params=c(method="loess", colour="red")))

  if (save2disk == T) {
    dir.create(dir.out, showWarnings=FALSE)
    pdf(file=paste0(dir.out, "/", fname, "m_scatter_plots.pdf"))
    print(corrMtrx)
    dev.off()
    save(corrMtrx, file=paste0(dir.out, "/", fname, "m_scatter_plots.rda"))
  }
  return(corrMtrx)
}


#' Create a table that summarises simulation parameters for each scenario
#'
#' \code{lookup_table} creates a table that summarises simulation parameters.
#' The final table will have a line for each scenario and one column for each
#' parameter requested with \code{SVs}. If the Vortex project has more than one
#' populations, the user has to indicate a populaiton ot be used as reference
#'
#'TODO examples
#'
#' @param data The output from \code {collate_dat}
#' @param project Vortex project name (used to name the output)
#' @param scenario Vortex scenario name (used to name the output)
#' @param pop The name of the pop to be used as reference
#' @param SVs The parameters to include in the table
#' @param save2disk Whether to save the output to disk, default: TRUE
#' @param dir.out The local path to store the output. Default: ProcessedData
#' @return a data.frame with scenario names and parameter values
#' @import data.table
#' @export
# the user may have set up his own SVs so he/she has to provide a list of SVs
lookup_table <-  function(data=NULL,
                          project=NULL,
                          scenario=NULL,
                          pop="Population 1",
                          SVs=c("SV1"),
                          save2disk=TRUE,
                          dir.out="ProcessedData") {
  require(data.table)

  fname <- if (is.null(scenario)) {
    project
  } else {
    paste0(project, "_", scenario)
  }

  LookUpT <- data.table(data)
  setkey(LookUpT, pop.name, Year)
  LookUpT <- LookUpT[J(pop, 0), c("scen.name", SVs), with=FALSE]
  setnames(LookUpT, "scen.name", "Scenario")

  if(save2disk == T) {
    df2disk(LookUpT, dir.out, fname, "_LookUpT")
  }
  return(LookUpT)
}

#' Calculate the effective population size (Ne)
#'
#' \code{Ne} calculates the effective population size (Ne) for several scenarios
#' based on the loss of genetic diversity (expected heterozygosity) using the
#' temporal approach.
#'
#' @param data The output from \code {collate_dat}
#' @param scenario A vector of scenario names for which Ne need to be calculated,
#' default: "all"
#' @param gen The generation time express in years
#' @param yr0,yrt The time window to be considered (first and last year respectively)
#' @param save2disk Whether to save the output to disk, deafult: TRUE
#' @param fname The name of the files where to save the output, defult: "Ne"
#' @param dir.out The local path to store the output. Default: DataAnalysis
#' @return a data.table with Ne values
#' @import data.table
#' @export
### Ne (effective pop size based on genetic loss)
Ne <-  function(data=NULL,
                scenarios="all",
                gen=1,
                yr0=1,
                yrt=2,
                save2disk=TRUE,
                fname="Ne",
                dir.out="DataAnalysis") {
  require(data.table)

  # Function definitions
  NeOne <- function(scenario) {
          stdatDT <- data.table(data)
          setkeyv(stdatDT, c("scen.name", "Year"))

          # From yr0+gen
          t0 <- round(yr0 + gen)
          tf <- yrt
          Gt0 <- stdatDT[J(scenario, t0), "Het", with=F]
          Gtf <- stdatDT[J(scenario, tf), "Het", with=F]
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

  if (save2disk == T) {
    df2disk(r.Ne, dir.out, fname)
  }

  message(paste("Effective population size based on loss of gene diversity from year",
                round(yr0 + gen), "to year", yrt))
  message("NOTE: The first year used in the calculations is adjusted using
the generation time provided (yr0 + gen) so that Ne:N ratio can be directly calculated
when passing the same argument values to the Ne and Nadults functions. See documentation for more information")
  return(r.Ne)
}

#' Calculate the harmonic mean of the total number of adults between two years
#'
#' \code{Nadults} calculates, for several scenarios, the harmonic mean of the total
#' number of adults between two years for several scenarios. These can be use to
#' calculate Ne/N ratios where relevant
#'
#' @param data The second element (census_means) of the output from \cose{collate_yr}
#' @param npops_noMeta The total number of populations excluding the metapopulation,
#' default: 1
#' @param appendMeta Whether to calculate data for the metapopulation,
#' default: FALSE
#' @param fname The name of the files where to save the output, defult: "Nadults"
#' @inheritParams Ne
#' @return a data.table with Nb values
#' @import plyr
#' @import data.table
#' @export
Nb <- function (data=NULL,
                scenarios="all",
                npops_noMeta=1,
                appendMeta=FALSE,
                gen=1,
                yr0=1,
                yrt=2,
                save2disk=TRUE,
                fname="Nadults",
                dir.out="DataAnalysis") {

  require(plyr)
  require(data.table)

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
  if (npops_noMeta > 1 & appendMeta == T) {
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
                "and number of breeders (sum of adults) from year", yr0,
                "to year", round(yrt-gen)))
  message("NOTE: The last year used in the calculations is adjusted using
the generation time provided (yrt - gen) so that Ne:N ratio can be directly
calculated when passing the same argument values to the Ne and Nadults functions.
See documentation for more information")

  if (scenarios == "all")
    scenarios <- data[ , unique(Scenario)]

  setkey(lcensusMeans, Scenario)
  slcensusMeans <- lcensusMeans[J(scenarios), ]
  setkey(slcensusMeans, Year)
  slcensusMeans <- slcensusMeans[.(yr0:round(yrt - gen)), ]

  slcensusMeans[ , Nb := sum(.SD), .SDcols=c("AM", "AF"),
                by=list(Scenario,Population,Year)]

  harm.means <- slcensusMeans[ , lapply(.SD, HarmMean), .SDcols=c("Nb", "N"),
                              by=list(Scenario,Population)]
  message("Done!")

  if (save2disk == T) {df2disk(harm.means, dir.out, fname)}
  return(harm.means)
}


#' Basic statistics to compare scenarios
#'
#' \code{Pairwise} conducts pairwise comparisons against a baseline scenario
#' uing sensitivity coefficients and strictly standardised mean difference. It also
#' ranks scenarios (and/or parameters when relevant) accordingly
#'
#' @param data A data.frame generated by \code{collate_dat}
#' @param project The Vortex project name
#' @param scenario The ST Vortex scenario name or the scenario that should be
#' used as baseline if simulations were not conducted with the ST module
#' @param params Parameters to be compared,
#'  default: c("PExtinct", "Nextant", "Het", "Nalleles")
#' @param yrs The year(s) to be analysed, default: "max"
#' @param ST Whether files are from sensitivity analysis (TRUE),
#'  or not (FALSE, default)
#' @param type Type of ST. Possible options are: "Sampled",
#' "Latin Hypercube Sampling", "Factorial" or "Single-Factor"
#' @param group.mean Whether calculate the mean of the statistics
#'  (SSMD and Sensitivity Coeffiecient) by group
#' @param SVs Parameters to be used to group scenarios, default: NA
#' @param save2disk Whether to save the output to disk, deafult: TRUE
#' @param dir.out The local path to store the output. Default: DataAnalysis/Pairwise
#' @return several output. See documentation for details
#' @import data.table irr
#' @export
#'
Pairwise <-  function(data=NULL,
                      project=NULL,
                      scenario=NULL,
                      params=c("PExtinct", "Nextant", "Het", "Nalleles"),
                      yrs="max",
                      ST=FALSE,
                      type=NA,
                      group.mean=FALSE,
                      SVs=NA,
                      save2disk=TRUE,
                      dir.out="DataAnalysis/Pairwise") {

  require(data.table)
  require(irr)

  # Function definitions
  SEname <- function(par) paste("SE.", par, ".", sep="")
  SDname <- function(parSD) paste("SD.", parSD, ".", sep="")
  naming.coef <- function(naming) paste("C", "_", naming, yr, sep="")
  naming.ssmd <- function(naming.ssmd)
                      paste("SSMD", "_", naming.ssmd, yr, sep="")
  pval <- function(x) pnorm(abs(x), lower.tail=FALSE)

  # Flag (false) columns where all entries are NA
  FColsAllNA <- function(lranks)
                        apply(lranks, 2, function(chk) !all(is.na(chk)))

  # Error handling
  suppressWarnings(if (!yrs == "max" & !is.numeric(yrs))
    stop("invalid value(s) for 'yrs' "))

  fname <- if (ST == TRUE) {
    paste(project, "_", scenario, sep="")
  } else {
    project
  }

  # set yrs to max
  suppressWarnings(if (yrs == "max") {yrs <- max(data$Year)})

  # set group.mean if needed
  if (ST == TRUE & type == "Single-Factor" & length(SVs) > 1) {group.mean <- T}

  # Set up headings for params and SE and SD
  params <- make.names(params)
  SE <- sapply(params, SEname)
  if ("r.stoch" %in% params) {SE["r.stoch"] <- "SE.r."}

  SD <- sapply(params, SDname)
  if ("r.stoch" %in% params) {SD["r.stoch"] <- "SD.r."}

  # Create a dataframe for the base scenario.
  if (ST == TRUE ) {
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
  kendall.out <- list(coef=NULL, SSMD=NULL)
  # NOTE : kendall function handles na listwise
  kendall.out$coef <- capture.output(
    cat("Rank comparison of sensitivity coefficients", "\n"),
    lapply(ranks.sc.fin, kendall, TRUE))
  kendall.out$SSMD <- capture.output(
    cat("Rank comparison of SSMD", "\n"), lapply(ranks.ssmd.fin, kendall, TRUE))
  if (save2disk == T) {
    # write results
    df2disk(table.coef, dir.out, fname, ".coef.table")
    df2disk(ssmd.table, dir.out, fname, ".SSMD.table")
    df2disk(ssmd.table.pvalues, dir.out, fname, ".SSMD.table.pvalues")
    df2disk(ranks.sc, dir.out, fname, ".ranks.sc")
    df2disk(ranks.ssmd, dir.out, fname, ".ranks.SSMD")
    capture.output(print(kendall.out, quote=F),
                   file=paste0(dir.out, "/", fname,".kendall.txt"))
  }
  # Collate results in a list
  r.OneWay<-list(coef.table=table.coef,
                 SSMD.table=ssmd.table,
                 SSMD.table.pvalues=ssmd.table.pvalues,
                 ranks.coef=ranks.sc,
                 ranks.SSMD=ranks.ssmd,
                 Kendall=kendall.out)
  # if group.mean == TRUE calculate the mean for each SV and rank SVs
  if (group.mean == TRUE) {

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
    kendall.mean.out <- list(coef=NULL, SSMD=NULL)
    kendall.mean.out$coef <- capture.output(
      print("Rank comparison of mean sensitivity coefficients"),
      lapply(ranks.msc.fin, kendall, TRUE))
    kendall.mean.out$SSMD <- capture.output(
      print("Rank comparison of mean SSMD"),
      lapply(ranks.mssmd.fin, kendall, TRUE))

    if (save2disk == T) {
      df2disk(mean.coef.table, dir.out, fname, ".mean.coef.table")
      df2disk(mean.ssmd.table, dir.out, fname, ".mean.SSMD.table")
      df2disk(mean.ssmd.table.pvalues, dir.out, fname, ".mean.SSMD.table.pvalues")
      df2disk(ranks.msc, dir.out, fname, ".ranks.msc")
      df2disk(ranks.mssmd, dir.out, fname, ".ranks.mSSMD")
      capture.output(print(kendall.mean.out, quote=F),
                     file=paste0(dir.out, "/", fname, ".Kendall.means.txt"))
    }

    # Collate results for means
    r.OneWay$mean.coef.table <- mean.coef.table
    r.OneWay$mean.SSMD.table <- mean.ssmd.table
    r.OneWay$mean.SSMD.table.pvalues <- mean.ssmd.table.pvalues
    r.OneWay$ranks.mean.coef <- ranks.msc
    r.OneWay$ranks.mean.SSMD <- ranks.mssmd
    r.OneWay$Kendall.means <- kendall.mean.out
  }
  return(r.OneWay)
}


#' Search for the best regression model(s) given a list of predictors
#'
#'
#' \code{fit_regression} fit either a Generilised Linear Model or a betareg model
#' to the data and search for the best model(s) using the R package glmulti.
#' \code{fit_regression} will conduct an exaustive search if ncand is less or
#' equal to the number of candidate models, otherwise it will use a genetic
#' search method (see glmulti documentations for more details about the search
#' methods)
#'
#' @param data The long format of census (from \code{conv_l_yr}) or run (lrun,
#' the second element) of the output from \code{collate_run}
#' @param lookup (Optional) A look-up table where the scenario names are listed
#' together with the (missing) variables needed to fit the regression models
#' @param census Whether the input is census data
#' @param yr The year that has to be used in the analysis if census=TRUE
#' @param project Vortex project name
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
#' The number of models to be looked for, i.e. the size of the returned confidence
#' set.)
#' @param save2disk Whether to save the output to disk, default: TRUE
#' @param dir.out The local path to store the output. Default: DataAnalysis/FitRegression
#' @return a glmulti object with the best models found
#' @import glmulti data.table plyr
#' @export
fit_regression <-  function(data=NULL,
                            lookup=NA,
                            census=T,
                            yr=NA,
                            project=NULL,
                            scenario=NULL,
                            popn=NA,
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
                            dir.out="DataAnalysis/FitRegression") {
  # Load required packages. Loading of betareg and R.utils are delayed after
  # it has been checked that they actually are needed to avoid wasting time
  # loading packages that are not used.
  require(glmulti)
  require(data.table)
  require(plyr)

  # Functions

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
    data <- join(data, lookup, by='Scenario', type="left")
  } )

  # convert data.table
  data <- data.table(data)

  # select yr if census and create pop vector
  if (census == T) {
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
  } else {
    require(betareg)
  }

  # Plot param
  # Tried DT[, plot(GeneDiv)] to call the plot directily within data.table
  # but when I use param rather than the column names DT return a datatable
  # rather than a vector (as it does with GeneDiv)
  paramvalues <- data[[param]]
  name <- paste(project, scenario, param, sep="_")
  if (save2disk == T) {
    dir.create(dir.out, showWarnings=FALSE)
    pdf(paste0(dir.out, "/", paste(name, "histogram.pdf", sep="_")))
  }

  hist(paramvalues,  main=paste("Histogram of", param), xlab=param)
  if (save2disk == T) {dev.off()}

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
      message("Loading library for package: R.utils...")
      require(R.utils)

      # Changed the c value for qaic search with glmulti
      setOption("glmulti-cvalue", summary(glm2)$dispersion)
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
                                               plotty=F, report=F, level=l,
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
                                                plotty=F, report=F, level=l,
                                                name=name, fitfunc=betareg,
                                                link=links[linkpos],
                                                na.action=na.omit)))
  }
  message("Done! Elapsed time:")
  print (tnp)

  if (save2disk == T) {
    message("Best models saved to disk in the file ...best.mod.rda")
    save(best.mod, file=paste0(dir.out, "/", name, "_best.mod.rda"))
    pdf(paste0(dir.out, "/", name, "_IC_plot.pdf"))
    plot(best.mod, type="p")
    dev.off()
  }
  return(best.mod)
}
