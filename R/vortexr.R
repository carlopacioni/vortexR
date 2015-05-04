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
#' \code{vortexR} has a number of functions that are useful during the development
#' of a Vortex project and for its analysis after completion. \code{vortexR}
#' makes it easy to automatise the creation of plots and computation of
#' basic statistics to inspect the effect of changes carried out in the Vortex
#' project. Once the project development is completed, the same framework used
#' in \code{vortexR} during the development of the project can be
#' refined and extended to include more advanced statistical analyses or can be
#' easily included in Markdown documents for the creation of reports (by
#' converting them into pdf) or update web-pages.
#'
#' The use of \code{vortexR} ensures reproducibility and standardises analytical
#' approaches in population viability analysis.
#'
#' @section Documentations:
#' Use \code{help(package = "vortexR")} for a list of \code{vortexR} functions and their
#' specific documentations.
#'
#' A more detailed description of the package and functions can be opened with:
#' \code{vignette(package="vortexR", topic="User-guide")}.
#'
#' More vignettes may be come available in the future. Use \code{vignette(package="vortexR")}
#' to see all the available vignettes.
#'
#' @section Citation:
#' If you use \code{vortexR}, please cite:
#' Pacioni, C., and Mayer, F. (in prep). vortexR: an R package for post Vortex
#' simulation analysis.
#'
#' @section Get in touch:
#' Please, use \url{https://github.org/carlopacioni/vortexR/issues} to report any
#' issue/bug with vortexR. If unsure, or for feedback, contact me at:
#' carlo.pacioni 'at' gmail.com.
#'
#' @section Publications:
#' \emph{Below there are listed a few publications that used \code{vortexR}.}
#'
#' Campbell et al. (in press). Assessing the economic benefits of starling
#'  detection and control to Western Australia. Australasian Journal of
#'  Environmental Management.
#'
#' Pacioni, C., Spencer, P.B.S., Lacy, R.C., and Wayne, A.F. (in prep)
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
#' @title Collated results from Vortex scenarios - Campbell et al (in press)
#' @description A dataset with the results from the main Vortex scenarios used
#'  in Campbell et al (in press). Vortex outputs, from the project named
#'  "Starlingv3PopBased" (.dat files), were collated with \code{collate_dat}.
#' @usage data(sta.main)
#' @format a \code{data.frame} with 1632 observations of 44 variables.
#' @source Campbell et al. In press. Assessing the economic benefits of starling
#'  detection and control to Western Australia. Australasian Journal of
#'  Environmental Management
NULL

#' @name sta.evy5
#' @title Collated results from Vortex scenarios - Campbell et al (in press)
#' @description A dataset with the results from Vortex scenarios used in
#'  Campbell et al (in press) to simulate major application of control measures
#'  in every 5 year cycle. Vortex outputs, from the project named
#'  "Starlingv3PopBased" and the sensitivity test scenario "MReductEvy5"
#'  (.stdat files), were collated with \code{collate_dat}.
#' @usage data(sta.evy5)
#' @format a \code{data.frame} with 1020 observations of 47 variables.
#' @source Campbell et al. In press. Assessing the economic benefits of starling
#'  detection and control to Western Australia. Australasian Journal of
#'  Environmental Management
NULL

#' @name sta.evy5.b11
#' @title Collated results from Vortex scenarios - Campbell et al (in press)
#' @description A dataset with the results from Vortex scenarios used in Campbell
#'  et al (in press) to simulate major application of control measures in every
#'  5 year cycle, maintaining 2011 levels of investment. Vortex outputs, from
#'  the project named "Starlingv3PopBased" and the sensitivity test scenario
#'  "MReduction_B11_09Evy5" (.stdat files), were collated with \code{collate_dat}.
#' @usage data(sta.evy5.b11)
#' @format a \code{data.frame} with 1020 observations of 47 variables.
#' @source Campbell et al. In press. Assessing the economic benefits of starling
#'  detection and control to Western Australia. Australasian Journal of
#'  Environmental Management
NULL

#' @name pac.clas
#' @title Collated results from Vortex scenarios - Pacioni et al. (in prep)
#' @description Data from Pacioni et al. (in prep) used to conduct a sensitivity
#'  analysis on demographic parameters. Vortex outputs, from the project named
#'  "Pacioni_et_al" and (Single-Factor) sensitivity test scenario "ST_Classic"
#'  (.stdat files), were collated with \code{collate_dat}.
#' @usage data(pac.clas)
#' @format a \code{data.frame} of 2904 observations of 68 variables.
#' @source Pacioni, C., Spencer, P.B.S., Lacy, R.C., and Wayne, A.F. (in prep)
#'  Predators and genetic fitness: key threatening factors for the conservation
#'  of bettong species. Pacific Conservation Biology.
NULL

#' @name pac.lhs
#' @title Collated results from Vortex scenarios - Pacioni et al. (in prep)
#' @description Data from Pacioni et al. (in prep) used to conduct a sensitivity
#'  analysis on demographic parameters. Vortex outputs, from the project named
#'  "Pacioni_et_al" and (Latin Hypercube Sampling) sensitivity test scenario
#'  "ST_LHS" (.stdat files), were collated with \code{collate_dat}.
#' @usage data(pac.lhs)
#' @format A \code{data.frame} of 6171 observations of 68 variables.
#' @source Pacioni, C., Spencer, P.B.S., Lacy, R.C., and Wayne, A.F. (in prep)
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
#' @source Pacioni, C., Spencer, P.B.S., Lacy, R.C., and Wayne, A.F. (in prep)
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
#' @source Pacioni, C., Spencer, P.B.S., Lacy, R.C., and Wayne, A.F. (in prep)
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
#' @source Pacioni, C., Spencer, P.B.S., Lacy, R.C., and Wayne, A.F. (in prep)
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
#' @source Pacioni, C., Spencer, P.B.S., Lacy, R.C., and Wayne, A.F. (in prep)
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
#' @source Pacioni, C., Spencer, P.B.S., Lacy, R.C., and Wayne, A.F. (in prep)
#'  Predators and genetic fitness: key threatening factors for the conservation
#'  of bettong species. Pacific Conservation Biology.
NULL

#' @name pac.clas.pairw
#' @title Results of pairwise comparisons of simulation scenarios
#' @description Results of pairwise comparisons of simulation scenarios included
#' in the sensitivity test scenario "ST_Classic" using \code{pairwise}.
#' @usage data(pac.clas.pairw)
#' @format A named list of 12 \code{element}s. See documentation for details.
#' @source Pacioni, C., Spencer, P.B.S., Lacy, R.C., and Wayne, A.F. (in prep)
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
#'  than those presented in the paper (only 3 runs for 120 'Vortex-years'.
#' @usage # To retrieve the path to the files use:
#'  system.file("extdata", "pacioni", package="vortexR")
#' @format One .xml file and several .run and .stdat files.
#' @source Pacioni, C., Spencer, P.B.S., Lacy, R.C., and Wayne, A.F. (in prep)
#'  Predators and genetic fitness: key threatening factors for the conservation
#'  of bettong species. Pacific Conservation Biology.
NULL

#' @name campbell
#' @title Raw Vortex outputs from Campbell et al (in press)
#' @description A folder with Vortex outputs from Campbell et al (in press).
#' @usage # To retrieve the path to the files use:
#'  system.file("extdata", "campbell", package="vortexR")
#' @format Several .dat and .stdat files.
#' @source Campbell et al. In press. Assessing the economic benefits of starling
#'  detection and control to Western Australia. Australasian Journal of
#'  Environmental Management
NULL

#------------------------------------------------------------------------------#
# Exported helper functions
#------------------------------------------------------------------------------#

#' Save a data.frame as both Rdata and CSV
#'
#' \code{df2disk} saves to disk a given data.frame as both Rdata and CSV with a
#' given name and optional name postfix to a given location.
#'
#' \code{df2disk} is used by the \code{collate_} functions when the operator
#' chooses to save2disk.
#'
#' @param df A data.frame
#' @param dirpath The destination path for written files, will be created if
#'   necessary
#' @param fname The file name
#' @param postfix An optional name postfix
#'
#' @examples
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
#' File names are sorted using gtools::mixedsort().
#'
#' @param path The directory to search in
#' @param pattern A pattern to match file names
#' @param fn_name The vortexR function name for verbose messages
#' @param verbose Progress messages, default: FALSE
#' @return A character vector of fully qualified file paths
get_file_paths <- function(path, pattern, fn_name, verbose=FALSE){
  files <- gtools::mixedsort(list.files(path=path,
                                        pattern=pattern,
                                        full.names=TRUE))

  if (length(files) == 0) {
    stop(paste0("ERROR vortexR::", fn_name," found no files",
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


#' Return a prefixed and repeated string of character
#'
#' @param chars A string of characters (popvalue)
#' @param times The number of repetitions (ncolpop), default: 1
#' @param prefix A text prefix, default: ""
PrefixAndRepeat <- function(chars, times=1, prefix="") {
  rep(paste0(prefix, chars), times)
}


#' Compile iterations from one .yr file
#'
#' Compile iterations from one .yr file and add a column with scenario names and
#' one with iteration number
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
# Data processing
#------------------------------------------------------------------------------#

#' Collate one local Vortex output file into a data.frame
#'
#' \code{collate_one_dat} parses one Vortex .dat or .stdat file, and returns the
#' data within as one data.frame.
#'
#' @param filename The fully qualified filename of a Vortex .dat or .stdat file
#' @param runs The number of simulation runs
#' @param verbose Progress messages, default: FALSE
#' @return A data.frame with data from one .dat or .stdat file and
#'  population/scenario names as factors
#' @export
#' @examples
#' # Pacioni et al. example files. See ?pacioni for more details.
#' pac.dir <- system.file("extdata", "pacioni", package="vortexR")
#' f <- paste0(pac.dir, "/", "Pacioni_et_al_ST_Classic(Base).stdat")
#' one.st.classic <- collate_one_dat(f, 3)
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

#' Collate Vortex .dat or .stdat output files into one data.frame
#'
#' \code{collate_dat} collates all Vortex output files matching a given project
#' name (and scenario name when relevant) in a given directory into one
#' data.frame using \code{collate_one_dat}.
#'
#' The number of Vortex simulation runs has to be specified, as it cannot be
#' inferred by vortexR.
#'
#' To read Vortex output files from Sensitivity Testing with the extension
#' ".stdat", specify the scenario name. If \code{scenario=NULL}, all files with
#' extension .dat, matching the project name will be imported.
#'
#' \code{dir_in} may contain other files; only files matching the project (and,
#' optionally, the scenario) name will be read.
#'
#' \code{dir_out} is created within the working directory unless a full path is
#' provided.
#'
#' If no matching files are found in the given directory, an error is reported.
#'
#' When \code{verbose=TRUE} the progress (i.e. the file being read) is reported
#' on screen.
#'
#' @param project The Vortex project name to be imported
#' @param scenario The scenario name if ST, default: NULL
#' @param runs The number of Vortex simulation runs
#' @param dir_in The local folder containing Vortex files, default: NULL. If
#'   not specified, will fall back to use current working directory.
#' @param save2disk Whether to save the data as rda and csv, default: TRUE
#' @param dir_out The local path to store the output. Default: ProcessedData
#' @param verbose Progress messages, default: TRUE
#' @return a data.frame with data from all matching Vortex files or NULL
#' @export
#' @examples
#' # Using Campbell et al. and Pacioni et al. example files.
#' # See ?pacioni and ?campbell for more details on example files.
#' camp.dir <- system.file("extdata", "campbell", package="vortexR")
#' pac.dir <- system.file("extdata", "pacioni", package="vortexR")
#'
#' # Campbell example, project "Starlingv3PopBased" (.dat)
#' starling <- collate_dat("Starlingv3PopBased", 10000,
#'             dir_in=camp.dir, save2disk=FALSE)
#'
#' # Read data from all .stdat of the project 'Pacioni_et_al' and the ST scenario
#' #   'ST_Classic' and store the output in the object 'woylie.st.classic'
#' woylie.st.classic <- collate_dat("Pacioni_et_al", 3, scenario = "ST_Classic",
#'                      dir_in = pac.dir, save2disk=FALSE)
#'
#' # Read data from all .stdat of the project 'Pacioni_et_al' and the ST scenario
#' #   'ST_LHS' and store the output in the object 'woylie.lhs'
#' woylie.lhs <- collate_dat("Pacioni_et_al", 3, scenario = "ST_LHS",
#'            dir_in = pac.dir, save2disk=FALSE)
#'
#' # Save collated data as .Rda and .txt
#' \dontrun{
#' # Read data from all .stdat of the project 'Pacioni_et_al' and the ST scenario
#' #   'ST_LHS' and store the output in the object 'woylie.lhs' and save output
#' # to disk
#' woylie.lhs <- collate_dat("Pacioni_et_al", 3, scenario = "ST_LHS",
#'            dir_in = pac.dir)
#' }
collate_dat <- function(project=NULL, runs,
                        scenario=NULL,
                        dir_in=NULL,
                        save2disk=TRUE,
                        dir_out="ProcessedData",
                        verbose=TRUE){

  if (is.null(scenario)) {
    fname <- project
    pat <- paste0("^", fname, ".*\\.dat$")
  } else {
    fname <- paste(project, scenario, sep="_")
    pat <- paste0("^", fname, ".*\\.stdat$")
  }

  if (is.null(dir_in)) {dir_in <- getwd()}

  files <- get_file_paths(path=dir_in,
                          pattern=pat,
                          fn_name="collate_dat",
                          verbose=verbose)

  d <- data.frame()
  if (verbose){message("vortexR::collate_dat is parsing:")}
    for (filename in files) {
      if (verbose){message(filename, "\r")}
      d <- rbind(d, collate_one_dat(filename, runs))
    }
    if (save2disk == T) {df2disk(d, dir_out, fname, "_data")}
    return(d)
}


#' Collate Vortex .run output files
#'
#' \code{collate_run} collates all Vortex output files with extension .run
#' matching the project and scenario name in a given directory  into one  named
#' list.
#'
#' \code{dir_in} may contain other files; only files matching the project and
#' the scenario name will be read.
#'
#' \code{dir_out} is created within the working directory unless a full path is
#' provided.
#'
#' If no matching files are found in the given directory, an error is reported.
#'
#' When \code{verbose=TRUE} the progress (i.e. the file being read) is reported
#' on screen.
#'
#' @param project The Vortex project name to be imported
#' @param scenario The scenario name, default: NULL
#' @param npops The total number of simulated populations including the
#'   metapopulation
#' @inheritParams collate_dat
#' @return a list with two elements: run, a data.frame with data from all
#'   Vortex files and lrun, where the same data are re-arranged in long format
#' @export
#' @examples
#' # Using Pacioni et al. example files. See ?pacioni for more details.
#' pac.dir <- system.file("extdata", "pacioni", package="vortexR")
#' # Run collate_run on all .run of the project 'Pacioni_et_al' and
#  # the ST scenario 'ST_LHS' in the selected folder and tore the result in 'run'
#' run <- collate_run("Pacioni_et_al", "ST_LHS", 1, dir_in=pac.dir,
#'                    save2disk=FALSE)
collate_run <- function(project=NULL,
                        scenario=NULL,
                        npops=1,
                        dir_in=NULL,
                        save2disk=TRUE,
                        dir_out="ProcessedData",
                        verbose=TRUE) {

  run <- data.frame()
  lrun <- data.frame()

  if (is.null(dir_in)) {dir_in <- getwd()}
  fname <- paste(project, scenario, sep="_")

  files <- get_file_paths(path=dir_in,
                          pattern=paste0("^", fname,".*\\.run$"),
                          fn_name="collate_run",
                          verbose=verbose)

  if (verbose){message("vortexR::collate_run is parsing:")}
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
    df2disk(run, dir_out, fname, "_run")
    df2disk(lrun, dir_out, fname, "_lrun")
  }

  r.RunST <- list(run=run, lrun=lrun)
  return(r.RunST)
}


#' Collate Vortex .yr output files
#'
#' \code{collate_yr} collates all the .yr output from Vortex matching the project
#' and scenario name into one R object and calculates the mean for each simulated
#' year across all iterations.
#'
#' \code{dir_in} may contain other files; only files matching the project and
#' the scenario name will be read.
#'
#' \code{dir_out} is created within the working directory unless a full path is
#' provided.
#'
#' If no matching files are found in the given directory, an error is reported.
#'
#' When \code{verbose=TRUE} the progress (i.e. the file being read) is reported
#' on screen.
#'
#' @param project The Vortex project name to be imported
#' @param scenario The scenario name, default: NULL
#' @param npops_noMeta The total number of populations excluding the metapopulation,
#' default: 1
#' @inheritParams collate_dat
#' @return a list with two elements: "census", a data.frame with data from all
#' Vortex files and "census_means", a data.table with the mean of each parameter
#' across all iterations for each simulated year
#' @import data.table
#' @export
#' @examples
#' # Using Pacioni et al. example files. See ?pacioni for more details.
#' pac.dir <- system.file("extdata", "pacioni", package="vortexR")
#'
#' # Run collate_yr on all .yr of the project 'Pacioni_et_al' and the ST scenario
#' # 'ST_Classic' in the selected folder and store the result in 'yr.st.classic'
#' yr.st.classic <- collate_yr(project="Pacioni_et_al", scenario="ST_Classic",
#'                             dir_in = pac.dir, save2disk=FALSE)
collate_yr <-  function(project=NULL,
                        scenario=NULL,
                        npops_noMeta=1,
                        dir_in=NULL,
                        save2disk=TRUE,
                        dir_out="ProcessedData",
                        verbose=TRUE) {

  if (is.null(dir_in)) {dir_in <- getwd()}
  fname <- paste0(project, "_", scenario)

  files <- get_file_paths(path=dir_in,
                          pattern=paste0("^", fname, ".*\\.yr$"),
                          fn_name="collate_yr",
                          verbose=verbose)

  censusData <- vector(mode="list", length=length(files))

  if (verbose){message("vortexR::collate_yr is parsing:")}
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

    censusData[[i]] <- rbindlist(one_yr)
  }
  censusAll <- rbindlist(censusData)

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


  setkey(censusAll, Scenario, Year)
  censusMeansDT <- censusAll[ , lapply(.SD, mean), by="Scenario,Year"]
  censusMeansDT <- censusMeansDT[ , Iteration := NULL]

  if (save2disk == TRUE) {
    df2disk(censusAll, dir_out, fname, "_census")
    df2disk(censusMeansDT, dir_out, fname, "_census_means")
  }

  return(list(census=censusAll, census_means=censusMeansDT))
}

#------------------------------------------------------------------------------#

#' Collate processed data generated by any of the 'collate' functions
#'
#' \code{collate_proc_data} collates multiple data frames generated by any of the
#' 'collate' functions. This may be useful when, for example, data generated by
#' different ST scenarios and/or 'standard' scenario runs are to be combined
#' into a unique dataframe that can then be passed to other functions.
#'
#' Only dfs generated by the same function can be combined together. Missing data
#' are filled with NA.
#'
#' When \code{save2disk=TRUE} the output file will be named "CombinedDB".
#' \code{dir_out} is created within the working directory unless a full path is
#' provided.
#'
#' @param data A list where each element is a df from \code{collate_[dat, yr, run]}
#' @inheritParams collate_dat
#' @return A data.frame with the collated data. Missing data are filled with NA.
#' @export
#' @examples
#' # Using Campbell et al. example data. See ?sta.main, ?sta.evy5, ?sta.evy5.b11
#' # for more details.
#' data(sta.main, sta.evy5, sta.evy5.b11)
#' dfs <- list(sta.main, sta.evy5, sta.evy5.b11)
#' combined <- collate_proc_data(dfs, save2disk=FALSE)

collate_proc_data <- function(data=NULL,
                              save2disk=TRUE,
                              dir_out="ProcessedData"){

  # TODO use data.table::rbindlist
  dfs <- plyr::rbind.fill(data)

  if (save2disk == T) {df2disk(dfs, dir_out, "CombinedDB")}

  return(dfs)
}

#' Convert 'census' data into long format
#'
#' \code{conv_l_yr} converts the first element of the output from \code{collate_yr}
#' (census) in long format. This can be then fed into downstream analysis
#' (e.g. fit_regression)
#'
#' \code{yrs} is used to indicate the years to be retained in the output. If more
#' than one year is required, these may be requested defining a numeric vector,
#' e.g. yrs=c(10, 20, 30). All simulated years can be included in the output by
#' passing a numeric vector with all years. For example, assuming that 100 years
#' were simulated, using \code{yrs=1:100} would retain all 100 years in the
#' output.
#'
#' \code{dir_out} is created within the working directory unless a full path is
#' provided.
#'
#' @param data The df 'census' from \code{collate_yr}
#' @param appendMeta Whether to calculate data for the metapopulation,
#'   default: FALSE
#' @param project Vortex project name (used to name the output)
#' @param scenario Vortex scenario name (used to name the output)
#' @param yrs The year(s) that need to be retained in the output
#' @inheritParams collate_yr
#' @return The census data.frame in long format
#' @import data.table
#' @export
#' @examples
#' # Using Pacioni et al. example data. See ?pac.yr for more details.
#' data(pac.yr)
#' lyr.classic <- conv_l_yr(pac.yr[[1]] , yrs=c(60, 120), save2disk=FALSE)
conv_l_yr <- function(data=NULL,
                      npops_noMeta=1,
                      appendMeta=FALSE,
                      project=NULL,
                      scenario=NULL,
                      yrs=c(1, 2),
                      save2disk=TRUE,
                      dir_out="ProcessedData") {

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
  data <- data[.(yrs), ] # Select the yrs

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
    lcensus <- plyr::rbind.fill(lcensus, meta)
    message("Done!")
  }

  if (save2disk == T) {
    df2disk(lcensus, dir_out, paste0(project, "_", scenario), "_lcensus")
  }
  return(lcensus)
}

#' Summary table of simulation parameters
#'
#' \code{lookup_table} creates a table that summarises simulation parameters.
#' The final table will have a line for each scenario and one column for each
#' parameter requested with \code{SVs}.
#'
#' If the name of the populations were changed, the user has to indicate a
#' population to be used as reference, otherwise \code{lookup_table} will look
#' for a population named "Population 1" (i.e. Vortex default name for the first
#' population).
#'
#' \code{lookup_table} reports the values of SVs at year zero. This is done
#' because parameters may take value 'zero' if the relevant population
#' goes extinct.There are cases where Vortex may not evaluate some parameters at
#' year 0. This may happen, for example, when a population is empty at
#' initialization (i.e. the initial population size is zero), or when K is set
#' to zero at the beginning of the simulation. The user should check the values
#' reported and check the Vortex input files if these do not look correct.
#'
#' \code{SVs} can be any variable included in the data, including GS or PS set up
#' in Vortex.
#'
#' @param data The output from \code{collate_dat}
#' @param project Vortex project name (used to name the output)
#' @param scenario Vortex scenario name (used to name the output)
#' @param pop The name of the pop to be used as reference
#' @param SVs The parameters to include in the table
#' @param save2disk Whether to save the output to disk, default: TRUE
#' @param dir_out The local path to store the output. Default: ProcessedData
#' @return A data.frame with scenario names and parameter values
#' @import data.table
#' @export
#' @examples
#' # Using Pacioni et al. example data. See ?pac.clas for more details.
#' data(pac.clas)
#' lkup.st.classic <- lookup_table(data=pac.clas, project="Pacioni_et_al",
#'                    scenario="ST_Classic", pop="Population 1",
#'                    SVs=c("SV1", "SV2", "SV3", "SV4", "SV5", "SV6", "SV7"),
#'                    save2disk=FALSE)
lookup_table <-  function(data=NULL,
                          project=NULL,
                          scenario=NULL,
                          pop="Population 1",
                          SVs=c("SV1"),
                          save2disk=TRUE,
                          dir_out="ProcessedData") {

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
    df2disk(LookUpT, dir_out, fname, "_LookUpT")
  }
  return(LookUpT)
}

#------------------------------------------------------------------------------#
# Data plotting
#------------------------------------------------------------------------------#

#' Line plots of Vortex parameters vs years
#'
#' \code{line_plot_year} generates line plots of the selected Vortex parameters
#' for the selected populations, for all simulated years.
#'
#' Plots are ggplot objects. When \code{save2disk=TRUE} these are saved as .rda
#' and .pdf files
#'
#' @param data A df from \code{collate_dat}
#' @param project Vortex project name (used to name the output)
#' @param scenario Vortex scenario name (used to name the output)
#' @param params Vortex parameters to be plotted,
#' default: c("PExtinct", "Nextant", "Het", "Nalleles")
#' @param plotpops The populations to be included in the plot, default: 'all'
#' @param save2disk Whether to save the output to disk, default: TRUE
#' @param dir_out The local path to store the output. Default: Plots
#' @return Line plot(s)
#' @import ggplot2
#' @import grid
#' @export
#' @examples
#' # Using Pacioni et al. example data. See ?pac.clas for more details.
#' data(pac.clas)
#' lineplot.st.classic <- line_plot_year(data=pac.clas, project="Pacioni_et_al",
#'                        scenario="ST_Classic",
#'                        params=c("PExtinct", "Nextant", "Het", "Nalleles"),
#'                        save2disk=FALSE)
line_plot_year <- function(data=NULL,
                           project=NULL,
                           scenario=NULL,
                           params=c("PExtinct", "Nextant", "Het", "Nalleles"),
                           plotpops=c("all"),
                           save2disk=TRUE,
                           dir_out="Plots") {

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
    dir.create(dir_out, showWarnings=FALSE)
    pdf(paste(dir_out, "/", fname_root, "_", "YearVsParams.pdf", sep=""))
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
         file=paste(dir_out, "/", fname_root, "_", "YearVsParams.rda", sep=""))
  }

  names(r.line_plot_year) <- ls(pattern=paste(fname_root, "_", ".*", "_",
                                            "plot", sep=""))
  return(r.line_plot_year)
}

#' Line plots of Vortex parameters vs years
#'
#' \code{line_plot_year_mid} generates line plots of the selected Vortex parameters
#' for the selected populations, from year zero to yrmid. The purpose of these
#' plots is to 'zoom' in the initial phase of the simulations to better
#' appreciate dynamics of the parameters of interest.
#'
#' Plots are ggplot objects. When \code{save2disk=TRUE} these are saved as .rda
#' and .pdf files
#'
#' @param yrmid The last year to plot
#' @inheritParams line_plot_year
#' @return Line plot(s)
#' @import ggplot2
#' @import grid
#' @export
#' @examples
#' # Using Pacioni et al. example data. See ?pac.clas for more details.
#' data(pac.clas)
#' lineMidPlot.st.classic <- line_plot_year_mid(data=pac.clas,
#'                           project="Pacioni_et_al",
#'                           scenario="ST_Classic",
#'                           yrmid=50,
#'                           params=c("PExtinct", "Nextant", "Het", "Nalleles"),
#'                           save2disk=FALSE)
line_plot_year_mid <-  function(data=NULL,
                                project=NULL,
                                scenario=NULL,
                                yrmid=1,
                                params=c("PExtinct", "Nextant", "Het", "Nalleles"),
                                plotpops=c("all"),
                                save2disk=TRUE,
                                dir_out="Plots") {

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
    dir.create(dir_out, showWarnings=FALSE)
    pdf(paste0(dir_out, "/", fname_root, "_", "YearMidVsParams.pdf"))
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
         file=paste0(dir_out, "/", fname_root, "_", "YearMidVsParams.rda"))
  }

  names(r.line_plot_year_mid) <- ls(pattern=pat)
  return(r.line_plot_year_mid)
}


#' Dot plots of mean Vortex parameters
#'
#' \code{dot_plot} generates dot plots of mean parameter values for each population
#' (row) at each year value requested with 'yrs' (columns). Bars represent
#' standard deviation.
#'
#' Plots are ggplot objects. When \code{save2disk=TRUE} these are saved as .rda
#' and .pdf files
#'
#' \code{yrs} can be a numeric vector of length >= 1 (e.g. \code{yrs=c(50,100)}).
#' Each point in time will be plotted in different columns.
#'
#' If a continuous variable is passed to \code{setcolour}, a continuous gradient
#' of colour will be assigned to the marker (e.g. for example, a scale from blue
#' to black). If a sharp change of colours between different values of a
#' continuous variable is desired, it has to be converted into a factor.
#'
#' @param yrs The years to be included in the plot
#' @param setcolour Variable to be used to set colours of data, default: scen.name
#' @inheritParams line_plot_year
#' @return Dot plots of mean parameter values with standard deviation
#' @import ggplot2
#' @import grid
#' @export
#' @examples
#' # Using Pacioni et al. example data. See ?pac.clas for more details.
#' data(pac.clas)
#' dot <- dot_plot(data=pac.clas, project="Pacioni_et_al", scenario="ST_Classic",
#'                yrs=c(80, 120),
#'                params=c("PExtinct", "Nextant", "Het", "Nalleles"),
#'                save2disk=FALSE)
dot_plot <- function(data=NULL,
                     project=NA,
                     scenario=NA,
                     yrs=c(1,2),
                     params=c("PExtinct", "Nextant", "Het", "Nalleles"),
                     setcolour="scen.name",
                     plotpops=c("all"),
                     save2disk=TRUE,
                     dir_out="Plots") {

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
    dir.create(dir_out, showWarnings=FALSE)
    pdf(paste(dir_out, "/", fname_root, "_", "dot_plots.pdf", sep=""))
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
         file=paste0(dir_out, "/", fname_root, "_", "dot_plots.rda"))
  }

  names(r.dot_plot) <- ls(pattern=pat)
  return(r.dot_plot)
}

#' Generates a matrix of scatter plots
#'
#' \code{m_scatter} generates a matrix of pairwise scatter plots to graphically
#' investigate possible associations between variables.
#'
#' The output from \code{collate_dat} is the preferred input for this function
#' as large datasets will require a long time to be plotted.
#'
#' It may be convenient to pass the dependent variable of a regression model
#' with \code{param} so that all the pairwise scatter plots of this variable
#' will be in one line.
#'
#' @param data  The output from \code{collate_dat}, the long format of the
#' output from \code{collate_run} or the output from \code{con_l_yr}
#' @param data_type The type of input data. Possible options are "dat", "yr" or "run"
#' @param lookup A table to add relevant variable matched using the scenarios
#' names
#' @param yr The year to be plotted
#' @param popn The sequential number of the population (in integer)
#' @param param The parameter to be plotted in the last raw
#' @param vs The parameters to be plotted
#' @param fname The name of the files where to save the output
#' @inheritParams line_plot_year
#' @return A matrix of scatter plots
#' @import data.table
#' @export
#' @examples
#' # Using Pacioni et al. example data. See ?pac.lhs for more details.
#' data(pac.lhs)
#' # Remove base scenario
#' pac.lhs.no.base <- pac.lhs[!pac.lhs$scen.name == "ST_LHS(Base)", ]
#'
#' # Use function lookup_table to obtain correct parameter values at year 0
#' lkup.ST_LHS <- lookup_table(data=pac.lhs.no.base, project="Pacioni_et_al",
#'                             scenario="ST_LHS",
#'                             pop="Population 1",
#'                             SVs=c("SV1", "SV2", "SV3", "SV4", "SV5", "SV6", "SV7"),
#'                             save2disk=FALSE)
#'
#' scatter.plot <- m_scatter(data=pac.lhs.no.base[1:33], data_type="dat",
#'                           lookup=lkup.ST_LHS, yr=120, popn=1, param="Nall",
#'                           vs=c("SV1", "SV2", "SV3"),
#'                           save2disk=FALSE)
m_scatter <- function (data=NULL,
                      data_type="dat", # possible options are "dat", "yr" or "run"
                      lookup=NA,
                      yr=1,
                      popn=1,
                      param="N",
                      vs=NA,
                      save2disk=TRUE,
                      fname=NULL,
                      dir_out="Plots") {

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
    data <- plyr::join(data, lookup, by='Scenario', type="left")
  } )
  corrMtrx <- GGally::ggpairs(data[ , c(vs, param), with=FALSE], alpha=0.2,
                      axisLabels='internal',
                      lower = list(continuous="smooth",
                                   params=c(method="loess", colour="red")))

  if (save2disk == T) {
    dir.create(dir_out, showWarnings=FALSE)
    pdf(file=paste0(dir_out, "/", fname, "m_scatter_plots.pdf"))
    print(corrMtrx)
    dev.off()
    save(corrMtrx, file=paste0(dir_out, "/", fname, "m_scatter_plots.rda"))
  }
  return(corrMtrx)
}

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
Nadults <- function (data=NULL,
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

  if (save2disk == T) {df2disk(harm.means, dir_out, fname)}
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
#' @return Several output. See vignette for details.
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
#' @export
#' @examples
#' # Using Pacioni et al. example data. See ?pac.clas for more details.
#' data(pac.clas)
#' pairw<-pairwise(data=pac.clas, project="Pacioni_et_al", scenario="ST_Classic",
#'                params=c("Nall", "Het"), yrs=c(60,120), ST=TRUE,
#'                type="Single-Factor",
#'                SVs=c("SV1", "SV2", "SV3", "SV4", "SV5", "SV6", "SV7"),
#'                save2disk=FALSE)
pairwise <-  function(data=NULL,
                      project=NULL,
                      scenario=NULL,
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
  kendall.out <- list(SC=NULL, SSMD=NULL)
  # NOTE : kendall function handles na listwise
  kendall.out$SC <- capture.output(
                      cat("Rank comparison of sensitivity coefficients", "\n"),
                      lapply(ranks.sc.fin, irr::kendall, TRUE))
  kendall.out$SSMD <- capture.output(
                      cat("Rank comparison of SSMD", "\n"),
                      lapply(ranks.ssmd.fin, irr::kendall, TRUE))
  if (save2disk == T) {
    # write results
    df2disk(table.coef, dir_out, fname, ".coef.table")
    df2disk(ssmd.table, dir_out, fname, ".SSMD.table")
    df2disk(ssmd.table.pvalues, dir_out, fname, ".SSMD.table.pvalues")
    df2disk(ranks.sc, dir_out, fname, ".ranks.sc")
    df2disk(ranks.ssmd, dir_out, fname, ".ranks.SSMD")
    capture.output(print(kendall.out, quote=F),
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
    kendall.mean.out <- list(SC=NULL, SSMD=NULL)
    kendall.mean.out$SC <- capture.output(
                      print("Rank comparison of mean sensitivity coefficients"),
                      lapply(ranks.msc.fin, irr::kendall, TRUE))
    kendall.mean.out$SSMD <- capture.output(
                      print("Rank comparison of mean SSMD"),
                      lapply(ranks.mssmd.fin, irr::kendall, TRUE))

    if (save2disk == T) {
      df2disk(mean.coef.table, dir_out, fname, ".mean.coef.table")
      df2disk(mean.ssmd.table, dir_out, fname, ".mean.SSMD.table")
      df2disk(mean.ssmd.table.pvalues, dir_out, fname, ".mean.SSMD.table.pvalues")
      df2disk(ranks.msc, dir_out, fname, ".ranks.mSC")
      df2disk(ranks.mssmd, dir_out, fname, ".ranks.mSSMD")
      capture.output(print(kendall.mean.out, quote=F),
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
#' with a Poisson error distribution and if \eqn{c^} (dispersion parameter,
#' calculated as \eqn{residual deviance / df} is larger than (the somewhat
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
#'  The number of models to be looked for, i.e. the size of the returned confidence
#' set.)
#' @param save2disk Whether to save the output to disk, default: TRUE
#' @param dir_out The local path to store the output.
#'  Default: DataAnalysis/FitRegression
#' @return A \code{glmulti} object with the best models found.
#' @references Calcagno, V., and C. de Mazancourt. 2010. glmulti: an R package
#' for easy automated model selection with (generalized) linear models. Journal
#' of Statistical Software 34:1-29.
#'
#' Cribari-Neto, F., and Zeileis, A. (2010) Beta regression in R. Journal of
#' Statistical Software 34(2).
#' @import glmulti data.table betareg
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
fit_regression <-  function(data=NULL,
                            lookup=NA,
                            census=TRUE,
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
  }

  # Plot param
  # Tried DT[, plot(GeneDiv)] to call the plot directily within data.table
  # but when I use param rather than the column names DT return a datatable
  # rather than a vector (as it does with GeneDiv)
  paramvalues <- data[[param]]
  name <- paste(project, scenario, param, sep="_")
  if (save2disk == T) {
    dir.create(dir_out, showWarnings=FALSE)
    pdf(paste0(dir_out, "/", paste(name, "histogram.pdf", sep="_")))
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
    save(best.mod, file=paste0(dir_out, "/", name, "_best.mod.rda"))
    pdf(paste0(dir_out, "/", name, "_IC_plot.pdf"))
    plot(best.mod, type="p")
    dev.off()
  }
  return(best.mod)
}
