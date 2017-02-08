
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
#' # Save collated data as .Rda and .txt
#' \dontrun{
#' # Read data from all .stdat of the project 'Pacioni_et_al' and the ST scenario
#' #   'ST_Classic'. Store the output in the object 'woylie.st.classic' and save
#' #   to disk
#' woylie.st.classic <- collate_dat("Pacioni_et_al", 3, scenario = "ST_Classic",
#'                      dir_in = pac.dir, save2disk=TRUE)
#' }

collate_dat <- function(project, runs,
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
    if (save2disk) {df2disk(d, dir_out, fname, "_data")}
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
#' @param scenario The scenario name
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
collate_run <- function(project,
                        scenario,
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
                               colClasses="character", comment.char="")
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

    if (save2disk) {
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
collate_yr <-  function(project,
                        scenario,
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

    if (save2disk) {
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

collate_proc_data <- function(data,
                              save2disk=TRUE,
                              dir_out="ProcessedData"){

    # TODO use data.table::rbindlist
    dfs <- plyr::rbind.fill(data)

    if (save2disk) {df2disk(dfs, dir_out, "CombinedDB")}

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
conv_l_yr <- function(data,
                      npops_noMeta=1,
                      appendMeta=FALSE,
                      project,
                      scenario,
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
    if (npops_noMeta > 1 & appendMeta) {
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

    if (save2disk) {
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
lookup_table <-  function(data,
                          project,
                          scenario,
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

    if(save2disk) {
        df2disk(LookUpT, dir_out, fname, "_LookUpT")
    }
    return(LookUpT)
}
