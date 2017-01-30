
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
        dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)
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
        dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)
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
    if (nLegItems < 32) {
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
    SDname <- function(parSD) paste("SD.", parSD, ".", sep="")
    SD <- sapply(params, SDname)
    if ("r.stoch" %in% params) SD["r.stoch"] <- "SD.r."

    # dot plots by pops & yrs of mean params with (SD) bars
    popstdat <- subset(data, pop.name %in% plotpops)

    r.dot_plot <- list()
    if (save2disk == TRUE) {
        dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)
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
    corrMtrx <- GGally::ggpairs(data[ , c(vs, param), with=FALSE],
                                axisLabels='internal',
                                lower = list(continuous="smooth"),
                                mapping = aes(colour="red", alpha=0.2))

    if (save2disk == T) {
        dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)
        pdf(file=paste0(dir_out, "/", fname, "m_scatter_plots.pdf"))
        print(corrMtrx)
        dev.off()
        save(corrMtrx, file=paste0(dir_out, "/", fname, "m_scatter_plots.rda"))
    }
    return(corrMtrx)
}

