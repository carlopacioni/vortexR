#' Cumulative probability of extinction at the end of the simulation
#'
#' \code{Pextinct} calculates the cumulative probability of extinction at the by
#' calculating the proportion of runs in which a population goes extinct for
#' each scenario.
#'
#' \code{Pextinct} then compares each scenario by calculating the strictly
#' standardised mean difference (SSMD, Zhang 2007) and reports this statistic
#' with its associated p values. Raw data are also reported.
#'
#' @param data The long format of run (lrun, the second element) of the output
#'   from \code{collate_run}
#' @inheritParams pairwise
#' @return A list with two elements, a table (\code{data.table}) with the mean
#'   Probability of extinction and its SD, the SSMD and its associated p-value
#'   for each scenario and population, and a table (\code{data.table}) with each
#'   iteration where extinction is coded as one (and zero for no extinction)
#' @references Zhang, X. D. 2007. A pair of new statistical parameters for
#'   quality control in RNA interference high-throughput screening assays.
#'   Genomics 89:552-561.
#' @import data.table
#' @export
#' @examples
#' # Using Pacioni et al. example data. See ?pac.run.lhs for more details.
#' data(pac.run.lhs)
#' Pext <- Pextinct(pac.run.lhs[[2]], project="Pacioni_et_al", scenario="ST_Classic",
#'              ST=TRUE, save2disk=FALSE, dir_out="DataAnalysis/Pextinct")


Pextinct <- function(data, project, scenario, ST=FALSE, save2disk=TRUE,
                     dir_out="DataAnalysis/Pextinct") {
    fname <- if (ST == TRUE) {
        paste(project, "_", scenario, sep="")
    } else {
        project
    }
    data <- data.table(data)

    if (ST == TRUE ) {
        scenario <- grep("(Base)", data[, unique(Scenario)], value=TRUE)
    }

    setkey(data, YrExt)
    data[!.(NA), Ext  := 1]
    data[.(NA), Ext  := 0]
    extTable <- data[, .(Pext=mean(Ext), SD=sd(Ext)), by=.(Scenario, Population)]
    setkey(extTable, Scenario)
    Base <- extTable[scenario, .(Pext, SD), by=Population]
    setnames(Base, c("Pext", "SD"), c("base", "SDbase"))
    extTable <- merge(extTable, Base, by="Population")
    extTable[, SSMD := (base - Pext) / sqrt(SD^2 + SDbase^2)]
    extTable[, pvalues := vortexR:::pval(SSMD)]

    if (save2disk == T) {
        # write results
        df2disk(extTable, dir_out, fname, ".PextTable")
        df2disk(data, dir_out, fname, ".withExt")
    }
    return(list(PextTable=extTable, datWithExt=data))
}

#' Calculates p-values from z-values
#' @param x z-values
pval <- function(x) pnorm(abs(x), lower.tail=FALSE)
