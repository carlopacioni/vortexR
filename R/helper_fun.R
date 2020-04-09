
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
#' @param row_names Whether to include row names inthe csv file
#'
#' @examples
#' my.df <- data.frame(1, 1:10, sample(LETTERS[1:3], 10, replace = TRUE))
#' my.folder <- file.path(tempdir(check=TRUE), 'test')
#' df2disk(df=my.df, dirpath=tempdir(check=TRUE), fname='testname')
#' df2disk(df=my.df, dirpath=my.folder, fname='testname', postfix='_testpostfix')
#' @export
df2disk <- function(df, dirpath, fname, postfix = "", row_names=FALSE) {

    dir.create(dirpath, showWarnings = FALSE, recursive = TRUE)

    save(df, file = paste0(dirpath, "/", fname, postfix, ".rda"), compress = "xz")

    write.table(df,
                file = paste0(dirpath, "/", fname, postfix, ".txt"), sep = ";",
                row.names = row_names)
}

#' Calculates p-values from z-values
#'
#' \code{pval} calculates one-tailed p values from a vector that contains
#'   z-values and it is geenrally used internally.
#'
#' @param x z-values
#' @return A numeric vector of length equal to length(x)
#' @export
#' @examples
#'  z <- c(1.645, 1.96, 3.09)
#'  pval(z)
pval <- function(x) pnorm(abs(x), lower.tail = FALSE)

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
#' @param fname Name of file name root
#' @importFrom gtools mixedsort
#' @return A character vector of fully qualified file paths
get_file_paths <- function(path, pattern, fn_name, fname, verbose = FALSE) {
    files <- gtools::mixedsort(list.files(path = path, pattern = pattern, full.names = TRUE))

    if (length(files) == 0) {
        stop(paste0("ERROR vortexR::", fn_name, " found no files",
                    " containing '", fname, "' in ", path))
    } else {
        if (verbose) {
            msg <- paste0("INFO vortexR::", fn_name, " found ", length(files),
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
se2sd <- function(se, no) {
    se * sqrt(no)
}


#' Return a prefixed and repeated string of character
#'
#' @param chars A string of characters (popvalue)
#' @param times The number of repetitions (ncolpop), default: 1
#' @param prefix A text prefix, default: ''
PrefixAndRepeat <- function(chars, times = 1, prefix = "") {
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
CompileIter <- function(iter, filename, n_rows, iter_ln, lines, header) {
    temp.df <- read.table(filename, header = FALSE, sep = ";", nrows = n_rows,
                          skip = iter_ln[iter], colClasses = "numeric",
                          comment.char = "")
    colnames(temp.df) <- header
    Iteration <- rep(iter, length = length(temp.df$Year))
    ScenNameStarts <- attr(
        regexpr(pattern = "\\$Scenario: ", lines[1]), "match.length") + 1
    ScenName <- substr(lines[1], ScenNameStarts, nchar(lines[1]))
    Scenario <- rep(ScenName, length = length(temp.df$Year))
    temp.df <- cbind(Scenario, Iteration, temp.df,  stringsAsFactors = TRUE)
    return(temp.df)
}
