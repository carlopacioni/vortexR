% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_fun.R
\name{df2disk}
\alias{df2disk}
\title{Save a data.frame as both Rdata and CSV}
\usage{
df2disk(df, dirpath, fname, postfix = "", row_names = FALSE)
}
\arguments{
\item{df}{A data.frame}

\item{dirpath}{The destination path for written files, will be created if
necessary}

\item{fname}{The file name}

\item{postfix}{An optional name postfix}

\item{row_names}{Whether to include row names inthe csv file}
}
\description{
\code{df2disk} saves to disk a given data.frame as both Rdata and CSV with a
given name and optional name postfix to a given location.
}
\details{
\code{df2disk} is used by the \code{collate_} functions when the operator
chooses to save2disk.
}
\examples{
my.df <- data.frame(1, 1:10, sample(LETTERS[1:3], 10, replace = TRUE))
my.folder <- file.path(tempdir(check=TRUE), 'test')
df2disk(df=my.df, dirpath=tempdir(check=TRUE), fname='testname')
df2disk(df=my.df, dirpath=my.folder, fname='testname', postfix='_testpostfix')
}
