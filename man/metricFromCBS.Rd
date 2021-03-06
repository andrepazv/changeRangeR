% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fdMetricCommMetafunc.R
\name{metricFromCBS}
\alias{metricFromCBS}
\title{title}
\usage{
metricFromCBS(
  cbsDir,
  scenario,
  env,
  spIndTable,
  cellAttributeTable = NULL,
  speciesAttributeTable = NULL,
  mc.cores = 1,
  FUN = phyloDiv,
  outputFUNnames,
  verbose = F,
  ...
)
}
\arguments{
\item{FUN}{function defining the metric to be calculated.}

\item{outputFUNnames}{desired output name of the metric}

\item{...}{other arguments to be passed to FUN}
}
\value{

}
\description{
title
}
\note{


most functions `FUN` will need to have the colnames to be the name of the species
}
\author{
Cory Merow <cory.merow@gmail.com>
}
