% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dentist.R
\name{plot.dentist}
\alias{plot.dentist}
\title{Plot the dented samples
This will show the univariate plots of the parameter values versus the likelihood as well as bivariate plots of pairs of parameters to look for ridges.}
\usage{
\method{plot}{dentist}(x, local.only = FALSE, ...)
}
\arguments{
\item{x}{An object of class dentist}

\item{local.only}{Boolean indicating whether to trim x and y lims to be near accepted points}

\item{...}{Other arguments to pass to plot}
}
\description{
Plot the dented samples
This will show the univariate plots of the parameter values versus the likelihood as well as bivariate plots of pairs of parameters to look for ridges.
}
