% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load.R
\name{simulation}
\alias{simulation}
\title{Semiparametric Causal Mediation Modeling of Semi-Competing Risks}
\usage{
simulation(
  simulation_type,
  hypo,
  sample_size = 1000,
  repeat_size = 1000,
  num_of_cores = 1,
  timer = TRUE,
  get_variance = c("a", "b")
)
}
\arguments{
\item{simulation_type}{Choices are \code{1} (unbiasedness) or \code{2} (coverage).}

\item{hypo}{Choices are \code{null} or \code{alter}.}

\item{sample_size}{a number only effective when \code{simulation_type = 2}. Default is \code{1000}.}

\item{repeat_size}{Default is \code{1000}.}

\item{num_of_cores}{the number of cores assigned. Default is \code{1}.}

\item{timer}{will show the progress. Default is \code{TRUE}.}
}
\value{
\code{simulation} returns 6 plots if \code{simulation_type == 1}; a data frame containing coverage if \code{simulation_type == 2}.
}
\description{
A function that simulates different situations.
}
