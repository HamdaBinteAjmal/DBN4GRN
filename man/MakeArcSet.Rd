% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{MakeArcSet}
\alias{MakeArcSet}
\title{Take a learned Bayesian network structure (bnlearn object) and
returns a data frame in which each row represents a pair of nodes conected in the network.
This function is used to retrieve the arcs from the network.}
\usage{
MakeArcSet(bn, data)
}
\arguments{
\item{bn}{Bnlearn object Bayesian network}

\item{data}{data.frame from which the Bayesian network structure is learned.}
}
\value{
arc_set a data.frame with two columns. First column shows the regulators, second column shows the targets.
}
\description{
Take a learned Bayesian network structure (bnlearn object) and
returns a data frame in which each row represents a pair of nodes conected in the network.
This function is used to retrieve the arcs from the network.
}
\keyword{internal}
