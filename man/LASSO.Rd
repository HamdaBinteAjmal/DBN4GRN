% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{LASSO}
\alias{LASSO}
\title{Internal function. Applies Lasso regression on the dataset using the package lars}
\usage{
LASSO(X, Y, use.Gram = TRUE)
}
\arguments{
\item{X}{Regressors}

\item{Y}{Response variable}
}
\value{
Lasso regression coefficients for each outcome variable Y
}
\description{
Internal function. Applies Lasso regression on the dataset using the package lars
}
\keyword{internal}
