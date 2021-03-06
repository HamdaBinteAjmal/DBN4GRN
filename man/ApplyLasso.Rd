% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{ApplyLasso}
\alias{ApplyLasso}
\title{This function applies L1 regularised Lasso regression on the dataset. The regularisation parameter is selected by cross validation.
It uses function from package "lars".}
\usage{
ApplyLasso(data, pred_pos = NULL)
}
\arguments{
\item{data}{a data.frame or a matrix with n rows (=time points) and p columns (=genes) containing the gene expression time series.}

\item{pred_pos}{indices of the potential predictor variables in the dataset. If the value is NULL, all variables are considered as potential predictors.}
}
\value{
a matrix of Lasso regression coefficients.
}
\description{
This function applies L1 regularised Lasso regression on the dataset. The regularisation parameter is selected by cross validation.
It uses function from package "lars".
}
