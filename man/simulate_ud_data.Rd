% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasim.R
\name{simulate_ud_data}
\alias{simulate_ud_data}
\title{Simulate Data from Ultimate Deconvolution Model}
\usage{
simulate_ud_data(n, w, U, V)
}
\arguments{
\item{n}{Number of data points to simulate.}

\item{w}{A numeric vector of length k specifying the prior mixture
weights. All entries must be non-negative, but need not sum to 1;
the mixture weights are automatically normalized to sum to 1. If
not provided, all the mixture weights are set to the same value.}

\item{U}{A list of length k specifying the covariance matrices in
the mixture prior; list element \code{U[[i]]} is the m x m
covariance matrix for the ith mixture component. For k = 1, U may
also be a matrix.}

\item{V}{The m x m residual covariance matrix. If missing, \code{V}
is set to the identity matrix.}
}
\value{
An n x m matrix in which each row is a draw from the
  Ultimate Deconvolution model.
}
\description{
Simulate data points from the Ultimate Deconvolution
  model. See \code{\link{ud_fit}} for the model definition.
}
\seealso{
\code{\link{ud_fit}}
}
