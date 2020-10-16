% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.R
\name{ud_fit}
\alias{ud_fit}
\alias{ud_fit_control_default}
\title{Fit Empirical Bayes Multivariate Normal Means Model}
\usage{
ud_fit(X, k, fit0, w, U, S = diag(ncol(X)), control = list(), verbose = TRUE)

ud_fit_control_default()
}
\arguments{
\item{X}{The n x m data matrix, in which each row is an
m-dimensional observation.}

\item{k}{An integer 2 or greater specifying the number of
components in the mixture-of-normals prior. This only needs to be
provided if neither \code{U} nor \code{w} are provided.}

\item{fit0}{A previous Ultimate Deconvolution fit. This is useful
for "re-fitting" the model. When a value for this argument is
given, the other inputs \code{k}, \code{w}, \code{U} and \code{S}
are not needed.}

\item{w}{A numeric vector of length k giving initial estimates of
the prior mixture weights. All entries must be non-negative, but
need not sum to 1; the mixture weights are automatically normalized
to sum to 1. If not provided, the mixture weights are set to
uniform.}

\item{U}{A list of length k giving initial estimates of the
covariance matrices in the mixture-of-multivariate-normals prior;
list element \code{U[[i]]} is the m x m covariance matrix for the
ith mixture component. If not provided, the initial estimates are
randomly generated by computing sample covariances on random
subsets of the data \code{X}.}

\item{S}{An m x m matrix giving the initial estimate of the
residual covariance matrix.}

\item{control}{A list of parameters controlling the behaviour of
the model fitting algorithm. See \sQuote{Details}.}

\item{verbose}{When \code{verbose = TRUE}, information about the
algorithm's progress is printed to the console at each
iteration. For interpretation of the columns, see the description
of the \code{progress} return value.}
}
\value{
A list object with the following elements:

\item{w}{A vector containing the estimated prior mixture
  weights. When \code{control$update.w = "none"}, this will be the
  same as the \code{w} provided at input.}

\item{U}{A list containing the estimated prior covariance matrices. When
  \code{control$update.U = "none"}, this will be the same as the \code{U}
  provided at input.}

\item{S}{The estimated residual covariance matrix. When
  \code{control$update.S = "none"}, this will be the same as the \code{S}
  provided at input.}

\item{loglik}{The log-likelihood at the current settings of the
  model parameters, \code{w}, \code{U} and \code{S}.}
  
\item{progress}{A data frame containing detailed information about
  the algorithm's progress. The columns of the data frame are:
  "iter", the iteration number; "loglik", the log-likelihood at the
  current estimates of the model parameters; "delta.w", the largest
  change in the mixture weights; "delta.u", the largest change in the
  prior covariance matrices; "delta.s", the largest change in the
  residual covariance matrix; and "timing", the elapsed time in
  seconds (recorded using \code{\link{proc.time}}).}
}
\description{
This function implements "Ultimate Deconvolution", an
empirical Bayes method for fitting a multivariate normal means
model. This method is closely related to approaches for
multivariate density deconvolution (Sarkar \emph{et al}, 2018), so
it can also be viewed as a method for multivariate density
deconvolution.
}
\details{
In the multivariate normal means model, each m-dimensional
observation \eqn{x} is drawn from a mixture of multivariate
normals, \eqn{x ~ w_1 N(0, T_1) + ... + w_k N(0, T_k)}, where
\eqn{k} is the number of mixture components, the \eqn{w_j}'s are
the mixture weights, and each \eqn{T_j = S + U_j} is a covariance
matrix. This is the marginal density derived from a model in which
each \eqn{x} is multivariate normal with mean \eqn{y} and
covariance \eqn{S}, and the underlying, or "latent", signal \eqn{y}
is in turn modeled by a mixture prior in which each mixture
component \eqn{j} is multivariate normal with zero mean and
covariance matrix \eqn{S_j}. The "Extreme Deconvolution" (ED) model
(Bovy \emph{et al}, 2011) is a slight generalization of this
multivariate normal means model; the ED model allows for the
mixture components to have nonzero means, sample-specific residual
covariances, and it allows one to specify a linear projection of
the underlying signal onto the observed signal. Therefore, this
method also implements a useful special case of Extreme
Deconvolution.

The multivariate normal means model is fit by
expectation-maximization (EM). The \code{control} argument is a
list in which any of the following named components will override
the default optimization algorithm settings (as they are defined by
\code{ud_fit_control_default}):

\describe{

\item{\code{update.w}}{When \code{update.w = "em"},
maximum-likelihood estimates of the mixture weights are computed
via an EM algorithm; when \code{update.w = "mixsqp"}, the mix-SQP
solver is used to update the mixture weights by maximizing the
likelihood with the other parameters fixed; when \code{update.w =
"none"}, the mixture weights are not updated.}

\item{\code{update.U}}{This setting determines the algorithm used
to estimate the prior covariance matrices. Two EM variants are
implemented: \code{update.U = "ed"}, the EM algorithm described by
Bovy \emph{et al} (2011); and \code{update.U = "teem"} (truncated
eigenvalue expectation-maximization), in which the M-step update
for each covariance matrix \code{U[[j]]} is solved by truncating
the eigenvalues in a spectral decomposition of the unconstrained
maximimum likelihood estimate (MLE) of \code{U[j]]}. The latter
method provides greater freedom in the updates for U, and is
expected to yield better fits, and converge more quickly to a good
fit.}

\item{\code{update.S}}{When \code{update.S = "em"},
maximum-likelihood estimate of the residual covariance matrix is
computed via EM; when \code{update.S = "none", the residual
covariance parameter is not updated.}}

\item{\code{version}}{R and C++ implementations of the model
fitting algorithm are provided; these are selected with
\code{version = "R"} and \code{version = "Rcpp"}.}

\item{\code{maxiter}}{The upper limit on the number of EM updates
to perform.}

\item{\code{tol}}{Convergence tolerance for the EM algorithm; the
updates are halted when the largest change in the model parameters
between two successive iterations of EM is less than \code{tol}.}

\item{\code{minval}}{A small, non-negative number specifying the
lower bound on the eigenvalues of the prior covariance matrices
\code{U}.}}

Using this function requires some care; currently only minimal
argument checking is performed. See the documentation and examples
for guidance.
}
\examples{
library(mvtnorm)
set.seed(1)
X   <- rmvt(1000,diag(2),df = 4)
fit <- ud_fit(X,k = 10,S = diag(2))

}
\references{
J. Bovy, D. W. Hogg and S. T. Roweis (2011). Extreme Deconvolution:
inferring complete distribution functions from noisy, heterogeneous
and incomplete observations. \emph{Annals of Applied Statistics},
\bold{5}, 1657–1677. doi:10.1214/10-AOAS439

A. Sarkar, D. Pati, A. Chakraborty, B. K. Mallick and R. J. Carroll
(2018). Bayesian semiparametric multivariate density
deconvolution. \emph{Journal of the American Statistical
Association} \bold{113}, 401–416. doi:10.1080/01621459.2016.1260467

J. Won, J. Lim, S. Kim and B. Rajaratnam
(2013). Condition-number-regularized covariance estimation.
\emph{Journal of the Royal Statistical Society, Series B} \bold{75},
427–450. doi:10.1111/j.1467-9868.2012.01049.x
}