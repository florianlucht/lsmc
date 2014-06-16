\docType{package}
\name{lsmc}
\alias{lsmc}
\alias{lsmc-package}
\title{Longstaff-Schwartz-Monte-Carlo-Algorithmus}
\usage{
lsmc(X, payoffFct, interestRate, maturity, payoffParameter = NULL,
  numExerciseDates = NULL, basisFct = NULL, numBasisFct = NULL,
  regressionCoef = NULL, regressionMethod = "qr", onlyInTheMoney = TRUE,
  pathGeneratorFct = NULL, pathInfo = NULL, verbose = FALSE)
}
\arguments{
  \item{X}{matrix where each column represents a price path
  of the underlying.  In the case of multiple assets, X is
  a list with length of number of assets. Each element
  corresponds to the paths of each asset as matrix.  If
  \code{pathGeneratorFct} is specified, \code{X} (scalar)
  will be used as initial value.}

  \item{payoffFct}{payoff function with arguments \code{(S,
  ...)} where \code{S} is a vector of prices. In the case
  of multiple assets, \code{S} is a list of vectors.}

  \item{interestRate}{interest rate (scalar)}

  \item{maturity}{maturity (scalar)}

  \item{payoffParameter}{optional list of arguments for the
  payoff function (list)}

  \item{numExerciseDates}{optional number of exercise dates
  (scalar). If the value is not specified, the number of
  rows of X minus one is used.}

  \item{basisFct}{optional basis function used for
  regression (function)}

  \item{numBasisFct}{number of basis functions (scalar).
  Required if \code{basisFct} is used.}

  \item{regressionCoef}{optional precalculated regression
  coefficients (\code{numExerciseDates-1} x
  \code{numBasisFct} matrix) (leads to low biased
  valuation)}

  \item{regressionMethod}{optional string to specify the
  method for solving the least squares problem. Currently
  \code{method = "qr"} for QR decomposition (default) and
  \code{method = "svd"} for singular value decomposition
  are supported.}

  \item{onlyInTheMoney}{logical. If \code{TRUE} (default)
  only paths which are in-the-money are used.}

  \item{pathGeneratorFct}{optional function to generate
  price paths. If this is specified, \code{X} is used as
  initial value.}

  \item{pathInfo}{an optional additional path information
  needed for valuation, e.g.  average price for asian
  option (scalar or vector)}

  \item{verbose}{logical. If \code{TRUE} extra information
  will be printed.}
}
\value{
list: \describe{ \item{V}{vector of option values}
\item{regressionCoef}{regression coefficients as
\code{numExerciseDates}x\code{numBasisFct} matrix}
\item{exerciseDates}{vector of exercise dates for each
path} \item{elapsedTime}{elapsed time for calculation} }
}
\description{
Longstaff-Schwartz-Monte-Carlo-Algorithmus

\coder(lsmc) implements the least-squares monte-carlo
alogrithm of Longstaff and Schwartz to valuing american
options.
}
\examples{
T <- 1
r <- 0.05
sigma <- 0.4
initalValue <- 100

N <- 100
M <- 10000

S <- blackscholes(initialValue, r, sigma, numDates=N, maturity=T, numPaths=M)

strike <- 100
payoffParameter <- list(strike=strike)
payoff <- put

result <- lsmc(S, payoff, r, T, payoffParameter)
cat("option value: ", result$V)
}

