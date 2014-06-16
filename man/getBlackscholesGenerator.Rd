\name{getBlackscholesGenerator}
\alias{getBlackscholesGenerator}
\title{Generator function of the Black-Scholes model}
\usage{
getBlackscholesGenerator(r, sigma, antithetic = FALSE)
}
\arguments{
  \item{r}{interest rate (scalar)}

  \item{sigma}{volatility (scalar)}

  \item{antithetic}{logical. If \code{TRUE} antitethic
  paths are generated.}
}
\value{
function with arguments: \describe{ \item{S}{current
value(s) (scalar or vector)} \item{dt}{step size (scalar)}
\item{notStopped}{optional index of paths which are not
exercised (vector of logical). This argument is used
internally to calculate the correct values if
\code{antithetic=TRUE} is given.} }
}
\description{
This method returns a function which can be used to
simulate one step in the Black-Scholes model for given
parameters.
}
\examples{
# get generator function
genBlackscholes <- getBlackScholesGenerator(0.05, 0.4)
# generate one step in the BS-model with step size 0.1 and current value 100
next <- genBlackscholes(100, 0.1)
# generate path with 100 steps
S <- c(100, rep(NaN, 100))
for (i in 2:100) {
  S[i] <- genBlackscholes(S[i-1], 0.1)
}
# blackscholes(100, 0.05, 0.4, 100) returns the same output but faster!
}

