\name{getMertonGenerator}
\alias{getMertonGenerator}
\title{Generator function of the Merton model}
\usage{
getMertonGenerator(r, sigma, lambda, mu, delta, antithetic = FALSE)
}
\arguments{
  \item{r}{interest rate (scalar)}

  \item{sigma}{volatility (scalar)}

  \item{lambda}{intensity/rate of the jumps (scalar)}

  \item{mu}{mean of the jumps (scalar)}

  \item{delta}{variance of the jumps (scalar)}

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
simulate one step in the Merton model for given parameters.
}
\examples{
# get generator function
genMerton <- getMertonGenerator(0.05, 0.4, 1, 0, 0.4)
# generate one step in the Merton model with step size 0.1 and current value 100
next <- genMerton(100, 0.1)
# reset volatility process
genMerton(NULL, NULL, reset=TRUE)
# generate path with 100 steps
S <- c(100, rep(NaN, 100))
for (i in 2:100) {
  S[i] <- genMerton(S[i-1], 0.1)
}
# merton(100,0.05, 0.4, 1, 0, 0.4, 100) returns the same output but faster!
}

