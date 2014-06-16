\name{getNigGenerator}
\alias{getNigGenerator}
\title{Generator function of the normal-inverse Gaussian model}
\usage{
getNigGenerator(r, alpha, beta, delta)
}
\arguments{
  \item{r}{interest rate (scalar)}

  \item{alpha}{tail heavyness (scalar)}

  \item{beta}{asymmetry parameter (scalar)}

  \item{delta}{scale parameter (scalar)}
}
\value{
function with arguments: \describe{ \item{S}{current
value(s) (scalar or vector)} \item{dt}{step size (scalar)}
}
}
\description{
This method returns a function which can be used to
simulate one step in the normal-inverse Gaussian model for
given parameters.
}
\examples{
# get generator function
genNig <- getNigGenerator(0.05, 20, -5, 1)
# generate one step in the NIG-model with step size 0.1 and current value 100
next <- genNig(100, 0.1)
# generate path with 100 steps
S <- c(100, rep(NaN, 100))
for (i in 2:100) {
  S[i] <- genNig(S[i-1], 0.1)
}
# nig(100, 0.05, 20, -5, 1, 100) returns the same output but faster!
}

