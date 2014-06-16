\name{getHestonGenerator}
\alias{getHestonGenerator}
\title{Generator function of the Heston model}
\usage{
getHestonGenerator(r, kappa, theta, rho, xi, sigmaZero, antithetic = FALSE,
  stepSize = FALSE)
}
\arguments{
  \item{r}{interest rate (scalar)}

  \item{kappa}{speed of adjustment (scalar)}

  \item{theta}{long run volatility (scalar)}

  \item{rho}{correlation (scalar)}

  \item{xi}{volatlility of volatility (scalar)}

  \item{sigmaZero}{initial volatility (scalar)}

  \item{antithetic}{logical. If \code{TRUE} antitethic
  paths are generated.}

  \item{stepSize}{step size in the Euler–Maruyama method
  for process simulation}
}
\value{
function with arguments: \describe{ \item{S}{current
value(s) (scalar or vector)} \item{dt}{step size (scalar)}
\item{notStopped}{optional index of paths which are not
exercised (vector of logical). This argument is used
internally to calculate the correct values if
\code{antithetic=TRUE} is given.} \item{reset}{logical. If
\code{TRUE} the volatility process will be reset to the
initial value on the next function call.}
}
\description{
This method returns a function which can be used to
simulate one step in the Heston model for given parameters.
}
\examples{
# get generator function
genHeston <- getHestonGenerator(0.05, 2, 0.4, -0.6, 0.4, 0.4)
# generate one step in the Heston model with step size 0.1 and current value 100
next <- genHeston(100, 0.1)
# reset volatility process
genHeston(NULL, NULL, reset=TRUE)
# genrate path with 100 steps
S <- c(100, rep(NaN, 100))
for (i in 2:100) {
  S[i] <- genHeston(S[i-1], 0.1)
}
# heston(100, 0.05, 2, 0.4, -0.6, 0.4, 0.4, 100) returns the same output but faster!
}

