\name{heston}
\alias{heston}
\title{Heston model}
\usage{
heston(initialValue, r, kappa, theta, rho, xi, sigmaZero, numDates,
  maturity = 1, numPaths = 1, stepSize = FALSE, antithetic = FALSE,
  plot = FALSE)
}
\arguments{
  \item{initialValue}{initial value (scalar)}

  \item{r}{interest rate (scalar)}

  \item{kappa}{speed of adjustment (scalar)}

  \item{theta}{long run volatility (scalar)}

  \item{rho}{correlation (scalar)}

  \item{xi}{volatlility of volatility (scalar)}

  \item{sigmaZero}{initial volatility (scalar)}

  \item{numDates}{number of dates (scalar)}

  \item{maturity}{maturity (scalar)}

  \item{numPaths}{number of paths (scalar)}

  \item{antithetic}{logical. If \code{TRUE} antitethic
  paths are generated.}

  \item{plot}{logical. If \code{TRUE} paths are plotted.}
}
\value{
simulated paths as matrix with \code{numDates+1} rows and
\code{numPaths} columns
}
\description{
This method simulates discrete paths in the Heston model.
The price paths satisfy the following stochastic
differential equations: \deqn{dS_t = S_t \left( r dt +
\sqrt{\nu_t} d W_t^{S} \right) \ ,} {dS = S * ( r * dt +
sqrt{nu} * dW ) ,} \deqn{d\nu_t = \kappa ( \theta - \nu_t )
dt + \xi \sqrt{\nu_t} dW_t^{\nu} \ ,} {d nu = kappa * (
theta - nu ) * dt + xi * sqrt{nu_t} * dB ,} \deqn{\rho =
Kor\left(W_t^{S}, W_t^{\nu}\right) \ .}{rho = Kor[ W, B ]
.}
}
\examples{
# generate one simple path
S <- heston(100, 0.05, 2, 0.4, -0.6, 0.4, 0.4, 100)

# generate 100 paths with maturity T=10
S <- heston(100, 0.05, 2, 0.4, -0.6, 0.4, 0.4, 100, maturity=10, numPaths=100)
}
\seealso{
Other models: \code{\link{blackscholes}};
\code{\link{merton}}; \code{\link{nig}}
}

