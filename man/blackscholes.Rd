\name{blackscholes}
\alias{blackscholes}
\title{Black-Scholes model}
\usage{
blackscholes(initialValue, r, sigma, numDates, maturity = 1, numPaths = 1,
  antithetic = FALSE, plot = FALSE)
}
\arguments{
  \item{initialValue}{initial value (scalar)}

  \item{r}{interest rate (scalar)}

  \item{sigma}{volatility (scalar)}

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
This method simulates discrete paths in the Black-Scholes
model. The price paths follow a geometric Brownian motion,
hence they satisfy the following stochastic differential
equation \deqn{dS = S ( r dt + \sigma dW )}{S = S * (r * dt
+ sigma * dW) } i.e. its explicit solution \deqn{S_t = S_0
\exp\left[ \left(r - \frac{1}{2}\sigma^2 \right)t + \sigma
W_t \right] } {S_t = S_0 exp[ (r - 1/2 * sigma^2) * t +
sigma * W_t ] .}
}
\examples{
# generate one simple path
S <- blackscholes(100, 0.05, 0.4, 100)

# generate 100 paths with maturity T=10
S <- blackscholes(100, 0.05, 0.4, 100, maturity=10, numPaths=100)
}
\seealso{
Other models: \code{\link{heston}}; \code{\link{merton}};
\code{\link{nig}}
}

