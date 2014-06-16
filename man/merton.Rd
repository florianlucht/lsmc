\name{merton}
\alias{merton}
\title{Merton model}
\usage{
merton(initialValue, r, sigma, lambda, mu, delta, numDates, maturity = 1,
  numPaths = 1, antithetic = FALSE, plot = FALSE)
}
\arguments{
  \item{initialValue}{initial value (scalar)}

  \item{r}{interest rate (scalar)}

  \item{sigma}{volatility (scalar)}

  \item{lambda}{intensity/rate of the jumps (scalar)}

  \item{mu}{mean of the jumps (scalar)}

  \item{delta}{variance of the jumps (scalar)}

  \item{numDates}{number of dates (scalar)}

  \item{maturity}{maturity (scalar)}

  \item{numPaths}{number of paths (scalar)}

  \item{antitheric}{logical. If \code{TRUE} antitethic
  paths are generated.}

  \item{plot}{logical. If \code{TRUE} paths are plotted.}
}
\value{
simulated paths as matrix with \code{numDates+1} rows and
\code{numPaths} columns
}
\description{
This method simulates discrete paths in the Merton model.
The price paths satisfy the following stochastic
differential equations: \deqn{S_t = S_0 \exp \left[ \left(
\mu-\frac{\sigma^2}{2} \right) t + \sigma W_t +
\sum\limits_{i=1}^{N_t} Y_i \right] \ ,} {S_t = S_0 exp[ (
mu- 1/2 * sigma^2 ) * t + sigma * W_t + sum[ Y_i | i =
1,...,N ] ] } \deqn{N_t \sim Pois(\lambda t) \ , \ Y_i \sim
N(\mu, \delta) \ . } {N ~ Pois(lambda*t) , Y_i ~ N(mu,
delta) . }
}
\examples{
# generate one simple path
S <- merton(100, 0.05, 0.4, 1, 0, 0.4, 100)

# generate 100 paths with maturity T=10
S <- merton(100, 0.05, 0.4, 1, 0, 0.4, 100, maturity=10, numPaths=100)
}
\seealso{
Other models: \code{\link{blackscholes}};
\code{\link{heston}}; \code{\link{nig}}
}

