\name{nig}
\alias{nig}
\title{Normal-inverse Gaussian model}
\usage{
nig(initialValue, r, alpha, beta, delta, numDates, maturity = 1,
  numPaths = 1, plot = FALSE)
}
\arguments{
  \item{initialValue}{initial value (scalar)}

  \item{r}{interest rate (scalar)}

  \item{alpha}{tail heavyness (scalar)}

  \item{beta}{asymmetry parameter (scalar)}

  \item{delta}{scale parameter (scalar)}

  \item{numDates}{number of dates (scalar)}

  \item{maturity}{maturity (scalar)}

  \item{numPaths}{number of paths (scalar)}

  \item{antithetic}{logical. If \code{TRUE} antitethic
  paths are generated.}

  \item{plot}{logical. If \code{TRUE} paths are plotted.}
}
\value{
Pfade als Matrix mit \eqn{numDates + 1}{numDates + 1}
Zeilen und \eqn{numPaths} Spalten
}
\description{
This method simulates discrete paths in the normal-inverse
Gaussian model. Dieses ist beschrieben durch \deqn{ S_t =
S_0 \exp\left[ (r-\omega) + \beta I_t + W_{I_t} \right] \ ,
}{S_t = S_0 exp[ (r - omega) + beta * I_t + W_{I_t} ] ,}
\deqn{ I_t \sim IG(\delta t, \gamma) \ , \ \omega = \lambda
\gamma - \lambda \sqrt{\alpha^2 (1+\beta)^2} \ .}{I_t ~
IG(delta*t , gamma) , omega = lambda * gamma - lambda *
sqrt{ alpha^2 * (1 + beta)^2 } .}
}
\examples{
# generate one simple path
S <- nig(100, 0.05, 20, -5, 1, 100)

# generate 100 paths with maturity T=10
S <- nig(100, 0.05, 20, -5, 1, 100, maturity=10, numPaths=100)
}
\seealso{
Other models: \code{\link{blackscholes}};
\code{\link{heston}}; \code{\link{merton}}
}

