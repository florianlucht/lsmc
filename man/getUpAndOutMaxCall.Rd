\name{getUpAndOutMaxCall}
\alias{getUpAndOutMaxCall}
\title{Generate Up-And-Out Max-Call option}
\usage{
getUpAndOutMaxCall(strike, barrier, forward = FALSE)
}
\arguments{
  \item{strike}{sktrike (scalar)}

  \item{barrier}{barrier (scalar)}

  \item{forward}{logical. If \code{TRUE} the function keeps
  the knock-out information.}
}
\value{
Methode mit folgenden Arguementen: \describe{
\item{x}{Basiswert(e) (Vektor oder einzelner Wert)}
\item{para}{TODO} \item{notStopped}{Index der
nichtausgeuebten Pfade (siehe lsmc Funktion)}
\item{reset}{falls TRUE wird beim naechsten Aufruf die
Knocked-out Information verworfen} }
}
\description{
This method returns a function which represents a
Up-And-Out Max-Call option with payoff \deqn{(K - \min_{i=
1,\dots,d} x_i)^+}{max[(K - min[x_i]), 0]}
}
\seealso{
Other payoff: \code{\link{bearspread}};
\code{\link{bullspread}}; \code{\link{call}};
\code{\link{getBearspread}}; \code{\link{getBullspread}};
\code{\link{getCall}}; \code{\link{getMinPut}};
\code{\link{getPut}}; \code{\link{minPut}};
\code{\link{put}}; \code{\link{upAndOutMaxCall}}
}

