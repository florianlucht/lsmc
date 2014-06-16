\name{upAndOutMaxCall}
\alias{upAndOutMaxCall}
\title{Up-And-Out Max-Call option}
\usage{
upAndOutMaxCall(x, strike, barrier, path)
}
\arguments{
  \item{x}{market value(s) of the underlying(s) (scalar or
  vector)}

  \item{strike}{strike (scalar)}

  \item{barrier}{barrier (scalar)}

  \item{path}{previous path of the underlying asset's price
  (vector)}
}
\value{
Optionswert(e)
}
\description{
Option with payoff \deqn{(\max\{S_t^{(1)},\dots,S_t^{(d)}\}
- K )^+ \cdot 1_{\max\{S_u^{(d)} \mid 0\leq u \leq t \ , \
1 \leq l \leq d \} \leq B }}{ max{ max{ S^1 , ... , S^d } -
K , 0 } * ( S_u^n < B) fuer alle 0 <= u <=t und alle S.}
}
\seealso{
Other payoff: \code{\link{bearspread}};
\code{\link{bullspread}}; \code{\link{call}};
\code{\link{getBearspread}}; \code{\link{getBullspread}};
\code{\link{getCall}}; \code{\link{getMinPut}};
\code{\link{getPut}}; \code{\link{getUpAndOutMaxCall}};
\code{\link{minPut}}; \code{\link{put}}
}

