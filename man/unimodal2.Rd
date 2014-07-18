\name{unimodal2}
\alias{unimodal2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
unimodal2(C, R, var, obs, RT_LIMIT)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{C}{
%%     ~~Describe \code{C} here~~
}
  \item{R}{
%%     ~~Describe \code{R} here~~
}
  \item{var}{
%%     ~~Describe \code{var} here~~
}
  \item{obs}{
%%     ~~Describe \code{obs} here~~
}
  \item{RT_LIMIT}{
%%     ~~Describe \code{RT_LIMIT} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (C, R, var, obs, RT_LIMIT) 
{
    N <- nrow(C)
    CP <- numeric(R)
    for (i in 1:R) {
        c <- matrix(C[, i], var, obs)
        out <- unimodal3(c, obs, RT_LIMIT)
        C[, i] <- matrix(out$C, N, 1)
        CP[i] <- out$cp
    }
    out <- list(C = C, CP = CP)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
