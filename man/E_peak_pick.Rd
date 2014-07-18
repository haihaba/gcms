\name{E_peak_pick}
\alias{E_peak_pick}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
E_peak_pick(x, NL, min_win)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{NL}{
%%     ~~Describe \code{NL} here~~
}
  \item{min_win}{
%%     ~~Describe \code{min_win} here~~
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
function (x, NL, min_win) 
{
    require(signal)
    xd1 <- sgolayfilt(x, 3, 11, 1)
    N1 <- which(x < NL)
    N2 <- numeric(length(x))
    for (i in 3:(length(x) - 2)) {
        if (xd1[i - 2] < 0 & xd1[i - 1] < 0 & xd1[i + 1] > 0 & 
            xd1[i + 2] > 0 & sum(N1 == i) == 1) 
            N2[i] <- 1
    }
    N <- intersect(N1, which(N2 == 1))
    if (length(N) == 0) 
        return()
    if (length(N) != 1) 
        while (min(diff(N)) < min_win) {
            p1 <- which.min(diff(N))
            p2 <- p1 + 1
            if (x[N[p1]] < x[N[p2]]) 
                N <- N[-p1]
            else N <- N[-p2]
            if (length(N) == 1) 
                break
        }
    xpeak <- numeric(length(x))
    xpeak[N] <- x[N]
    x <- list(peak = xpeak, out = x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
