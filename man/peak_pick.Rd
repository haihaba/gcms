\name{peak_pick}
\alias{peak_pick}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
peak_pick(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x) 
{
    xout <- x <- as.matrix(x)
    xd1 <- sgolayfilt(xout, 3, 11, 1)
    NOISE_LIMIT <- median(x)
    N1 <- which(x > NOISE_LIMIT)
    N2 <- matrix(0, nrow(x), ncol(x))
    for (i in 3:(ncol(x) - 2)) if (xd1[i - 2] > 0) 
        if (xd1[i - 1] > 0) 
            if (xd1[i + 1] < 0) 
                if (xd1[i + 2] < 0) 
                  if (sum(N1 == i) == 1) 
                    N2[i] <- TRUE
    N <- which(as.logical(N2))
    if (length(N) > 1) {
        while (min(diff(N)) < 3) {
            p1 <- min(which(diff(N) < 3))
            p2 <- p1 + 1
            if (xout[N[p1]] < xout[N[p2]]) 
                N <- N[-p1]
            else N <- N[-p2]
            if (length(N) == 1) 
                break
        }
    }
    xpeak <- matrix(0, nrow(x), ncol(x))
    xpeak[N] <- xout[N]
    out <- list(xpeak = xpeak, xout = xout)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
