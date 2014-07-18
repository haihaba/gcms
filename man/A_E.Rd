\name{A_E}
\alias{A_E}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
A_E(DATA, Nmin, Nmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DATA}{
%%     ~~Describe \code{DATA} here~~
}
  \item{Nmin}{
%%     ~~Describe \code{Nmin} here~~
}
  \item{Nmax}{
%%     ~~Describe \code{Nmax} here~~
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
function (DATA, Nmin, Nmax) 
{
    edges <- c(1, length(DATA))
    N <- 0.5
    while (max(diff(edges)) > Nmax) {
        n <- which.max(diff(edges))
        range <- (edges[n] + Nmin):(edges[n + 1] - Nmin)
        x <- E_peak_pick(DATA[range], median(DATA[range]), Nmin)
        if (!any(x$peak > 0)) {
            N <- N * 1.5
            if (N > 7) {
                n <- which.min(DATA[range])
                if (length(n) > 1) {
                  k <- round(length(n)/2)
                  n <- n[k]
                }
                N <- 0.5
                x$peak[n] <- 1
            }
        }
        else N <- 0.5
        edges <- sort(c(edges, range[which(x$peak > 0)]))
    }
    edges <- edges[!(edges == 1 | edges == length(DATA) | edges > 
        max(which(DATA > 0)) | edges < min(which(DATA > 0)))]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
