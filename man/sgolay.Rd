\name{sgolay}
\alias{sgolay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sgolay(p, n, m = 0, ts = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{ts}{
%%     ~~Describe \code{ts} here~~
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
function (p, n, m = 0, ts = 1) 
{
    library(MASS)
    if (n\%\%2 != 1) 
        stop("sgolay needs an odd filter length n")
    if (p >= n) 
        stop("sgolay needs filter length n larger than polynomial order p")
    F = matrix(0, n, n)
    k = floor(n/2)
    for (row in 1:(k + 1)) {
        C = (((1:n) - row) \%*\% matrix(1, 1, p + 1))^(matrix(1, 
            n) \%*\% (0:p))
        A = ginv(C)
        F[row, ] = A[1 + m, ]
    }
    F[(k + 2):n, ] = (-1)^m * F[k:1, n:1]
    if (m > 0) 
        F = F * prod(1:m)/(ts^m)
    class(F) = "sgolayFilter"
    F
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
