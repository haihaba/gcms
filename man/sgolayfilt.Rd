\name{sgolayfilt}
\alias{sgolayfilt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sgolayfilt(x, p = 3, n = p + 3 - p%%2, m = 0, ts = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
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
function (x, p = 3, n = p + 3 - p\%\%2, m = 0, ts = 1) 
{
    len = length(x)
    if (class(p) == "sgolayFilter" || (!is.null(dim(p)) && dim(p) > 
        1)) {
        F = p
        n = nrow(F)
    }
    else F = sgolay(p, n, m, ts)
    k = floor(n/2)
    filt <- F[k + 1, n:1]
    z <- na.omit(stats:::filter(c(rep(0, length(filt) - 1), x), 
        filt, sides = 1))
    c(F[1:k, ] \%*\% x[1:n], z[n:len], F[(k + 2):n, ] \%*\% x[(len - 
        n + 1):len])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
