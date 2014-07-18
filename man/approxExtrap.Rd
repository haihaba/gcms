\name{approxExtrap}
\alias{approxExtrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
approxExtrap(x, y, xout, method = "linear", n = 50, rule = 2, f = 0, ties = "ordered", na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{xout}{
%%     ~~Describe \code{xout} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{rule}{
%%     ~~Describe \code{rule} here~~
}
  \item{f}{
%%     ~~Describe \code{f} here~~
}
  \item{ties}{
%%     ~~Describe \code{ties} here~~
}
  \item{na.rm}{
%%     ~~Describe \code{na.rm} here~~
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
function (x, y, xout, method = "linear", n = 50, rule = 2, f = 0, 
    ties = "ordered", na.rm = FALSE) 
{
    if (is.list(x)) {
        y <- x[[2]]
        x <- x[[1]]
    }
    if (na.rm) {
        d <- !is.na(x + y)
        x <- x[d]
        y <- y[d]
    }
    d <- !duplicated(x)
    x <- x[d]
    y <- y[d]
    d <- order(x)
    x <- x[d]
    y <- y[d]
    w <- approx(x, y, xout = xout, method = method, n = n, rule = 2, 
        f = f, ties = ties)$y
    r <- range(x)
    d <- xout < r[1]
    if (any(is.na(d))) 
        stop("NAs not allowed in xout")
    if (any(d)) 
        w[d] <- (y[2] - y[1])/(x[2] - x[1]) * (xout[d] - x[1]) + 
            y[1]
    d <- xout > r[2]
    n <- length(y)
    if (any(d)) 
        w[d] <- (y[n] - y[n - 1])/(x[n] - x[n - 1]) * (xout[d] - 
            x[n - 1]) + y[n - 1]
    list(x = xout, y = w)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
