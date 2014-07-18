\name{pure}
\alias{pure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pure(d, nr, f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{nr}{
%%     ~~Describe \code{nr} here~~
}
  \item{f}{
%%     ~~Describe \code{f} here~~
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
function (d, nr, f) 
{
    w <- p <- s <- matrix(0, nrow = max(nr, 1), ncol = ncol(d))
    f <- f/100
    s[1, ] <- apply(d, 2, sd)
    m <- colMeans(d)
    ll <- s[1, ]^2 + m^2
    f <- max(m) * f
    p[1, ] <- s[1, ]/(m + f)
    imp <- numeric()
    mp <- max(p[1, ])
    imp[1] <- which.max(p[1, ])
    l <- sqrt(s[1, ]^2 + (m + f)^2)
    for (i in 1:nrow(d)) d[i, ] <- d[i, ]/l
    c <- (t(d) \%*\% d)/nrow(d)
    w[1, ] <- ll/l^2
    p[1, ] <- w[1, ] * p[1, ]
    s[1, ] <- w[1, ] * s[1, ]
    if (nr > 1) {
        for (i in 2:nr) {
            for (j in 1:ncol(d)) {
                dm <- wmat(c, imp, i, j)
                w[i, j] <- det(dm)
                p[i, j] <- p[1, j] * w[i, j]
                s[i, j] <- s[1, j] * w[i, j]
            }
            mp[i] <- max(p[i, ])
            imp[i] <- which.max(p[i, ])
        }
    }
    sp <- normv2(t(d[, imp[1:nr]]))
    purelist <- list(sp = sp, imp = imp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
