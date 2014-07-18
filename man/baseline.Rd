\name{baseline}
\alias{baseline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
baseline(X, projectpath, FL = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{projectpath}{
%%     ~~Describe \code{projectpath} here~~
}
  \item{FL}{
%%     ~~Describe \code{FL} here~~
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
function (X, projectpath, FL = 0) 
{
    N <- nrow(X)
    K <- ncol(X)
    if (FL) {
        Xbc <- X * 0
        Xbc[1:FL, ] <- X[1:FL, ]
        Xbc[(nrow(X) - FL + 1):nrow(X), ] <- X[(nrow(X) - FL + 
            1):nrow(X), ]
        start <- FL + 1
        end <- dim(X)[1] - FL
        Xbc[start:end, ] <- filter(X, filter = rep(1/(FL * 2 + 
            1), (FL * 2 + 1)))[start:end, ]
        if (any(Xbc < 0)) 
            Xbc[Xbc < 0] <- 0
    }
    else Xbc <- X
    return(Xbc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
