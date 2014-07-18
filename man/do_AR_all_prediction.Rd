\name{do_AR_all_prediction}
\alias{do_AR_all_prediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
do_AR_all_prediction(x, S, CP, RT_LIMIT)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{CP}{
%%     ~~Describe \code{CP} here~~
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
function (x, S, CP, RT_LIMIT) 
{
    var <- nrow(x)
    mz <- ncol(x)
    C <- x \%*\% S \%*\% ginv(t(S) \%*\% S)
    C[C < 0] <- 0
    Cnew <- C * 0
    for (i in 1:ncol(C)) {
        c <- as.matrix(C[, i], nrow = var)
        Cnew[, i] <- unimodal_3(c, 1, CP[i], RT_LIMIT)
    }
    C <- Cnew
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
