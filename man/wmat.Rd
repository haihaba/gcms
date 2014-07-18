\name{wmat}
\alias{wmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wmat(tempmat, imp, irank, jvar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tempmat}{
%%     ~~Describe \code{tempmat} here~~
}
  \item{imp}{
%%     ~~Describe \code{imp} here~~
}
  \item{irank}{
%%     ~~Describe \code{irank} here~~
}
  \item{jvar}{
%%     ~~Describe \code{jvar} here~~
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
function (tempmat, imp, irank, jvar) 
{
    dm <- matrix(0, irank, irank)
    dm[1, 1] <- tempmat[jvar, jvar]
    dm[1, 2:irank] <- tempmat[jvar, imp[1:(irank - 1)]]
    dm[2:irank, 1] <- tempmat[imp[1:(irank - 1)], jvar]
    dm[2:irank, 2:irank] <- tempmat[imp[1:(irank - 1)], imp[1:(irank - 
        1)]]
    dm
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
