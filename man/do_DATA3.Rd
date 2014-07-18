\name{do_DATA3}
\alias{do_DATA3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
do_DATA3(INT, MZ, TIME, count, MZmin, MZmax, MZP)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{INT}{
%%     ~~Describe \code{INT} here~~
}
  \item{MZ}{
%%     ~~Describe \code{MZ} here~~
}
  \item{TIME}{
%%     ~~Describe \code{TIME} here~~
}
  \item{count}{
%%     ~~Describe \code{count} here~~
}
  \item{MZmin}{
%%     ~~Describe \code{MZmin} here~~
}
  \item{MZmax}{
%%     ~~Describe \code{MZmax} here~~
}
  \item{MZP}{
%%     ~~Describe \code{MZP} here~~
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
function (INT, MZ, TIME, count, MZmin, MZmax, MZP) 
{
    MZ <- round(MZ, digits = MZP) * 10^MZP
    INT <- round(INT)
    lenINT <- length(INT)
    lencount <- length(count)
    lenTIME <- length(TIME)
    DATA <- matrix(0, lenTIME, max(MZ))
    scans <- 0
    for (i in 1:lenTIME) {
        if (count[i] != 0) {
            mz <- MZ[(sum(count[0:(i - 1)]) + 1):sum(count[1:i])]
            int <- INT[(sum(count[0:(i - 1)]) + 1):sum(count[1:i])]
            DATA[i, mz] <- int
        }
    }
    return(DATA)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
