\name{do_DATA1}
\alias{do_DATA1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
do_DATA1(INT, MZ, TIME, count, MZmin, MZmax, MZP)
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
    if (length(abs(diff(MZ))[abs(diff(MZ)) != 1 & abs(diff(MZ)) != 
        (MZmax - MZmin)])) {
        cat("Using alternative method to arrange data matrix (slower)\n")
        DATA <- matrix(0, lenTIME, max(MZ))
        scans <- min(which(count != 0))
        for (i in 1:lenINT) {
            DATA[scans, MZ[i]] <- DATA[scans, MZ[i]] + INT[i]
            if (i != lenINT) 
                if (MZ[i + 1] < MZ[i]) 
                  scans <- scans + min(which(count[(scans + 1):lencount] != 
                    0))
        }
    }
    else DATA <- cbind(matrix(0, nrow = lenTIME, ncol = max((min(MZ) - 
        1), 0)), matrix(INT, nrow = lenTIME, ncol = (max(MZ) - 
        min(MZ) + 1), byrow = TRUE))
    return(DATA)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
