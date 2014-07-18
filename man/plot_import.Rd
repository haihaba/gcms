\name{plot_import}
\alias{plot_import}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_import(projectpath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{projectpath}{
%%     ~~Describe \code{projectpath} here~~
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
function (projectpath) 
{
    require(tcltk)
    if (file.exists(file.path(projectpath, "Filtered"))) 
        setwd(file.path(projectpath, "Filtered"))
    rdatafiles <- tk_choose.files(caption = "Select imported Rdata files", 
        multi = TRUE)
    setwd(projectpath)
    load(rdatafiles[1])
    plot(apply(Xbc, 1, sum), type = "l", ylab = "signal", xlab = "scans")
    rm(list = c("Xbc", "SCAN_INFO", "SCAN_RANGE"))
    if (length(rdatafiles) > 1) {
        for (i in 2:length(rdatafiles)) {
            load(rdatafiles[i])
            points(apply(Xbc, 1, sum), type = "l")
            rm(list = c("Xbc", "SCAN_INFO", "SCAN_RANGE"))
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
