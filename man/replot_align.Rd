\name{replot_align}
\alias{replot_align}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
replot_align(DATA, ADATA, start, stop, zoomwidth = c(0, 0), datasource, maintext, targetfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DATA}{
%%     ~~Describe \code{DATA} here~~
}
  \item{ADATA}{
%%     ~~Describe \code{ADATA} here~~
}
  \item{start}{
%%     ~~Describe \code{start} here~~
}
  \item{stop}{
%%     ~~Describe \code{stop} here~~
}
  \item{zoomwidth}{
%%     ~~Describe \code{zoomwidth} here~~
}
  \item{datasource}{
%%     ~~Describe \code{datasource} here~~
}
  \item{maintext}{
%%     ~~Describe \code{maintext} here~~
}
  \item{targetfile}{
%%     ~~Describe \code{targetfile} here~~
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
function (DATA, ADATA, start, stop, zoomwidth = c(0, 0), datasource, 
    maintext, targetfile) 
{
    zoomwidth[1] <- max(0, zoomwidth[1])
    zoomwidth[2] <- min(zoomwidth[2], ncol(DATA))
    if (sum(zoomwidth) == 0) {
        screen(1)
        for (i in 1:nrow(DATA)) {
            plot(1:length(DATA[i, ]), DATA[i, ], type = "l", 
                col = i, xlim = c(0, ncol(DATA) * 1.01), ylim = c(0, 
                  max(DATA[, min(start):max(stop)]) * 1.01), 
                main = c(maintext, paste("Target file: ", targetfile)), 
                xlab = "", ylab = datasource)
            par(new = TRUE)
        }
        if (!missing(ADATA)) {
            screen(2)
            for (i in 1:nrow(ADATA)) {
                plot(1:length(ADATA[i, ]), ADATA[i, ], type = "l", 
                  col = i, xlim = c(0, ncol(ADATA) * 1.01), ylim = c(0, 
                    max(DATA[, min(start):max(stop)]) * 1.01), 
                  main = "Aligned data", xlab = "Time (scans)", 
                  ylab = datasource)
                par(new = TRUE)
            }
        }
    }
    else {
        xlimit <- c(max(0, zoomwidth[1]), min(zoomwidth[2], ncol(DATA)))
        ylimit <- c(0, max(ADATA[, zoomwidth[1]:zoomwidth[2]]) * 
            1.01)
        screen(1)
        for (i in 1:nrow(DATA)) {
            plot(1:length(DATA[i, ]), DATA[i, ], type = "l", 
                col = i, xlim = xlimit, ylim = ylimit, main = c(maintext, 
                  paste("Target file: ", targetfile)), xlab = "", 
                ylab = datasource)
            par(new = TRUE)
        }
        if (!missing(ADATA)) {
            screen(2)
            for (i in 1:nrow(ADATA)) {
                plot(1:length(ADATA[i, ]), ADATA[i, ], type = "l", 
                  col = i, xlim = xlimit, ylim = ylimit, main = "Aligned data", 
                  xlab = "Time (scans)", ylab = datasource)
                par(new = TRUE)
            }
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
