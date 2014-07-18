\name{replot2}
\alias{replot2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
replot2(edges, DATA, zoomwidth, maintext)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edges}{
%%     ~~Describe \code{edges} here~~
}
  \item{DATA}{
%%     ~~Describe \code{DATA} here~~
}
  \item{zoomwidth}{
%%     ~~Describe \code{zoomwidth} here~~
}
  \item{maintext}{
%%     ~~Describe \code{maintext} here~~
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
function (edges, DATA, zoomwidth, maintext) 
{
    if (sum(zoomwidth) == 0) {
        zoomwidth[1] <- 1
        zoomwidth[2] <- ncol(DATA)
    }
    else {
        zoomwidth[1] <- max(1, zoomwidth[1])
        zoomwidth[2] <- min(zoomwidth[2], ncol(DATA))
    }
    ymax <- max(DATA[, zoomwidth[1]:zoomwidth[2]]) * 1.05
    par(new = FALSE)
    for (i in 1:nrow(DATA)) {
        plot(zoomwidth[1]:zoomwidth[2], DATA[i, zoomwidth[1]:zoomwidth[2]], 
            xlim = c(zoomwidth[1], zoomwidth[2]), ylim = c(0, 
                ymax), type = "l", col = i, main = maintext, 
            xlab = "", ylab = "")
        par(new = TRUE)
    }
    abline(v = edges, col = "red", lwd = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
