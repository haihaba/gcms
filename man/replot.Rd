\name{replot}
\alias{replot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
replot(plotlimits = FALSE, plotedges = FALSE, plotfigure = FALSE, ymax, edges, DATA, num, zoomwidth, maintext, min_win, max_win)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plotlimits}{
%%     ~~Describe \code{plotlimits} here~~
}
  \item{plotedges}{
%%     ~~Describe \code{plotedges} here~~
}
  \item{plotfigure}{
%%     ~~Describe \code{plotfigure} here~~
}
  \item{ymax}{
%%     ~~Describe \code{ymax} here~~
}
  \item{edges}{
%%     ~~Describe \code{edges} here~~
}
  \item{DATA}{
%%     ~~Describe \code{DATA} here~~
}
  \item{num}{
%%     ~~Describe \code{num} here~~
}
  \item{zoomwidth}{
%%     ~~Describe \code{zoomwidth} here~~
}
  \item{maintext}{
%%     ~~Describe \code{maintext} here~~
}
  \item{min_win}{
%%     ~~Describe \code{min_win} here~~
}
  \item{max_win}{
%%     ~~Describe \code{max_win} here~~
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
function (plotlimits = FALSE, plotedges = FALSE, plotfigure = FALSE, 
    ymax, edges, DATA, num, zoomwidth, maintext, min_win, max_win) 
{
    if (num == 0 & !missing(DATA)) 
        ymax <- max(DATA) * 1.05
    if (plotfigure) {
        zoomwidth <- ifelse(num == 0, length(edges) - 1, zoomwidth)
        num <- max(num, 1)
        L <- edges[num + zoomwidth] - edges[num]
        xmin <- edges[num] - L/100
        xmax <- edges[num + zoomwidth] + L/100
        ymax <- max(DATA[, edges[num]:edges[num + zoomwidth]]) * 
            1.05
        par(new = FALSE)
        for (i in 1:nrow(DATA)) {
            plot(edges[num]:edges[num + zoomwidth], DATA[i, edges[num]:edges[num + 
                zoomwidth]], xlim = c(xmin, xmax), ylim = c(0, 
                ymax), type = "l", col = i, main = c(maintext, 
                paste("Window: ", paste("[", num, "-", num + 
                  zoomwidth - 1, "] (", length(edges) - 1, ")", 
                  sep = ""))), xlab = "", ylab = "")
            par(new = TRUE)
        }
    }
    if (plotlimits) {
        m <- edges_move_limit(edges, min_win, max_win)
        segments(x0 = edges, y0 = rep(ymax, length(edges)), x1 = m$right, 
            y1 = rep(ymax, length(edges)), lwd = 2.5, col = "blue")
        segments(x1 = m$left, y1 = rep(0, length(edges)), x0 = edges, 
            y0 = 0, lwd = 2.5, col = "green")
    }
    if (plotedges) 
        abline(v = edges, col = "red", lwd = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
