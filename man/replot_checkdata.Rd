\name{replot_checkdata}
\alias{replot_checkdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
replot_checkdata(x, y, scale_text, maintext, OBS, ssX, s, coords = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{scale_text}{
%%     ~~Describe \code{scale_text} here~~
}
  \item{maintext}{
%%     ~~Describe \code{maintext} here~~
}
  \item{OBS}{
%%     ~~Describe \code{OBS} here~~
}
  \item{ssX}{
%%     ~~Describe \code{ssX} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
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
function (x, y, scale_text, maintext, OBS, ssX, s, coords = NULL) 
{
    xtext <- paste("Component 1 (", round(t(x) \%*\% x/ssX * 100, 
        1), "\% of variation in X)", sep = "")
    ytext <- paste("Component 2 (", round(t(y) \%*\% y/ssX * 100, 
        1), "\% of variation in X)", sep = "")
    color <- s \%in\% s[grep("[OUT]", s)] + 1
    if (missing(coords)) {
        plot(x, y, main = c(maintext, scale_text), xlab = xtext, 
            ylab = ytext, col = color, pch = 16)
        text(x, y, OBS, pos = 4, col = color)
    }
    else {
        plot(x, y, main = c(maintext, scale_text), xlab = xtext, 
            ylab = ytext, col = color, pch = 16, xlim = coords$x, 
            ylim = coords$y)
        text(x, y, OBS, pos = 4, col = color)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
