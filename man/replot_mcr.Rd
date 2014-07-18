\name{replot_mcr}
\alias{replot_mcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
replot_mcr(x, y, maintext, incl, excl, numfiles, set1, set2, seltype = (!missing(set1) & !missing(set2)), plotlabels = TRUE, pred)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{maintext}{
%%     ~~Describe \code{maintext} here~~
}
  \item{incl}{
%%     ~~Describe \code{incl} here~~
}
  \item{excl}{
%%     ~~Describe \code{excl} here~~
}
  \item{numfiles}{
%%     ~~Describe \code{numfiles} here~~
}
  \item{set1}{
%%     ~~Describe \code{set1} here~~
}
  \item{set2}{
%%     ~~Describe \code{set2} here~~
}
  \item{seltype}{
%%     ~~Describe \code{seltype} here~~
}
  \item{plotlabels}{
%%     ~~Describe \code{plotlabels} here~~
}
  \item{pred}{
%%     ~~Describe \code{pred} here~~
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
function (x, y, maintext, incl, excl, numfiles, set1, set2, seltype = (!missing(set1) & 
    !missing(set2)), plotlabels = TRUE, pred) 
{
    if (seltype == 1) {
        color <- 2 * (1:numfiles \%in\% set1) + 4 * (1:numfiles \%in\% 
            set2) + 3 * (1:numfiles \%in\% pred) + (1:numfiles \%in\% 
            setdiff(1:numfiles, c(set1, set2, pred)))
    }
    else {
        color <- 2 * (1:numfiles \%in\% incl) + (1:numfiles \%in\% 
            excl) + 3 * (1:numfiles \%in\% pred)
    }
    plot(x, y, xlim = c(range(x)), ylim = range(y), main = maintext, 
        xlab = "t[1]", ylab = "t[2]", col = color, cex = 0.9, 
        pch = 16)
    if (plotlabels) 
        text(x, y, pos = 4, cex = 0.7)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
