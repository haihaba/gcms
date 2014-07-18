\name{scaling}
\alias{scaling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scaling(Xin, K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xin}{
%%     ~~Describe \code{Xin} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
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
function (Xin, K) 
{
    if (missing(K)) 
        K <- menu(c("Centering", "UV-Scaling", "Pareto", "None"), 
            title = "Select scaling function for X-matrix")
    if (K == 1) {
        Xout <- cent(Xin)
        scale_text <- "DATA is mean centered"
    }
    else if (K == 2) {
        Xout <- uv_sc(Xin)
        scale_text <- "DATA is centered and UV-scaled"
    }
    else if (K == 3) {
        Xout <- par_sc(Xin)
        scale_text <- "DATA is centered and Pareto-scaled"
    }
    else if (K == 4 || K == 0) {
        Xout <- Xin
        scale_text <- "DATA neither centered nor scaled"
    }
    scale <- list(Xout = Xout, scale_text = scale_text)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
