\name{unimodal_3}
\alias{unimodal_3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
unimodal_3(C, obs, cp, RT_LIMIT)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{C}{
%%     ~~Describe \code{C} here~~
}
  \item{obs}{
%%     ~~Describe \code{obs} here~~
}
  \item{cp}{
%%     ~~Describe \code{cp} here~~
}
  \item{RT_LIMIT}{
%%     ~~Describe \code{RT_LIMIT} here~~
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
function (C, obs, cp, RT_LIMIT) 
{
    for (i in 1:obs) {
        xpeak <- peak_pick(t(C[, i]))$xpeak
        mp <- which(as.logical(xpeak))
        if (!length(mp)) {
            C[, i] <- C[, i] * 0
        }
        else {
            if (which.min(abs(mp - cp))) 
                mp <- min(mp[which.min(abs(mp - cp))])
            else mp <- 1
            if (abs(mp - cp) < RT_LIMIT) {
                D <- diff(C[1:mp, i])
                if (length(which(D < 0))) 
                  poss <- max(which(D < 0))
                else poss <- 1
                for (j in poss:1) C[j, i] <- min(C[j:mp, i])
                D <- diff(C[mp:nrow(C), i])
                poss <- if (length(which(D > 0))) 
                  poss <- min(which(D > 0))
                else poss <- FALSE
                if (poss) {
                  for (j in (mp - 1 + poss):nrow(C)) C[j, i] <- min(C[mp:j, 
                    i])
                }
                knew <- C[(diff(C[, i]) != 0), i]
                mpnew <- which.max(knew)
                k <- C[, i] * 0
                k[1:length(knew) + mp - mpnew] <- knew
                C[, i] <- k
            }
            else C[, i] <- C[, i] * 0
        }
    }
    C
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
