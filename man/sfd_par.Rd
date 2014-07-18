\name{sfd_par}
\alias{sfd_par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sfd_par(X, num, exclude = nrow(as.matrix(X)) + 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{num}{
%%     ~~Describe \code{num} here~~
}
  \item{exclude}{
%%     ~~Describe \code{exclude} here~~
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
function (X, num, exclude = nrow(as.matrix(X)) + 1) 
{
    N <- nrow(as.matrix(X))
    K <- ncol(as.matrix(X))
    if (num == (N - length(exclude))) {
        sel <- (1:N)[-exclude]
        return(list(object = sel))
    }
    D <- obj_eu(X)
    diag(D) <- 2 * max(D)
    MD <- j <- 0
    for (run in 1:500) {
        object <- sample((1:N)[-exclude])[1:num]
        OK <- 3
        md <- 0
        number <- 1
        while (OK) {
            XD <- D[object, object]
            OK <- OK - 1
            nk <- which.min(XD)
            k <- col(XD)[nk]
            if (OK) 
                n <- row(XD)[nk] * (OK != 2) + k * (OK == 2)
            else n <- number
            number <- number + (OK == 0)
            if (!OK) 
                OK <- (number <= num)
            obs <- object[-n]
            ex_object <- setdiff((1:N)[-exclude], obs)
            I <- which.max(apply(D[obs, ex_object], 2, min))
            for (jj in I) {
                obs[num] <- ex_object[jj]
                if ((md <= min(D[obs, obs])) & sum(abs(sort(obs) - 
                  sort(object))) > 0) {
                  XD <- D[obs, obs]
                  diag(XD) <- 0
                  if (md == min(D[obs, obs])) {
                    if (j < sum(XD)/2) {
                      md <- min(D[obs, obs])
                      object <- sort(obs)
                      j <- sum(XD)/2
                      number <- number * (OK != 1) + (OK == 1)
                      if (OK > 1) 
                        OK <- 3
                    }
                  }
                  else {
                    md <- min(D[obs, obs])
                    object <- sort(obs)
                    j <- sum(XD)/2
                    number <- number * (OK != 1) + (OK == 1)
                    if (OK > 1) 
                      OK <- 3
                  }
                }
            }
        }
        XD <- D[object, object]
        md <- min(D[object, object])
        diag(XD) <- 0
        j <- sum(XD)/2
        object <- sort(object)
        if (MD < md) {
            J <- j
            sel <- object
            MD <- md
        }
    }
    list(object = sel, j = J, md = MD)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
