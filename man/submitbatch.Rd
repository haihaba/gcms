\name{submitbatch}
\alias{submitbatch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
submitbatch(projectpath, cores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{projectpath}{
%%     ~~Describe \code{projectpath} here~~
}
  \item{cores}{
%%     ~~Describe \code{cores} here~~
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
function (projectpath, cores) 
{
    if (!file.exists(file.path(projectpath, "batch.R"))) {
        cat("library(batch)\n", file = "batch.R")
        cat("projectpath<-getwd()\n", file = "batch.R", append = T)
        cat("library(GCMS)\n", file = "batch.R", append = T)
        cat("load(file.path(projectpath,\"HMCR\",\"MCR.Rdata\"))\n", 
            file = "batch.R", append = T)
        cat("wincores<-0\n", file = "batch.R", append = T)
        cat("parseCommandArgs()\n", file = "batch.R", append = T)
        cat("samples<-MCR$reg$samples\n", file = "batch.R", append = T)
        cat("incl<-MCR$reg$incl\n", file = "batch.R", append = T)
        cat("excl<-MCR$reg$excl\n", file = "batch.R", append = T)
        cat("pred<-MCR$reg$pred\n", file = "batch.R", append = T)
        cat("MCR$reg$excluded<-excl\n", file = "batch.R", append = T)
        cat("win<-find_spectrum(projectpath,incl,pred,wincores)\n", 
            file = "batch.R", append = T)
    }
    load(file.path(projectpath, "HMCR", "MCR.Rdata"))
    winlist <- as.numeric(substr(MCR$windowlist[grep("[P]", MCR$windowlist)], 
        8, 10))
    if (length(winlist) < cores) 
        cores <- length(winlist)
    batches <- list()
    corecounter <- 1
    for (i in winlist) {
        eval(parse(text = paste("batches$core", corecounter, 
            "<-c(batches$core", corecounter, ",i)", sep = "")))
        corecounter <- corecounter + 1
        if (corecounter > cores) 
            corecounter <- 1
    }
    seed <- 1000
    for (i in 1:cores) {
        eval(parse(text = paste("k<-batches$core", i, sep = "")))
        a <- as.character()
        for (j in 1:length(k)) {
            a <- paste(a, ",", k[j], sep = "")
        }
        system(paste("R --vanilla --args wincores \"c(", substr(a, 
            2, nchar(a)), ")\" < \"batch.R\"", sep = ""), wait = F, 
            ignore.stdout = T)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
