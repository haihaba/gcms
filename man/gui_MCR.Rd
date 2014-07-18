\name{gui_MCR}
\alias{gui_MCR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gui_MCR(projectpath)
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
    if (file.exists(file.path(projectpath, "Aligned", "files.Rdata")) | 
        file.exists(file.path(projectpath, "Aligned", "shift.Rdata"))) {
        load(file.path(projectpath, "Aligned", "files.Rdata"))
        load(file.path(projectpath, "Aligned", "shift.Rdata"))
    }
    else {
        cat("Error! You must align data first!\n")
        Sys.sleep(1)
        return(NULL)
    }
    load(files[which.min(shift)])
    EDGES_TIME <- SCAN_INFO[, 2]
    save(EDGES_TIME, file = file.path(projectpath, "EDGES_TIME.Rdata"))
    out <- 2
    temp <- 5
    while (out == 2 | temp == 5) {
        if (file.exists(file.path(projectpath, "HMCR", "MCR.Rdata"))) {
            load(file.path(projectpath, "HMCR", "MCR.Rdata"))
            cat("Current settings:\n\n")
            cat("Number of samples: ", length(files), "\n")
            cat("Windows to be processed: ", grep("[P]", MCR$windowlist), 
                "\n\n")
            cat("Regular method settings:\n")
            cat("Included files: ", MCR$reg$incl, "\n")
            cat("Prediction files: ", MCR$reg$pred, "\n")
            cat("Excluded files: ", MCR$reg$excl, "\n\n")
            cat("Similarity criterion:\n")
            cat("Chrom. profile: ", MCR$cp, "\n")
            cat("Spec. profile: ", MCR$sp, "\n\n")
            out <- 1
            if (out == 3) 
                temp <- 6
            if (out == 1) 
                temp <- 1
        }
        if (out == 2) {
            MCR <- MCR_settings(projectpath)
            if (is.null(MCR)) 
                return(NULL)
        }
    }
    if (temp == 1) {
        cat("======== Regular method ========\n\n")
        dir.create(file.path(projectpath, "HMCR", "REG"), recursive = TRUE, 
            showWarnings = FALSE)
        samples <- MCR$reg$samples
        incl <- MCR$reg$incl
        excl <- MCR$reg$excl
        pred <- MCR$reg$pred
        MCR$reg$excluded <- excl
        save(MCR = MCR, file = file.path(projectpath, "HMCR", 
            "MCR.Rdata"))
        win <- find_spectrum(projectpath, incl, pred, grep("[P]", 
            MCR$windowlist))
        window_data <- read_win(projectpath, "REG", win)
        warnings()
    }
    if (temp == 1) {
    }
    cat("Ending Method 2.\n\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
