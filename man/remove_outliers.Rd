\name{remove_outliers}
\alias{remove_outliers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
remove_outliers(projectpath, outliers)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{projectpath}{
%%     ~~Describe \code{projectpath} here~~
}
  \item{outliers}{
%%     ~~Describe \code{outliers} here~~
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
function (projectpath, outliers) 
{
    OK <- menu(c("Yes", "No"), title = "Remove outliers from further analysis?")
    if (OK == 1 & length(outliers) & is.numeric(outliers)) {
        cat("Removing outliers...\n\n")
        load(file.path(projectpath, "Aligned", "TIC.Rdata"))
        load(file.path(projectpath, "Aligned", "BASEPEAK.Rdata"))
        load(file.path(projectpath, "Aligned", "SUM_MZ.Rdata"))
        load(file.path(projectpath, "Aligned", "shift.Rdata"))
        load(file.path(projectpath, "Aligned", "files.Rdata"))
        load(file.path(projectpath, "Aligned", "COLUMNID1.Rdata"))
        files <- files[-outliers]
        shift <- shift[-outliers]
        BASEPEAK <- BASEPEAK[-outliers, ]
        TIC <- TIC[-outliers, ]
        SUM_MZ <- SUM_MZ[-outliers, ]
        COLUMNID1 <- COLUMNID1[-outliers]
        save(TIC, file = file.path(projectpath, "Aligned", "TIC.Rdata"))
        save(BASEPEAK, file = file.path(projectpath, "Aligned", 
            "BASEPEAK.Rdata"))
        save(SUM_MZ, file = file.path(projectpath, "Aligned", 
            "SUM_MZ.Rdata"))
        save(shift, file = file.path(projectpath, "Aligned", 
            "shift.Rdata"))
        save(files, file = file.path(projectpath, "Aligned", 
            "files.Rdata"))
        save(COLUMNID1, file = file.path(projectpath, "Aligned", 
            "COLUMNID1.Rdata"))
    }
    else cat("Aborted by user or no outliers selected!\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
