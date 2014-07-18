\name{aligndatapred}
\alias{aligndatapred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
aligndatapred(projectpath, predpath, newsamples, target, minshift, datasource, ion)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{projectpath}{
%%     ~~Describe \code{projectpath} here~~
}
  \item{predpath}{
%%     ~~Describe \code{predpath} here~~
}
  \item{newsamples}{
%%     ~~Describe \code{newsamples} here~~
}
  \item{target}{
%%     ~~Describe \code{target} here~~
}
  \item{minshift}{
%%     ~~Describe \code{minshift} here~~
}
  \item{datasource}{
%%     ~~Describe \code{datasource} here~~
}
  \item{ion}{
%%     ~~Describe \code{ion} here~~
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
function (projectpath, predpath, newsamples, target, minshift, 
    datasource, ion) 
{
    load(file.path(projectpath, "maxMZ.Rdata"))
    alignfiles <- newsamples[order(newsamples)]
    A_DATA <- list()
    n <- length(alignfiles)
    dir.create(file.path(predpath, "Aligned"), showWarnings = FALSE)
    COLUMNID1 <- cbind(sub("[.][^.]*$", "", basename(alignfiles)))
    A_DATA$files <- files <- alignfiles
    A_DATA$fsize <- matrix(0, n, 2)
    for (i in 1:n) {
        cat("Loading ", paste(basename(alignfiles[i]), " (", 
            i, "/", n, ")\n", sep = ""))
        load(alignfiles[i])
        A_DATA$fsize[i, ] <- dim(Xbc)
        num <- which(SCAN_RANGE == ion)
        if (i == 1) {
            NUM_MZ <- ncol(Xbc)
            NUM_scans <- nrow(Xbc)
            A_DATA$tic <- matrix(0, n, NUM_scans)
            A_DATA$bp <- matrix(0, n, NUM_scans)
            A_DATA$ic <- matrix(0, n, NUM_scans)
            SUM_MZ <- matrix(0, n, NUM_MZ)
        }
        if (nrow(Xbc) > NUM_scans) {
            A_DATA$tic <- cbind(A_DATA$tic, matrix(0, n, nrow(Xbc) - 
                NUM_scans))
            A_DATA$bp <- cbind(A_DATA$bp, matrix(0, n, nrow(Xbc) - 
                NUM_scans))
            A_DATA$ic <- cbind(A_DATA$ic, matrix(0, n, nrow(Xbc) - 
                NUM_scans))
            NUM_scans <- nrow(Xbc)
        }
        if (ncol(Xbc) > NUM_MZ) 
            SUM_MZ <- cbind(SUM_MZ, matrix(0, n, ncol(Xbc) - 
                NUM_MZ))
        A_DATA$tic[i, 1:nrow(Xbc)] <- rowSums(Xbc)
        A_DATA$bp[i, 1:nrow(Xbc)] <- apply(Xbc, 1, max)
        A_DATA$ic[i, 1:nrow(Xbc)] <- t(Xbc[, num])
        SUM_MZ[i, 1:ncol(Xbc)] <- colSums(Xbc)
        if (ncol(Xbc) < maxMZ) {
            Xbc <- cbind(Xbc, matrix(0, nrow = nrow(Xbc), ncol = maxMZ - 
                ncol(Xbc)))
            SCAN_RANGE <- c(SCAN_RANGE, (max(SCAN_RANGE) + 1):maxMZ)
            save(Xbc, SCAN_INFO, SCAN_RANGE, file = alignfiles[i])
        }
        else if (ncol(Xbc) >= maxMZ) {
            Xbc <- Xbc[, 1:maxMZ]
            SCAN_RANGE <- SCAN_RANGE[1:maxMZ]
            save(Xbc, SCAN_INFO, SCAN_RANGE, file = alignfiles[i])
        }
    }
    A_DATA$ion <- ion
    A_DATA$ok <- 1
    if (file.exists(file.path(predpath, "SETTINGS.Rdata"))) {
        load(file.path(predpath, "SETTINGS.Rdata"))
        max_shift <- SETTINGS$MS
    }
    else {
        cat("No settings for maximum shift found! Using default value (500)")
        max_shift <- 500
    }
    if (datasource == "TIC") {
        DATA <- A_DATA$tic
    }
    else if (datasource == "Basepeak") {
        DATA <- A_DATA$bp
    }
    else if (datasource == "IC") 
        DATA <- A_DATA$ic
    if (length(target) < ncol(DATA)) 
        target <- c(target, rep(0, ncol(DATA) - length(target)))
    shift <- numeric()
    for (i in 1:length(alignfiles)) {
        CO_VAR <- numeric()
        A <- target[(max_shift + 1):(ncol(DATA) - max_shift)]
        for (j in -max_shift:max_shift) {
            B <- DATA[i, (max_shift + 1 + j):(ncol(DATA) - max_shift + 
                j)]
            CO_VAR <- c(CO_VAR, sum(A * B))
        }
        shift <- c(shift, (which.max(CO_VAR) - max_shift - 1))
    }
    shift <- shift - minshift
    cat("Min shift: ", min(shift), "\n")
    cat("Max shift: ", max(shift), "\n")
    if (min(shift) < 0) 
        cat("Adjusting negative shifts...\n")
    shift <- shift - min(shift)
    start <- 1 + shift[1:nrow(A_DATA$tic)]
    stop <- ncol(A_DATA$tic) - max(shift) + shift[1:nrow(A_DATA$tic)]
    BASEPEAK <- matrix(0, nrow(A_DATA$bp), ncol(A_DATA$bp))
    TIC <- matrix(0, nrow(A_DATA$tic), ncol(A_DATA$tic))
    for (i in 1:nrow(A_DATA$bp)) {
        TIC[i, 1:(stop[i] - start[i] + 1)] <- A_DATA$tic[i, start[i]:stop[i]]
        BASEPEAK[i, 1:(stop[i] - start[i] + 1)] <- A_DATA$bp[i, 
            start[i]:stop[i]]
    }
    TIC_RAW <- A_DATA$tic
    cat("Saving variables..\n")
    save(TIC_RAW, file = file.path(predpath, "Aligned", "TIC_RAW.Rdata"))
    save(shift, file = file.path(predpath, "Aligned", "shift.Rdata"))
    save(TIC, file = file.path(predpath, "Aligned", "TIC.Rdata"))
    save(BASEPEAK, file = file.path(predpath, "Aligned", "BASEPEAK.Rdata"))
    save(SUM_MZ, file = file.path(predpath, "Aligned", "SUM_MZ.Rdata"))
    save(COLUMNID1, file = file.path(predpath, "Aligned", "COLUMNID1.Rdata"))
    save(files, file = file.path(predpath, "Aligned", "files.Rdata"))
    save(SCAN_RANGE, file = file.path(predpath, "Aligned", "SCAN_RANGE.Rdata"))
    save(NUM_scans, file = file.path(predpath, "Aligned", "NUM_scans.Rdata"))
    textstr <- paste(files, " Shift: ", shift, " Scans: ", A_DATA$fsize[, 
        1], " MZ: ", A_DATA$fsize[, 2], sep = "\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
