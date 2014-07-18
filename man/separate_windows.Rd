\name{separate_windows}
\alias{separate_windows}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
separate_windows(projectpath, windowlist, files, edges, model, prediction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{projectpath}{
%%     ~~Describe \code{projectpath} here~~
}
  \item{windowlist}{
%%     ~~Describe \code{windowlist} here~~
}
  \item{files}{
%%     ~~Describe \code{files} here~~
}
  \item{edges}{
%%     ~~Describe \code{edges} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{prediction}{
%%     ~~Describe \code{prediction} here~~
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
function (projectpath, windowlist, files, edges, model, prediction) 
{
    if (length(windowlist)) {
        load(file.path(projectpath, "maxMZ.Rdata"))
        load(file.path(projectpath, "Aligned", "shift.Rdata"))
        dir.create(file.path(projectpath, "Edges", "dat", "Model samples"), 
            showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(projectpath, "Edges", "dat", "Prediction samples"), 
            showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(projectpath, "Edges", "dat", "Metadata"), 
            showWarnings = FALSE, recursive = TRUE)
        cat("Checking integrity of existing window data files..\n")
        modelfilespath <- list.files(file.path(projectpath, "Edges", 
            "dat", "Model samples"), full.names = TRUE)
        predictionfilespath <- list.files(file.path(projectpath, 
            "Edges", "dat", "Prediction samples"), full.names = TRUE)
        if (length(files) == length(c(model, prediction))) {
            excludedfiles <- numeric()
        }
        else {
            excludedfiles <- (1:length(files))[-c(model, prediction)]
        }
        newlist <- logical(max(windowlist))
        metadata <- list.files(file.path(projectpath, "Edges", 
            "dat", "Metadata"), full.names = TRUE)
        metadata <- metadata[as.numeric(substr(basename(metadata), 
            4, 6)) \%in\% windowlist]
        newlist[windowlist[!(windowlist \%in\% as.numeric(substr(basename(metadata), 
            4, 6)))]] <- TRUE
        newlist[windowlist[!(windowlist \%in\% as.numeric(substr(basename(modelfilespath), 
            4, 6)))]] <- TRUE
        newlist[windowlist[!(windowlist \%in\% as.numeric(substr(basename(predictionfilespath), 
            4, 6)))]] <- TRUE
        if (length(metadata)) {
            n <- 0
            for (k in as.numeric(substr(basename(metadata), 4, 
                6))) {
                n <- n + 1
                load(metadata[n])
                if (!(all.equal(modelfiles, model) == TRUE & 
                  all.equal(predictionfiles, prediction) == TRUE)) 
                  newlist[k] <- TRUE
            }
        }
        windowlist <- which(newlist)
        if (any(as.numeric(substr(basename(modelfilespath), 4, 
            6)) \%in\% windowlist)) 
            file.remove(modelfilespath[which(as.numeric(substr(basename(modelfilespath), 
                4, 6)) \%in\% windowlist)])
        if (any(as.numeric(substr(basename(predictionfilespath), 
            4, 6)) \%in\% windowlist)) 
            file.remove(predictionfilespath[which(as.numeric(substr(basename(predictionfilespath), 
                4, 6)) \%in\% windowlist)])
        cat("Windows to separate: [ ", windowlist, " ]\n")
        cat("In total", length(windowlist), "windows\n")
        filenum <- 0
        modelfiles <- model
        predictionfiles <- prediction
        for (i in sort(c(model, prediction))) {
            filenum <- filenum + 1
            cat("Loading ", basename(files[i]), paste("(", filenum, 
                "/", length(c(model, prediction)), ") ", sep = ""))
            load(files[i])
            cat("... Exporting windows... \n")
            for (j in windowlist) {
                if (nrow(Xbc) < (edges[j] + shift[i])) {
                  a <- 0
                }
                else if (nrow(Xbc) < (edges[j + 1] + shift[i])) {
                  a <- (edges[j] + shift[i]):nrow(Xbc)
                }
                else a <- edges[j]:edges[j + 1] + shift[i]
                Xbc <- cbind(Xbc, matrix(0, nrow = nrow(Xbc), 
                  ncol = ifelse(maxMZ > ncol(Xbc), maxMZ - ncol(Xbc), 
                    0)))
                x <- rbind(Xbc[a, ], matrix(0, nrow = sum((0:(edges[j + 
                  1] + shift[i] - nrow(Xbc))) > 0), ncol = ncol(Xbc)))
                write.table(t(x), file = file.path(projectpath, 
                  "Edges", "dat", ifelse(i \%in\% model, "Model samples", 
                    "Prediction samples"), paste("win", ifelse(j <= 
                    99, ifelse(j <= 9, paste("00", j, sep = ""), 
                    paste("0", j, sep = "")), j), ".dat", sep = "")), 
                  append = TRUE, row.names = FALSE, col.names = FALSE)
            }
        }
        for (j in windowlist) save(modelfiles, predictionfiles, 
            excludedfiles, file = file.path(projectpath, "Edges", 
                "dat", "Metadata", paste("win", ifelse(j <= 99, 
                  ifelse(j <= 9, paste("00", j, sep = ""), paste("0", 
                    j, sep = "")), j), "_metadata.Rdata", sep = "")))
    }
    else cat("No windows selected!\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
