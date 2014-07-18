\name{spec2NIST}
\alias{spec2NIST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
spec2NIST(projectpath, type, all = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{projectpath}{
%%     ~~Describe \code{projectpath} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{all}{
%%     ~~Describe \code{all} here~~
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
function (projectpath, type, all = FALSE) 
{
    load(file.path(projectpath, "Aligned", "SCAN_RANGE.Rdata"))
    load(file.path(projectpath, "HMCR", type, "MVA_DATA.Rdata"))
    load(file.path(projectpath, "HMCR", "MCR.Rdata"))
    excludedfiles <- eval(parse(text = paste("MCR$", tolower(type), 
        "$excluded", sep = "")))
    load(file.path(projectpath, "Aligned", "files.Rdata"))
    filenames <- sub("(.+)[.][^.]+$", "\\1", basename(files))
    if (!all) {
        cat("Select compounds to export.\n")
        choice <- which(VARID1 \%in\% select.list(VARID1, multiple = TRUE))
    }
    else choice <- 1:length(VARID1)
    if (!length(choice)) {
        cat("Error! No compounds selected!\n")
    }
    else {
        S <- t(SPECTRUM[, choice])
        RT <- VARID2[choice]
        NAME <- as.matrix(VARID1[choice])
        ID <- cbind(round(RT, 2), NAME)
        if (file.exists(file.path(projectpath, "ALKANE_SERIE", 
            "RT_INFO.Rdata"))) {
            load(file.path(projectpath, "ALKANE_SERIE", "RT_INFO.Rdata"))
            RI <- approxExtrap(RT_INFO[, 2], RT_INFO[, 1], RT)$y
        }
        else RI <- rep(-99, length(VARID1[choice]))
        if (length(excludedfiles)) 
            row.names(DATA) <- filenames[-excludedfiles]
        else row.names(DATA) <- filenames
        out <- data.frame(Window = VARID1[choice], RT_s = round(VARID2[choice], 
            2), RI, Annotated = character(length(VARID1[choice])), 
            Classified = character(length(VARID1[choice])), Comments = character(length(VARID1[choice])), 
            t(DATA)[choice, ])
        write.table(out, file = file.path(projectpath, "HMCR", 
            type, "Spectra.txt"), quote = FALSE, sep = "\t")
        if (file.exists(file.path(projectpath, "HMCR", type, 
            "Spectrum.txt"))) 
            file.remove(file.path(projectpath, "HMCR", type, 
                "Spectrum.txt"))
        fid <- file(file.path(projectpath, "HMCR", type, "Spectrum.txt"), 
            "a+")
        for (j in 1:nrow(ID)) {
            s <- S[j, ]
            s <- s/max(s) * 999
            s[round(s) < 1.5] <- 0
            N <- which(s == 0)
            if (length(N)) {
                s <- s[-N]
                MZ <- SCAN_RANGE[-N]
            }
            else MZ <- SCAN_RANGE
            DATA <- cbind(MZ, round(s))
            cat("Name:", paste(ID[j, ], collapse = " "), "\n", 
                file = fid, sep = "")
            cat("Synon:Annotated:\n", file = fid, sep = "")
            cat("Synon:Classified:\n", file = fid, sep = "")
            cat("Synon:RI:", round(RI[j], 1), "\n", file = fid, 
                sep = "")
            cat("Synon:RT(s):", round(RT[j], 1), "\n", file = fid, 
                sep = "")
            cat("Synon:MCR-Reference:", VARID1[choice][j], "\n", 
                file = fid, sep = "")
            cat("Synon:Date:", as.character(Sys.Date()), "\n", 
                file = fid, sep = "")
            cat("Comments:", "\n", file = fid, sep = "")
            cat("DB#:", j, "\n", file = fid, sep = "")
            cat("Num peaks:", length(MZ), "\n", file = fid, sep = "")
            cat(file = fid, paste(apply(DATA, 1, function(DATA) paste(DATA, 
                collapse = " ")), ";", rep(c(rep("", 3), "\n"), 
                len = length(MZ)), sep = "", collapse = " "), 
                "\n")
        }
        close(fid)
        cat("\nSpectrum exported!\n\n")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
