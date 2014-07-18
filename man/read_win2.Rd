\name{read_win2}
\alias{read_win2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read_win2(projectpath, type, win)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{projectpath}{
%%     ~~Describe \code{projectpath} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{win}{
%%     ~~Describe \code{win} here~~
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
function (projectpath, type, win) 
{
    SPECTRUM <- DATA <- VARID1 <- VARID2 <- numeric()
    load(file.path(projectpath, "Edges", "edges.Rdata"))
    num_windows <- length(edges) - 1
    load(file.path(projectpath, "Aligned", "files.Rdata"))
    if (!missing(win)) {
        if (length(list.files(file.path(projectpath, "HMCR", 
            type, "win"))) > length(win)) {
            cat("\nRecently processed windows: ", win, "\n")
            cat("Previously processed windows: ", setdiff(as.numeric(substr(basename(list.files(file.path(projectpath, 
                "HMCR", type, "win"))), 4, 6)), win), "\n")
            if (1 == menu(c("Yes", "No"), title = "There are some processed windows that are not selected, do you wish to include some of them as well?")) 
                while (!length(win <- sort(as.numeric(substr(basename(tk_choose.files(caption = "Select windows to process")), 
                  4, 6))))) cat("You must select at least one window!\n")
        }
    }
    else {
        cat("Select windows.")
        while (!length(win <- sort(as.numeric(substr(basename(tk_choose.files(caption = "Select windows to process")), 
            4, 6))))) cat("You must select at least one window!\n")
    }
    cat("\n")
    for (i in win) {
        if (file.exists(file.path(projectpath, "HMCR", type, 
            "win", paste("win", ifelse(i <= 99, ifelse(i <= 9, 
                paste("00", i, sep = ""), paste("0", i, sep = "")), 
                i), ".Rdata", sep = "")))) {
            load(file.path(projectpath, "HMCR", type, "win", 
                paste("win", ifelse(i <= 99, ifelse(i <= 9, paste("00", 
                  i, sep = ""), paste("0", i, sep = "")), i), 
                  ".Rdata", sep = "")))
            if (length(S)) {
                cat(nrow(as.matrix(C)), "\n")
                if (nrow(as.matrix(C)) == length(files)) {
                  cat("--------------------------------------------------------\n")
                  SPECTRUM <- cbind(SPECTRUM, S)
                  DATA <- cbind(DATA, C)
                  cat("Window num: ", i, "\nNumber of comps: ", 
                    ncol(S), "\nTotal number of comps: ", ncol(SPECTRUM), 
                    "\n")
                  VARID2 <- c(VARID2, TIME)
                  for (j in 1:ncol(S)) VARID1 <- rbind(VARID1, 
                    ifelse(i > 99, ifelse(j < 10, paste("Win", 
                      i, "_C0", j, sep = ""), paste("Win ", i, 
                      "_C", j, sep = "")), ifelse(i < 10, ifelse(j < 
                      10, paste("Win00", i, "_C0", j, sep = ""), 
                      paste("Win00", i, "_C", j, sep = "")), 
                      ifelse(j < 10, paste("Win0", i, "_C0", 
                        j, sep = ""), paste("Win0", i, "_C", 
                        j, sep = "")))))
                }
                else cat("Number of files and rows in data matrix does not match!\n")
            }
        }
    }
    load(file.path(projectpath, "Aligned", "COLUMNID1.Rdata"))
    OBSID1 <- COLUMNID1
    save(DATA, VARID1, VARID2, OBSID1, SPECTRUM, file = file.path(projectpath, 
        "HMCR", type, "MVA_DATA.Rdata"))
    save(SPECTRUM = SPECTRUM, file = file.path(projectpath, "HMCR", 
        type, "SPECTRUM.Rdata"))
    save(VARID1 = VARID1, file = file.path(projectpath, "HMCR", 
        type, "VARID1.Rdata"))
    save(VARID2 = VARID2, file = file.path(projectpath, "HMCR", 
        type, "VARID2.Rdata"))
    write.table(data.frame(ID = OBSID1, DATA), file = file.path(projectpath, 
        "HMCR", type, "MVA_DATA.txt"), sep = "\t", row.names = FALSE, 
        col.names = c("Primary ID", VARID1), quote = FALSE)
    out <- list(DATA = DATA, SPECTRUM = SPECTRUM, VARID1 = VARID1, 
        OBSID1 = OBSID1, VARID2 = VARID2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
