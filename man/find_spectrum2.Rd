\name{find_spectrum2}
\alias{find_spectrum2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
find_spectrum2(predpath, projectpath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predpath}{
%%     ~~Describe \code{predpath} here~~
}
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
function (predpath, projectpath) 
{
    require(MASS)
    X11.options(type = "cairo")
    datamenu <- character()
    if (file.exists(file.path(projectpath, "HMCR", "REG", "MVA_DATA.Rdata"))) 
        datamenu <- "REG H-MCR DATA"
    if (file.exists(file.path(projectpath, "HMCR", "CV", "MVA_DATA.Rdata"))) 
        datamenu <- c(datamenu, "CV H-MCR DATA")
    if (length(datamenu)) {
        dataexport <- menu(c(datamenu, "Cancel"), title = "Select data")
        if (dataexport == 0 | dataexport == length(c(datamenu, 
            "Cancel"))) 
            datamenu <- character()
    }
    if (!length(datamenu)) {
        cat("Error! (Aborted by user or no data found)\n\n")
        return(character())
    }
    else {
        type = ifelse(length(datamenu) == 1, strsplit(datamenu, 
            " ")[[1]][1], c("REG", "CV")[dataexport])
        load(file.path(predpath, "Aligned", "files.Rdata"))
        load(file.path(predpath, "Edges", "edges.Rdata"))
        load(file.path(predpath, "Aligned", "shift.Rdata"))
        load(file.path(predpath, "Aligned", "SCAN_RANGE.Rdata"))
        load(file.path(predpath, "maxMZ.Rdata"))
        load(files[which.min(shift)])
        load(file.path(projectpath, "HMCR", "MCR.Rdata"))
        windowlist <- as.numeric(substr(basename(list.files(file.path(projectpath, 
            "HMCR", type, "win"))), 4, 6))
        EDGES_TIME <- SCAN_INFO[, 2]
        load(file.path(projectpath, "SETTINGS.Rdata"))
        NL <- SETTINGS$NL
        RT_LIMIT <- SETTINGS$MPS
        DO_BL <- SETTINGS$BC2
        color <- cbind("red", "green", "blue", "black", "purple", 
            "grey", "yellow4", "red", "green", "blue", "black", 
            "purple", "grey", "yellow4", "red", "green", "blue", 
            "black", "purple", "grey", "yellow4", "red", "green", 
            "blue", "black", "purple", "grey", "yellow4", "red", 
            "green", "blue", "black", "purple", "grey", "yellow4")
        rm(Xbc, SETTINGS)
        dir.create(file.path(predpath, "Edges", "dat"), recursive = TRUE, 
            showWarnings = FALSE)
        temp <- list.files(file.path(projectpath, "Edges", "dat"), 
            full.names = TRUE)
        dir.create(file.path(predpath, "Edges", "dat", "Model samples_bg_corr"), 
            showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(predpath, "HMCR", type, "win_png"), 
            showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(predpath, "HMCR", type, "win"), 
            showWarnings = FALSE, recursive = TRUE)
        gc()
        Scores <- numeric()
        if (length(windowlist)) {
            for (i in 1:length(files)) {
                load(files[i])
                cat("Resolving ", basename(files[i]), paste(" (", 
                  i, "/", length(files), ")\n", sep = ""))
                for (win in windowlist) {
                  if (edges[win] + shift[i] > 0 & edges[win + 
                    1] + shift[i] < nrow(Xbc)) {
                    x <- Xbc[edges[win]:edges[win + 1] + shift[i], 
                      ]
                    if (ncol(x) < maxMZ) 
                      x <- cbind(x, matrix(0, nrow = nrow(x), 
                        ncol = maxMZ - ncol(x)))
                    if (DO_BL) {
                      BL <- matrix(apply(x, 2, function(x) approx(c(1, 
                        length(x)), x[c(1, length(x))], 1:length(x))$y), 
                        nrow = length(edges[win]:edges[win + 
                          1]), ncol = maxMZ)
                      x <- x - BL
                      x[x < 0] <- 0
                      rm(BL)
                    }
                    load(file.path(projectpath, "HMCR", "REG", 
                      "win", paste("win", ifelse(win <= 99, ifelse(win <= 
                        9, paste("00", win, sep = ""), paste("0", 
                        win, sep = "")), win), ".Rdata", sep = "")))
                    if (!length(S)) 
                      c2 <- numeric()
                    else c2 <- do_AR_all_prediction(x, S, CP, 
                      RT_LIMIT)
                    if (exists("PCApara")) 
                      scores <- (colSums(x) <- PCApara$mean) \%*\% 
                        PCApara$loadings[, 1:2]
                    else scores <- numeric()
                    ok <- 1
                  }
                  else {
                    load(file.path(projectpath, "HMCR", "REG", 
                      "win", paste("win", ifelse(win <= 99, ifelse(win <= 
                        9, paste("00", win, sep = ""), paste("0", 
                        win, sep = "")), win), ".Rdata", sep = "")))
                    c2 <- matrix(0, length(edges[win]:edges[win + 
                      1]), ncol(as.matrix(C)))
                    ok <- 0
                    if (exists("PCApara")) 
                      scores <- c(0, 0)
                    else scores <- numeric()
                  }
                  if (i == 1) {
                    if (ok) {
                      C <- matrix(colSums(as.matrix(c2, ncol = ncol(as.matrix(C)))), 
                        nrow = 1)
                      CC <- c2
                    }
                    else {
                      C <- matrix(colSums(as.matrix(c2, ncol = ncol(as.matrix(C)))), 
                        nrow = 1) * NA
                      CC <- c2
                    }
                    Scores <- rbind(Scores, scores)
                    save(C, CC, S, TIME, Scores, file = file.path(predpath, 
                      "HMCR", type, "win", paste("win", ifelse(win <= 
                        99, ifelse(win <= 9, paste("00", win, 
                        sep = ""), paste("0", win, sep = "")), 
                        win), ".Rdata", sep = "")))
                  }
                  else {
                    load(file.path(predpath, "HMCR", type, "win", 
                      paste("win", ifelse(win <= 99, ifelse(win <= 
                        9, paste("00", win, sep = ""), paste("0", 
                        win, sep = "")), win), ".Rdata", sep = "")))
                    Scores <- rbind(Scores, scores)
                    if (ok) {
                      C <- rbind(C, colSums(as.matrix(c2, ncol = ncol(as.matrix(C)))))
                      CC <- rbind(CC, c2)
                    }
                    else {
                      C <- rbind(C, colSums(as.matrix(c2, ncol = ncol(as.matrix(C)))) * 
                        NA)
                      CC <- rbind(CC, c2)
                    }
                    save(C, CC, S, TIME, Scores, file = file.path(predpath, 
                      "HMCR", type, "win", paste("win", ifelse(win <= 
                        99, ifelse(win <= 9, paste("00", win, 
                        sep = ""), paste("0", win, sep = "")), 
                        win), ".Rdata", sep = "")))
                  }
                  gc()
                }
                gc()
            }
        }
        else cat("No windows processed.\n")
        return(list(type = type, windowlist = windowlist))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
