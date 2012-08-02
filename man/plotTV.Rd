\name{plotTV}
\alias{plotTV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plotTV
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotTV(tvScore, tvScore2, narrays, increm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tvScore}{
%%     ~~Describe \code{tvScore} here~~
}
  \item{tvScore2}{
%%     ~~Describe \code{tvScore2} here~~
}
  \item{narrays}{
%%     ~~Describe \code{narrays} here~~
}
  \item{increm}{
%%     ~~Describe \code{increm} here~~
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
function (tvScore, tvScore2, narrays, increm) 
{
    colinit = brewer.pal(8, "Set1")
    pal <- colorRampPalette(colinit)
    cols = pal(narrays)
    forind1 = 1:(nrow(tvScore)/2)
    revind1 = ((nrow(tvScore)/2) + 1):nrow(tvScore)
    forind2 = ((nrow(tvScore2)/2) + 1):nrow(tvScore2)
    revind2 = 1:(nrow(tvScore2)/2)
    formaxinds1 = tvScoreFirstHalfMax(tvScore)[, 1]
    revmaxinds1 = tvScoreSecondHalfMax(tvScore)[, 1]
    par(mfrow = c(2, 1))
    xlimuse = c(min(-forind1 * increm), max(forind1 * increm))
    ylimuse = (c(0, max(c(tvScore2, tvScore))))
    matplot(forind1 * increm, tvScore[forind1, ], xlab = "Window Length", 
        ylab = "TV Score", cex = 0.8, col = cols, type = "p", 
        pch = ">", main = "TV Scores For Different Windows Extending In Both Directions From 0", 
        xlim = xlimuse, ylim = ylimuse)
    matpoints(-forind1 * increm, tvScore[revind1, ], xlab = "Window Length", 
        ylab = "TV Score", cex = 0.8, col = cols, type = "p", 
        pch = "<")
    legend("topleft", legend = 1:narrays, fill = cols, cex = 0.5)
    matplot(matrix(rep(forind1 * increm, narrays), ncol = narrays) - 
        rep(revmaxinds1 * increm, each = length(forind1)), tvScore2[forind2, 
        ], xlab = "Window Length", ylab = "TV Score", main = "TV Scores For  Different Windows Extending From Optimal Round 1 Window", 
        cex = 0.8, col = cols, type = "p", pch = ">", xlim = xlimuse, 
        ylim = ylimuse)
    matpoints(matrix(rep(-forind1 * increm, narrays), ncol = narrays) + 
        rep(formaxinds1 * increm, each = length(forind1)), tvScore2[revind2, 
        ], xlab = "Window Length", ylab = "TV Score", cex = 0.8, 
        col = cols, type = "p", pch = "<")
    legend("topleft", legend = 1:narrays, fill = cols, cex = 0.5)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
