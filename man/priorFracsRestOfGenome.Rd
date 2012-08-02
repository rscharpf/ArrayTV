\name{priorFracsRestOfGenome}
\alias{priorFracsRestOfGenome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{priorFracsRestOfGenome
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
priorFracsRestOfGenome(forwardExtend, reverseExtend, chrToInvestigate, locs, chrs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{forwardExtend}{
%%     ~~Describe \code{forwardExtend} here~~
}
  \item{reverseExtend}{
%%     ~~Describe \code{reverseExtend} here~~
}
  \item{chrToInvestigate}{
%%     ~~Describe \code{chrToInvestigate} here~~
}
  \item{locs}{
%%     ~~Describe \code{locs} here~~
}
  \item{chrs}{
%%     ~~Describe \code{chrs} here~~
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
function (forwardExtend, reverseExtend, chrToInvestigate, locs, 
    chrs) 
{
    gcFrac = rep(0, length(locs))
    winlen = reverseExtend + forwardExtend
    print(paste("calculating gc content for the rest of the genome for window", 
        forwardExtend, "bp"))
    for (x in chrToInvestigate) {
        li = which(chrs == x)
        ll = locs[li]
        schr = ifelse(grep("chr", x), x, paste("chr", x, sep = ""))
        schr = gsub("23", "X", schr)
        schr = gsub("24", "Y", schr)
        startinds = ll - reverseExtend
        startinds[startinds < 1] = 1
        startinds[startinds > (length(Hsapiens[[schr]]) - winlen)] = (length(Hsapiens[[schr]]) - 
            winlen)
        if (diff(range(startinds)) >= winlen) {
            gcIndivChr = letterFrequencyInSlidingView(unmasked(Hsapiens[[schr]]), 
                view.width = winlen, "CG", as.prob = T)[startinds]
            gcIndivChr[is.na(gcIndivChr)] = tail(gcIndivChr, 
                1)
        }
        else {
            gcIndivChr = rep(letterFrequencyInSlidingView(unmasked(Hsapiens[[schr]]), 
                view.width = length(Hsapiens[[schr]]), "CG", 
                as.prob = T), length(startinds))
        }
        cat(".")
        if (length(gcFrac[li]) != length(ll)) {
            print("something is wrong with gc dnastring set")
        }
        quants = sort(unique(quantile(gcIndivChr, probs = seq(0, 
            1, 0.01))))
        gcFrac[li] = quants[findInterval(gcIndivChr, quants)]
    }
    gcFrac
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
