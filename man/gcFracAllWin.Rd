\name{gcFracAllWin}
\alias{gcFracAllWin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate gc Fraction in windows of size increment
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gcFracAllWin(maxwin, increm, chr, allstarts, samplechr, uniqchrs, strategyuse)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxwin}{
%%     ~~Describe \code{maxwin} here~~
}
  \item{increm}{
%%     ~~Describe \code{increm} here~~
}
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{allstarts}{
%%     ~~Describe \code{allstarts} here~~
}
  \item{samplechr}{
%%     ~~Describe \code{samplechr} here~~
}
  \item{uniqchrs}{
%%     ~~Describe \code{uniqchrs} here~~
}
  \item{strategyuse}{
%%     ~~Describe \code{strategyuse} here~~
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
function (maxwin, increm, chr, allstarts, samplechr, uniqchrs, 
    strategyuse) 
{
    gcFracBoth = foreach(vv = 1:4, .combine = list, .multicombine = T) \%:\% 
        foreach(schr = uniqchrs[uniqchrs \%in\% samplechr], .combine = c) \%dopar\% 
        {
            starts = allstarts[chr \%in\% schr]
            if (strategyuse == 1) {
                print(paste("vv is", vv))
                if (vv == 1) {
                  ranges = IRanges(start = rep(starts, each = maxwin/increm/2) + 
                    rep(seq(0, (maxwin/2) - increm, increm), 
                      length(starts)), width = increm)
                }
                else if (vv == 2) {
                  ranges = IRanges(start = rep(starts, each = maxwin/increm/2) + 
                    rep(seq((maxwin/2), maxwin - increm, increm), 
                      length(starts)), width = increm)
                }
                else if (vv == 3) {
                  ranges = IRanges(start = rep(starts, each = maxwin/increm/2) - 
                    rep(seq(increm, maxwin/2, increm), length(starts)), 
                    width = increm)
                }
                else if (vv == 4) {
                  ranges = IRanges(start = rep(starts, each = maxwin/increm/2) - 
                    rep(seq((maxwin/2) + increm, maxwin, increm), 
                      length(starts)), width = increm)
                }
                ranges[which(start(ranges) < 1)] = IRanges(start = 1, 
                  end = increm)
                print("before All views")
                inDNAset = DNAStringSet()
                preinDNAset2 = DNAStringSet(Views(Hsapiens[[schr]], 
                  ranges))
                preinDNAset2[which(width(preinDNAset2) < increm)] = DNAStringSet((paste(rep("N", 
                  increm), collapse = "")))
                inDNAset = unlist(preinDNAset2)
                dnasetLen = length(inDNAset)
                gcFrac = letterFrequencyInSlidingView(inDNAset, 
                  view.width = increm, "CG", as.prob = T)[seq(1, 
                  dnasetLen, increm)]
            }
            if (strategyuse == 2) {
                pregcFrac = letterFrequencyInSlidingView(unmasked(Hsapiens[[schr]]), 
                  view.width = increm, "CG", as.prob = T)
                if (vv == 1) {
                  startinds = rep(starts, each = maxwin/increm/2) + 
                    rep(seq(0, (maxwin/2) - increm, increm), 
                      length(starts))
                }
                else if (vv == 2) {
                  startinds = rep(starts, each = maxwin/increm/2) + 
                    rep(seq((maxwin/2), maxwin - increm, increm), 
                      length(starts))
                }
                else if (vv == 3) {
                  startinds = rep(starts, each = maxwin/increm/2) - 
                    rep(seq(increm, maxwin/2, increm), length(starts))
                }
                else if (vv == 4) {
                  startinds = rep(starts, each = maxwin/increm/2) - 
                    rep(seq((maxwin/2) + increm, maxwin, increm), 
                      length(starts))
                }
                print(paste("gc window chunk", vv))
                startinds[which(startinds < 1)] = 1
                startinds[which(startinds > length(pregcFrac))] = length(pregcFrac)
                gcFrac = pregcFrac[startinds]
                gcFrac[is.na(gcFrac)] = 0
            }
            gcFrac
        }
    return(gcFracBoth)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
