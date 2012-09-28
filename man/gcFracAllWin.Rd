\name{gcFracAllWin}
\alias{gcFracAllWin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gcFracAllWin
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gcFracAllWin(maxwin, maxwin2, increm, increm2, chr, allstarts, samplechr, uniqchrs, strategyuse)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxwin}{
%%     ~~Describe \code{maxwin} here~~
}
  \item{maxwin2}{
%%     ~~Describe \code{maxwin2} here~~
}
  \item{increm}{
%%     ~~Describe \code{increm} here~~
}
  \item{increm2}{
%%     ~~Describe \code{increm2} here~~
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
function (maxwin, maxwin2, increm, increm2, chr, allstarts, samplechr, 
    uniqchrs, strategyuse) 
{
    gcFracBoth = foreach(schr = uniqchrs[uniqchrs \%in\% samplechr], 
        .combine = rbind) \%dopar\% {
        starts = allstarts[chr \%in\% schr]
        if (strategyuse == 1) {
            print(paste("vv is", vv))
            if (vv == 1) {
                ranges = IRanges(start = rep(starts, each = maxwin/increm/2) + 
                  rep(seq(0, (maxwin/2) - increm, increm), length(starts)), 
                  width = increm)
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
            gcFrac = letterFrequencyInSlidingView(inDNAset, view.width = increm, 
                "CG", as.prob = T)[seq(1, dnasetLen, increm)]
        }
        if (strategyuse == 2) {
            schr = gsub("23", "X", schr)
            schr = gsub("24", "X", schr)
            print("Getting gc content From BS genome Object")
            pregcFrac = letterFrequencyInSlidingView(unmasked(Hsapiens[[schr]]), 
                view.width = increm, "CG", as.prob = T)
            pregcFrac2 = letterFrequencyInSlidingView(unmasked(Hsapiens[[schr]]), 
                view.width = increm2, "CG", as.prob = T)
            print("gc content stored")
            startinds = rep(starts, each = maxwin/increm) + rep(seq(0, 
                maxwin - increm, increm), length(starts))
            startinds2 = rep(starts, each = maxwin2/increm2) + 
                rep(seq(0, maxwin2 - increm2, increm2), length(starts))
            startindbackwards = rep(starts, each = maxwin/increm) - 
                rep(seq(increm, maxwin, increm), length(starts))
            startindbackwards2 = rep(starts, each = maxwin2/increm2) - 
                rep(seq(increm2, maxwin2, increm2), length(starts))
            startinds[which(startinds < 1)] = 1
            startinds[which(startinds > length(pregcFrac))] = length(pregcFrac)
            gcFrac1 = pregcFrac[startinds]
            gcFrac1[is.na(gcFrac1)] = 0
            startinds2[which(startinds2 < 1)] = 1
            startinds2[which(startinds2 > length(pregcFrac2))] = length(pregcFrac2)
            gcFrac2 = pregcFrac2[startinds2]
            gcFrac2[is.na(gcFrac2)] = 0
            startindbackwards[which(startindbackwards < 1)] = 1
            startindbackwards[which(startindbackwards > length(pregcFrac))] = length(pregcFrac)
            gcFracbackwards1 = pregcFrac[startindbackwards]
            gcFracbackwards1[is.na(gcFracbackwards1)] = 0
            startindbackwards2[which(startindbackwards2 < 1)] = 1
            startindbackwards2[which(startindbackwards2 > length(pregcFrac2))] = length(pregcFrac2)
            gcFracbackwards2 = pregcFrac2[startindbackwards2]
            gcFracbackwards2[is.na(gcFracbackwards2)] = 0
            gcFrac = cbind((gcFrac1 + gcFracbackwards1)/2, (gcFrac2 + 
                gcFracbackwards2)/2)
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
