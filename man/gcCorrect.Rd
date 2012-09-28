\name{gcCorrect}
\alias{gcCorrect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gcCorrect
}
\description{
methods to call gcCorrectMain for array data with varying formats
}
\usage{
gcCorrect(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
The array data 
}
  \item{\dots}{
all other arguments needed by gcCorrectMain
}
}
\details{
see gcCorrectMain
}
\value{
The default is to return corrected values for the arrays passed in. This
result will have the same dimensions as the first argument (the
uncorrected values).if returnOnlyTV is TRUE this will be a matrix of tv score values, one row per window, one column per array. The rownames of the matrix will indicate the size of the window. If onlyGC=TRUE the gc fraction of each window will be returned.
}
\references{
  Benjamini Y, Speed TP. Summarizing and cor\
  recting the GC content bias in high-throughput sequencing. Nucleic
  Acids Res 2012;40:e72
}
\author{
 Eitan Halper-Stromberg
}
\seealso{
gcCorrectMain
}
\examples{
nim2.1rda=system.file("extdata","nim2.1example.rda",package="ArrayTV")
load(nim2.1rda)

starts=loc-29.5 ## We want our start locations to be at the beginning of the probe, not in the middle of the probe
samplechr=c('chr14','chr15')
increms=c(10,1000,100e3)
wins=c(100,10e3,1e6)
initnodes=1
nodes=as.numeric(readline('How many cpus do you have available? (Enter the number)'))
nodes=max(c(initnodes,nodes),na.rm=TRUE)
if(nodes>1){
library(doMC)
registerDoMC(nodes)
}
tvScores=gcCorrect(Ms,chr, starts,samplechr,nodes,increms,wins,
jittercorrection=FALSE,returnOnlyTV=TRUE,onlyGC=FALSE,providedGC=0,'hg18',verbose=TRUE)

tvScoresSplit=split(data.frame(tvScores),rep(1:3,each=10))

print(tvScoresSplit)

dev.new()
par(mfrow=c(1,3))
for(ii in 1:3){
    plot(rownames(tvScoresSplit[[ii]]),tvScoresSplit[[ii]][,1],col='blue',ylim=c(0,max(tvScoresSplit[[ii]])),xlab='window extension from center',
         ylab='tv score')
    points(rownames(tvScoresSplit[[ii]]),tvScoresSplit[[ii]][,2],col='red')
}

### Now correct TVs based on best: we'll do it twice with one medium sized window (9000) and one large window (600,000)

cM1=gcCorrect(Ms,chr,starts,samplechr,nodes,6e5,6e5,
jittercorrection=FALSE,returnOnlyTV=FALSE,onlyGC=FALSE,providedGC=0,'hg18',verbose=TRUE)

cM2=gcCorrect(cM1,chr,starts,samplechr,nodes,9e3,9e3,
jittercorrection=FALSE,returnOnlyTV=FALSE,onlyGC=FALSE,providedGC=0,'hg18',verbose=TRUE)




### view results Whole Chromosome 14
dev.new()
par(mfrow=c(2,1))
## first array
plot(starts[chr=='chr14'],Ms[chr=='chr14',1],ylim=c(-.7,.7),main='Array 1 Uncorrected chr14',xlab='position',ylab='Signal intensity')
plot(starts[chr=='chr14'],cM2[chr=='chr14',1],ylim=c(-.7,.7),main='Array 1 Corrected chr14',xlab='position',ylab='Signal intensity')
dev.new()
par(mfrow=c(2,1))
## second array
plot(starts[chr=='chr14'],Ms[chr=='chr14',2],ylim=c(-.7,.7),main='Array 2 Uncorrected chr14',xlab='position',ylab='Signal intensity')
plot(starts[chr=='chr14'],cM2[chr=='chr14',2],ylim=c(-.7,.7),main='Array 2 Corrected chr14',xlab='position',ylab='Signal intensity')

### view Whole Chromosome 15
dev.new()
par(mfrow=c(2,1))
## first array
plot(starts[chr=='chr15'],Ms[chr=='chr15',1],ylim=c(-.7,.7),main='Array 1 Uncorrected chr15',xlab='position',ylab='Signal intensity')
plot(starts[chr=='chr15'],cM2[chr=='chr15',1],ylim=c(-.7,.7),main='Array 1 Corrected chr15',xlab='position',ylab='Signal intensity')
dev.new()
par(mfrow=c(2,1))
## second array
plot(starts[chr=='chr15'],Ms[chr=='chr15',2],ylim=c(-.7,.7),main='Array 2 Uncorrected chr15',xlab='position',ylab='Signal intensity')
plot(starts[chr=='chr15'],cM2[chr=='chr15',2],ylim=c(-.7,.7),main='Array 2 Corrected chr15',xlab='position',ylab='Signal intensity')


## view cbs line through a spike in (expected copy number 2.5)
cnaGR1=list();cnaGR2=list()
for(ii in 1:2){
    MsUse=get(c('Ms','cM2')[ii])
cna1=CNA(MsUse,chr,starts,data.type="logratio",sampleid=c(1,2))
smoothcna1=smooth.CNA(cna1)
segment.smoothed.CNA.object=segment(smoothcna1,verbose=1)
cnaGR=GRanges(seqnames=as.character(segment.smoothed.CNA.object$output$chrom),ranges=IRanges(start=segment.smoothed.CNA.object$output$loc.start,
end=segment.smoothed.CNA.object$output$loc.end))
ind1=which(as.character(segment.smoothed.CNA.object$output$ID)=='X1')
ind2=which(as.character(segment.smoothed.CNA.object$output$ID)=='X2')
cnaGR1[[ii]]=cnaGR[ind1]
cnaGR2[[ii]]=cnaGR[ind2]

elementMetadata(cnaGR1[[ii]])=segment.smoothed.CNA.object$output$seg.mean[ind1]
elementMetadata(cnaGR2[[ii]])=segment.smoothed.CNA.object$output$seg.mean[ind2]
}

## Plot Region near Spike In, the SpikeIn may be seen in the middle of the window as a small amplification
for(ii in 1:2){
    cnaGRuse=get(c('cnaGR1','cnaGR2')[ii])
dev.new()
par(mfrow=c(2,1))
spikeinStart=45396506;spikeinEnd=45537976
## uncorrected near spike in
plot(starts[chr=='chr15' & starts>(spikeinStart-3e6) & starts<(spikeinEnd+3e6)],
            Ms[chr=='chr15' & starts>(spikeinStart-3e6) & starts<(spikeinEnd+3e6),ii],
     main=paste('Uncorrected Array',ii), xlab='position',ylab='signal intensity',ylim=c(-.7,.5))
cbsXs=as.vector(rbind(start(cnaGRuse[[1]][seqnames(cnaGRuse[[1]])=='chr15']),end(cnaGRuse[[1]][seqnames(cnaGRuse[[1]])=='chr15']),NA))
cbsYs=rep(elementMetadata(cnaGRuse[[1]])[seqnames(cnaGRuse[[1]])=='chr15',1],each=3)
points(cbsXs,cbsYs,col='red',type='l')
## corrected near spikein
plot(starts[chr=='chr15' & starts>(spikeinStart-3e6) & starts<(spikeinEnd+3e6)],
            cM2[chr=='chr15' & starts>(spikeinStart-3e6) & starts<(spikeinEnd+3e6),ii],
     main=paste('Corrected Array',ii),xlab='position',ylab='signal intensity',ylim=c(-.7,.5))
cbsXs=as.vector(rbind(start(cnaGRuse[[2]][seqnames(cnaGRuse[[2]])=='chr15']),end(cnaGRuse[[2]][seqnames(cnaGRuse[[2]])=='chr15']),NA))
cbsYs=rep(elementMetadata(cnaGRuse[[2]])[seqnames(cnaGRuse[[2]])=='chr15',1],each=3)

points(cbsXs,cbsYs,col='red',type='l')
}
}

