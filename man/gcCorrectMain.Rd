\name{gcCorrectMain}
\alias{gcCorrectMain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gcCorrectMain
}
\description{
Correct gc biases in array data. The user must provide a range of windows and this function will compute the tv score for each window and correct the arrays based upon the window with the highest tv score (i.e with window showing the most gc bias)
}
\usage{
   gcCorrectMain(Ms, chr, starts, samplechr, nodes, increms,
                          maxwins, jittercorrection=FALSE,
                          returnOnlyTV=FALSE, onlyGC=FALSE,providedGC=0,build, verbose=FALSE)
}
\arguments{
  \item{Ms}{
A matrix holding the array values
}
  \item{chr}{
a vector indicating the chromosome of each probe (should be the same length as the number of rows in the array matrix)
}
  \item{starts}{
a vector holding the genomic coordinates of the start positions the probes (should be the same length as the number of rows inthe array matrix)
}
  \item{samplechr}{
a subset (or the entire set) of the chromosomes represented in the arrays. The probes on these chromosomes will be used to compute the tvScore
}
  \item{nodes}{
the number of cpus available (more cpus=faster results)
}
  \item{increms}{
a vector of integers
}
  \item{maxwins}{
a vector ofintegers with each value >= the corresponding value in increms. Corresponding values of maxwins and increms must have the same quotient. For instance, if increms is c(10,1000) then an acceptable value of maxwins would be(500,5000) because the quotients of the corresponding elements are both 5.
}
\item{jittercorrection}{
if TRUE amplify the correction values by a small amount to decrease the autocorrelation of proximal probes
}
\item{returnOnlyTV}{
if TRUE only return the tvScores rather than applying corrections
}
\item{onlyGC}{
if TRUE return the gc fraction of each window (the result may be subsequently passed into the providedGC argument to be used in correction)
}
\item{providedGC}{
a vector of gc values (fractions from zero to 1) to use for computing corrections. This should be the same length as the number of rows in the array, if provided
}
 \item{build}{
the build matching the array. This is used to load the BSgenome object. If the arrays are annotated based on hg18 this would be 'hg18', if the annotation is based on hg19 this would be 'hg19'.
}
\item{verbose}{
provide output throughout computation
}
}
\details{
see example
}
\value{
The default is to return corrected values for the arrays passed in. This result will have the same dimensions as the first argument (the uncorrected values).if returnOnlyTV is TRUE this will be a matrix of tv score values, one row per window, one column per array. The rownames of the matrix will indicate the size of the window. If onlyGC=TRUE the gc fraction of each window will be returned.
}
\references{
Optimal window finding was inspired by a similar concept applied to gc correction for next generation sequencing experiments:Benjamini Y, Speed TP. Summarizing and correcting the GC content bias in high-throughput sequencing. Nucleic Acids Res 2012;40:e72
}
\author{
Eitan Halper-Stromberg
}
\examples{
nim2.1rda <- system.file("extdata","nim2.1example.rda",package="ArrayTV")
load(nim2.1rda)

starts=loc-29.5 ## We want our start locations to be at the beginning of the probe, not in the middle of the probe
samplechr=c('chr14','chr15')
increms=c(10,1000,100e3)
wins=c(100,10e3,1e6)
initnodes=1
nodes=as.numeric(readline('How many cpus do you have available? (Enter
the number)'))
nodes=max(c(initnodes,nodes),na.rm=TRUE)
if(nodes>1){
library(doMC)
registerDoMC(nodes)
}
tvScores=gcCorrectMain(Ms,chr, starts,samplechr,nodes,increms,wins,
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

cM1=gcCorrectMain(Ms,chr,starts,samplechr,nodes,6e5,6e5,
jittercorrection=FALSE,returnOnlyTV=FALSE,onlyGC=FALSE,providedGC=0,'hg18',verbose=TRUE)

cM2=gcCorrectMain(cM1,chr,starts,samplechr,nodes,9e3,9e3,
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
