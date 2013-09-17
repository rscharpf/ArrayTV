gcFracOneRange<-function(schr,starts,increm,maxwin,verbose=FALSE){
    if(verbose) message('Getting gc content From BS genome Object')
    Hsapiens <- get("Hsapiens")
    pregcFrac <- letterFrequencyInSlidingView(unmasked(Hsapiens[[schr]]),view.width=increm,'CG',as.prob=T)
    if(verbose) message('gc content stored')
    startinds <- rep(starts, each=maxwin/increm) + rep(seq(0, maxwin-increm, increm), length(starts))
    startindbackwards <- rep(starts,each=maxwin/increm)-rep(seq(increm,maxwin,increm),length(starts))


    startinds[which(startinds<1)] <- 1
    startinds[which(startinds>length(pregcFrac))] <- length(pregcFrac)
    gcFrac1 <- pregcFrac[startinds]
    gcFrac1[is.na(gcFrac1)] <- 0

    startindbackwards[which(startindbackwards<1)] <- 1
    startindbackwards[which(startindbackwards>length(pregcFrac))] <- length(pregcFrac)
    gcFracbackwards1 <- pregcFrac[startindbackwards]
    gcFracbackwards1[is.na(gcFracbackwards1)] <- 0

    return(round((gcFrac1+gcFracbackwards1)/2,3))
}
