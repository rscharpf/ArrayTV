test_gcCorrectMain <- function() {
    ## this is the same correction we apply in the vignette
    path <- system.file("extdata", package="ArrayTV")
    load(file.path(path, "array_logratios.rda"))
    nimblegen[, "M"] <- nimblegen[, "M"]/1000
    max.window <- c(100,10e3,1e6)
    increms <- c(20, 2000, 200e3)
    if(require(doParallel)){
        cl <- makeCluster(1)
        registerDoParallel(cl)
    }

    nimcM1List <- gcCorrect(object=nimblegen[, "M", drop=FALSE],
			chr=rep("chr15", nrow(nimblegen)),
                        starts=nimblegen[, "position"],
                        increms=increms,
                        maxwins=max.window,
                        build='hg18')


    ## Do the correction values follow gc content? plot these two to visualize
    checkTrue(cor.test(nimblegen[,"M"]-nimcM1List[['correctedVals']],nimcM1List[['GCvals']])$estimate > .8)
    ## Our corrected points should be closer together than uncorrected
    checkTrue(mad(nimblegen[,"M"])-mad(nimcM1List[['correctedVals']]) > .013)
    ## Our TVscore should be pretty large for this chromosome
    checkTrue(nimcM1List[["maxTVscore"]] > 0.05)
    save(nimcM1List, file=file.path(path,"nimcM1List.rda"))
    stopCluster(cl)
}
