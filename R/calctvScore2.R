calctvScore2 <- function(useM, gcFracBoth, tvScore, nparts, narrays){
	tvScore2 <- foreach(vv=1:2, .combine='rbind', .packages="ArrayTV") %:%
		foreach(ar=seq_along(narrays), .combine='cbind') %dopar% {
			print(vv)
			print(ar)
			tvScorepart <- as(optimalwin(useM,gcFracBoth,tvScore,vv,ar,nparts),"matrix")
		}
	return(tvScore2)
}
