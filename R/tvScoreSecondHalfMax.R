tvScoreSecondHalfMax <-
function(tvScore){
  tvScoreReverses=as(tvScore[((nrow(tvScore)/2)+1):nrow(tvScore),],"matrix")
  revmaxs=apply(tvScoreReverses,2,which.max)
  revmaxvals=apply(tvScoreReverses,2,max)  
  cbind(revmaxs,revmaxvals)
}
