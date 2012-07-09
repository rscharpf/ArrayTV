tvScoreSecondHalfMax <-
function(tvScore){
  tvScoreReverses=tvScore[((nrow(tvScore)/2)+1):nrow(tvScore),]
  revmaxs=apply(tvScoreReverses,2,which.max)
  revmaxvals=apply(tvScoreReverses,2,max)  
  cbind(revmaxs,revmaxvals)
}
