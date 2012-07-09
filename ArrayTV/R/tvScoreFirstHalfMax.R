tvScoreFirstHalfMax <-
function(tvScore){
tvScoreForwards=tvScore[1:(nrow(tvScore)/2),]
formaxs=apply(tvScoreForwards,2,which.max)
formaxvals=apply(tvScoreForwards,2,max)
cbind(formaxs,formaxvals)
}
