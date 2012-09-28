tvScoreFirstHalfMax <-
function(tvScore){
tvScoreForwards=as(tvScore[1:(nrow(tvScore)/2),],"matrix")
formaxs=apply(tvScoreForwards,2,which.max)
formaxvals=apply(tvScoreForwards,2,max)
cbind(formaxs,formaxvals)
}
