calctvScore2 <-
function(useM,gcFracBoth,tvScore,nparts,narrays){
tvScore2=foreach(vv=1:2,.combine='rbind') %:% 
  foreach(ar=1:narrays,.combine='cbind') %dopar% {
    print(vv)  
  print(ar)
  tvScorepart=optimalwin(useM,gcFracBoth,tvScore,vv,ar,nparts)
}
return(tvScore2)
}
