priorFracs <-
function(gcFracBoth,maxuse,vv,nparts,v1,v2){
optnum=4
    extras=maxuse/(nparts/optnum)
    gcAdd=0
    if(extras>1){
      mo=maxuse %% (nparts/optnum)
      if(mo==0){
        gcAdd=gcFracBoth[[v2[-vv]]]
      }else{  
      gcAdd=gcFracBoth[[v2[-vv]]]
      mods=(1:length(gcAdd)) %% (nparts/optnum) 
      gcAdd[mods > mo | mods ==0]=0
    }
    }
     gcFracUse=gcFracBoth[[v1[-vv]]]+gcAdd
      
    priorGC=gcFracUse
    offset=ifelse(maxuse > (nparts/optnum),nparts/optnum,maxuse)
    priorFrac=(cumsum(priorGC)[offset:length(priorGC)]-c(0,cumsum(priorGC)[1:(length(priorGC)-offset)]))[seq(1,(length(priorGC)-offset+1),nparts/optnum)]/maxuse
  }
