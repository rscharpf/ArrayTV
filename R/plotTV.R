plotTV <-
function(tvScore,tvScore2,narrays,increm){
colinit=brewer.pal(8,'Set1')
pal <- colorRampPalette(colinit)
cols=pal(narrays)

## First Half of tvScore is extending window in forward direction
forind1=1:(nrow(tvScore)/2)
## Second Half of tvScore is extending window in reverse direciton
revind1=((nrow(tvScore)/2)+1):nrow(tvScore)
## The forward and reverse extension is flipped in tvScore2
forind2=((nrow(tvScore2)/2)+1):nrow(tvScore2)
revind2=1:(nrow(tvScore2)/2)


formaxinds1=tvScoreFirstHalfMax(tvScore)[,1]
revmaxinds1=tvScoreSecondHalfMax(tvScore)[,1]

#dev.new()

par(mfrow=c(2,1))
  xlimuse=c(min(-forind1*increm),max(forind1*increm))
  ylimuse=(c(0,max(c(tvScore2,tvScore))))
  matplot(forind1*increm,tvScore[forind1,],xlab='Window Length',ylab='TV Score',cex=.8,col=cols,type='p',pch='>',main='TV Scores For Different Windows Extending In Both Directions From 0',
       xlim=xlimuse,ylim=ylimuse)
  matpoints(-forind1*increm,tvScore[revind1,],xlab='Window Length',ylab='TV Score',cex=.8,col=cols,type='p',pch='<')

legend('topleft',legend=1:narrays,fill=cols,cex=.5)

  matplot(matrix(rep(forind1*increm,narrays),ncol=narrays)-rep(revmaxinds1 * increm,each=length(forind1)),tvScore2[forind2,],xlab='Window Length',ylab='TV Score',
  main='TV Scores For  Different Windows Extending From Optimal Round 1 Window',        cex=.8,col=cols,type='p',pch='>',xlim=xlimuse,ylim=ylimuse)
  matpoints(matrix(rep(-forind1*increm,narrays),ncol=narrays)+rep(formaxinds1 * increm,each=length(forind1)),tvScore2[revind2,],xlab='Window Length',ylab='TV Score',cex=.8,col=cols,type='p',pch='<')

legend('topleft',legend=1:narrays,fill=cols,cex=.5)
#abline(v=c(forind1*increm,-revind1*increm,forind2*increm - revmaxinds1[jj]*increm,-revind2*increm +formaxinds1[jj]*increm)[which.max(c(tvScore[forind1,jj],tvScore[revind1,jj],tvScore2[forind2,jj],tvScore2[revind2,jj]))],col=cols[jj])
}
