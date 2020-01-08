plot.phylodt <- function(x,legend.position='topleft',edge.alpha=0.2,...){
  par(mfrow=c(2,1))
  mx <- max(x$SI_dynamics$SI)
  mn <- min(x$SI_dynamics$SI)
  plot(x$SI_dynamics$time,x$SI_dynamics$SI[,'S'],col='black',type = 'l',lwd=2,ylim=c(0.5*mn,1.1*mx),
       xlab='time',ylab='# Individuals',main='SI dynamics')
  lines(x$SI_dynamics$time,x$SI_dynamics$SI[,'I'],type = 'l',pch=NA,col='red',lwd=2)
  legend(x=legend.position,lty=c(1,1),col=c('black','red'),legend = c('S','I'),lwd=c(2,2))
  
  plot(x$tree,show.tip.label = F,edge.color = rgb(0,0,0,edge.alpha))
}
