computeDists<-function(HistoryMap,D){
  ix <- which(HistoryMap!=0,arr.ind = T)
  if (nrow(ix)==0){
    return(NULL)
  } else {
    DF <- data.table('y'=rep(0,nrow(ix)),'x'=0)
    DF$y <- match(sp(HistoryMap[as.matrix(ix)]),colnames(D))
    DF$x <- match(rownames(HistoryMap)[ix[,'row']],colnames(D))
    return(D[as.matrix(DF)])
  }
}
