minDists <- function(x,y,D){
  if (length(y)==1){
    return(rep(as.numeric(gsub('s','',y)),length(x)))
  } else {
    return(as.numeric(gsub('s','',y[sapply(x,FUN=function(xx,y,D) which.min(D[xx,y]),y,D)])))
  }
}
