getOrientation <- function(np,tree){
    db <- data.table('node'=np,'orientation'=NA)
    for (i in 2:length(np)){
      db$orientation[i]=c('pos','neg')[which(tree$edge[tree$edge[,1]==np[i],2]==np[i-1])]
    }
    return(db)
}