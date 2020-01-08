pdDynamics <- function(tree){
  D <- ape::node.depth.edgelength(tree)
  PD <- ape::dist.nodes(tree)
  tt <- sort(unique(D))
  n=ape::Ntip(tree)
  output <- data.table('time'=c(0,tt),'pd'=NA,'species'=NA)
  output$pd[1] <- 0
  output$species[1] <- 1
  no_descendants <- function(nds,tree){
    for (nd in nds){
      rp=setdiff(ape::nodepath(tree,nd,ape::Ntip(tree)+1),nd)
      nds <- setdiff(nds,rp)
    }
    return(nds)
  }
  interior_nodes <- NULL
  for (d in 1:length(tt)){
    ix <- which(D<=tt[d])
    pd <- PD[ix,ix,drop=F]
    possible_tips <- setdiff(ix,interior_nodes)
    lineages <- no_descendants(possible_tips,tree)
    interior_nodes <- unique(c(interior_nodes,setdiff(ix,lineages)))
    output$species[d+1] <- length(lineages)
    output$pd[d+1] <- sum(pd)
  }
  output[,avg_pd:=pd/species]
  return(output)
}
