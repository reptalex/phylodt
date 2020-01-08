evolve <- function(tree,Active_Strains=NULL,dt,rate=1){ #grow all active tips
  if (is.null(Active_Strains)){
    Active_Strains <- 1:ape::Ntip(tree)
  }
  growing.edges <- match(Active_Strains,tree$edge[,2])
  tree$edge.length[growing.edges] <- tree$edge.length[growing.edges]+rate*dt
  return(tree)
}