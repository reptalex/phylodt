evolveD <- function(D,Active_Strains,dt,rate=1){
  x <- sp(Active_Strains)
  D[x,]=D[x,]+rate*dt
  D[,x]=D[,x]+rate*dt
  diag(D) <- 0
  return(D)
}
