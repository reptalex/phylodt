library(phylofactor)
fls <- list.files('R/')
sapply(paste('R/',fls,sep=''),source)

set.seed(1)
tr <- rtree(10)
tr$tip.label <- gsub('t','s',tr$tip.label)
plot(tr)
nodelabels()

nodebank <- makeNodeMap(tr)

t_8=speciate(tr,nodebank,which(tr$tip.label=='s3'))
t_8$tree <- evolve(t_8$tree,dt=0.1)


t_5=speciate(t_8$tree,t_8$nodebank,5)
t_5$tree <- evolve(t_5$tree,dt=0.1)

par(mfrow=c(1,3))
plot(tr)
nodelabels()
plot(t_8$tree)
nodelabels()
plot(t_5$tree)
nodelabels()
