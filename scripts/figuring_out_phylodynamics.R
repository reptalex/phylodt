library(ape)
library(data.table)
library(magrittr)

N=10
Population <- data.table('Individual'=1:N,
                         'Status'=sample(c('S','I'),N,T),
                         'history'=c('2','1;3','4',rep(0,N-3)))
Population$strain <- as.integer(0)

SI <- table(Population$Status)[c('S','I')]
Population[Status=='I']$strain <- 5:(4+sum(Population$Status=='I'))
# tr <- rtree(N,tip.label = 1:N)
Active_Strains <- unique(Population[strain>0,strain])


# how to grow a tree ------------------------------------------------------

###### AND efficiently change distance matrix
# 
# evolve <- function(tree,dt.=dt,Active_Strains.=Active_Strains,rate=1){ #grow all active tips
#   growing.edges <- match(Active_Strains,tree$edge[,2])
#   tree$edge.length[growing.edges] <- tree$edge.length[growing.edges]+rate*dt
#   return(tree)
# }
# evolveD <- function(D,dt.=dt,Active_Strains.=Active_Strains,rate=1){
#   x <- sp(Active_Strains)
#   D[x,]=D[x,]+rate*dt
#   D[,x]=D[,x]+rate*dt
#   diag(D) <- 0
#   return(D)
# }
# 
# speciate <- function(tree,splitting.tip.=splitting.tip){ #split tip
#   N=length(tree$tip.label)
#   x <- max(tree$edge[tree$edge>N])+2 #number of our new node
#   tree$edge[tree$edge>N] <- tree$edge[tree$edge>N]+1
#   edge <- which(tree$edge[,2]==splitting.tip)
#   tree$edge[edge,2] <- x
#   tree$edge <- rbind(tree$edge,matrix(c(x,splitting.tip,x,N+1),byrow=T,nrow=2,ncol=2))
#   tree$tip.label <- c(tree$tip.label,paste('s',N+1,sep=''))
#   tree$Nnode <- tree$Nnode+1
#   tree$edge.length <- c(tree$edge.length,0,0)
#   return(tree)
# }

# sp <- function(x) paste('s',x,sep='')
# 
# 
# # initialize
# tree <- NULL
# tree$edge <- matrix(c(3,1,3,2),nrow=2,ncol=2,byrow=T)
# tree$tip.label <- c('s1','s2')
# tree$edge.length <- c(0,0)
# tree$Nnode <- 1
# class(tree) <- 'phylo'
# plot(tree)
# 
# Active_Strains <- c(1,2)
# 
# D <- matrix(0,nrow=2,ncol=2)
# rownames(D) <- paste('s',1:2,sep='')
# colnames(D) <- rownames(D)
# 
# ## Assume strain 1 goes extinct in 5 days
# dt=5
# tree <- evolve(tree)
# D <- evolveD(D)
# Active_Strains <- setdiff(Active_Strains,1)
# 
# 
# ## Now assume strain 2 speciates in 1 day
# dt=1
# splitting.tip=2
# ix=sp(splitting.tip)
# D <- evolveD(D)
# tree <- evolve(tree) %>% speciate
# D <- rbind(D,D[ix,])
# D <- cbind(D,rbind(D[,ix,drop=F]))
# rownames(D)[nrow(D)] <- data.table::last(tree$tip.label)
# colnames(D)[ncol(D)] <- data.table::last(tree$tip.label)
# Active_Strains <- c(Active_Strains,tree$edge[nrow(tree$edge),2])
# 
# 
# ## Now assume strain 2 dies in 1 day
# dt=1
# tree <- evolve(tree)
# D <- evolveD(D)
# Active_Strains <- setdiff(Active_Strains,2)
# D <- D[sp(Active_Strains),sp(Active_Strains),drop=F]
# 
# ## Species 3 speciates in 3 days
# dt=3
# splitting.tip <- 3
# ix=sp(splitting.tip)
# D <- evolveD(D)
# tree <- evolve(tree) %>% speciate
# D <- rbind(D,D[ix,])
# D <- cbind(D,rbind(D[,ix,drop=F]))
# rownames(D)[nrow(D)] <- data.table::last(tree$tip.label)
# colnames(D)[ncol(D)] <- data.table::last(tree$tip.label)
# Active_Strains <- c(Active_Strains,tree$edge[nrow(tree$edge),2])
# 
# # Species 3 speciates in 1 day
# 
# dt=1
# splitting.tip <- 3
# ix=paste('s',splitting.tip,sep='')
# D <- evolveD(D)
# tree <- evolve(tree) %>% speciate
# D <- rbind(D,D[ix,])
# D <-  cbind(D,rbind(D[,ix,drop=F]))
# rownames(D)[nrow(D)] <- data.table::last(tree$tip.label)
# colnames(D)[ncol(D)] <- data.table::last(tree$tip.label)
# Active_Strains <- c(Active_Strains,tree$edge[nrow(tree$edge),2])
# 
# # Species 5 goes extinct in 0.5 day
# dt=0.5
# tree <- evolve(tree)
# D <- evolveD(D)
# Active_Strains <- setdiff(Active_Strains,5)
# D <- D[sp(Active_Strains),sp(Active_Strains),drop=F]


######### NOTE: we'll need to include species with history in Population$history
#########       in our distance matrix, D



# Gillespie SIR ----------------------------------------------------------
rm(list=ls())
fls <- list.files('R/')
sapply(paste('R/',fls,sep=''),source)

set.seed(1)
N_avg=200 ## average population size at equilibrium without infections
mu=(1/(50*365)) ## life expectancy of 50 years
lambda=N_avg*mu
# lambda=0.01  ## rate of migration of new susceptibles in days - on average 100 days
nu <- 1/7 # clearance rate
I0=3
d0=1
R0=2.6
bmax=R0*(nu+mu)/(N_avg-I0)
Tmax <- 50

beta <- function(d,bmax.=bmax,r=log(2)/(30*365)) bmax*(1-exp(-r*d)) #default memory half-life of 15 years
sp <- function(x) paste('s',x,sep='')
h <- function(x) paste('h',x,sep='')


########### INITIALIZE
Population <- data.table('Individual'=1:N_avg,
                         'Status'=c(rep('I',I0),rep('S',N_avg-I0)),
                         'strain'=c(1:I0,rep(0,N_avg-I0)))
SI <- table(Population$Status)[c('S','I')]
SI
N=sum(SI)
output<- list('time'=0,'SI'=SI)

histories <- vector(mode='list',length=N_avg)
names(histories) <- h(1:N_avg)
histories[1:I0] <- as.list(1:I0)

tree <- rtree(I0)
tree$tip.label <- sp(1:I0)
# tree <- phytools::starTree(sp(1:I0)) # initialize w/ star tree
# tree$edge.length <- rep(d0,nrow(tree$edge))
nodebank <- makeNodeMap(tree)
Active_Strains <- 1:I0         #strains still evolving - determines +dt updates of D
Memory_Strains <- 1:I0         #strains still in serological memory - determines rows/cols of D
D <- matrix(d0,nrow=I0,ncol=I0)
diag(D) <- 0
rownames(D) <- sp(1:I0)
colnames(D) <- rownames(D)
B <- matrix(Inf,nrow=I0,ncol=N-I0) %>% beta
rownames(B) <- sp(1:I0)
colnames(B) <- h((I0+1):N)
## need HistoryMap(strain,S,D) which outputs distance from S history to strain.
HistoryMap <- B
HistoryMap[HistoryMap!=0] <- 0
# HistoryMap['s2','h12'] will be the strain k 'sk' nearest to the history of 'h12' and 's2'
## updates will be: 
## infection --> add row, lose column
## recovery --> lose row, add column
## No element of HistoryMap will ever be a row/circulating strain (S don't have circulating strains).
## ==> elements of HistoryMap can be searched efficiently by focusing on setdiff(Memory_Strains,Active_Strains)
### For the time being, we assume the tips where strains went extinct are 
## suitable approximate place marker for an individual's serological memory
## as opposed to the basal node where the individual became infected.

Tmax=300
tt <- 0
events <- c('Birth','Death','Recovery','Infection')
propensities <- rep(NA,length(events))
names(propensities) <- events

while(tt<=Tmax){
  propensities['Birth'] <- lambda
  propensities['Death'] <- mu*sum(SI)
  propensities['Recovery'] <- nu*SI['I']
  b=rowSums(B)
  propensities['Infection'] <- sum(b)
  
  ## tau=1 #units: time 
  Rate=sum(propensities) #units: per-time
  ## N_events=max(rpois(1,tau*Rate),1)
  event <- sample(events,1,prob=propensities)  ## insert tau-leaping here
  
  ### evolve
  dt <- rexp(1,Rate)
  tree <- evolve(tree,Active_Strains,dt)
  D <- evolveD(D,Active_Strains,dt)
  tt=tt+dt
  ### update
  if (event=='Birth'){
    
    ################################### BIRTH #####################################
    Population <- rbind(Population,
                        data.table('Individual'=max(Population$Individual)+1,
                                   'Status'='S','strain'=0))
    x=Population[.N,Individual]
    HistoryMap <- cbind(HistoryMap,rep(0,nrow(HistoryMap)))
    colnames(HistoryMap)[ncol(HistoryMap)] <- h(x)
    B <- cbind(B,rep(bmax,nrow(B)))
    colnames(B)[ncol(B)] <- h(x)
    
    SI['S']=SI['S']+1
    N=N+1
  } else if (event=='Death'){
    
    ################################### DEATH #####################################
    ix <- sample(nrow(Population),1)
    x=Population[ix,Individual]
    is_infected <- Population[ix,Status]=='I'
    
    if (is_infected){
      ## eradicate unique serological history & current strain
      ## must update historymap, B
      strain <- Population[ix,strain]
      forgotten_strains <- setdiff(unique(c(histories[[h(x)]],strain)),
              unique(unlist(histories[setdiff(names(histories),h(x))])))
      Memory_Strains <- setdiff(Memory_Strains,forgotten_strains)
      Active_Strains <- setdiff(Active_Strains,strain)
      
      histories <- histories[setdiff(names(histories),h(x))]
      HistoryMap <- HistoryMap[sp(Memory_Strains),,drop=F]
      B <- B[sp(Memory_Strains),,drop=F]
      D <- D[sp(Memory_Strains),sp(Memory_Strains)]
      Population <- Population[setdiff(1:nrow(Population),ix)]
      
      SI['I']=SI['I']-1
      N=N-1
    } else {
      ## eradicate unique serological history
      forgotten_strains <- setdiff(histories[[h(x)]],
                                   unlist(histories[setdiff(names(histories),h(x))]))
      Memory_Strains <- setdiff(Memory_Strains,forgotten_strains)
      histories <- histories[setdiff(names(histories),h(x))]
      
      HistoryMap <- HistoryMap[intersect(rownames(HistoryMap),sp(Memory_Strains)),
                               setdiff(colnames(HistoryMap),h(x)),drop=F]
      B <- B[intersect(rownames(B),sp(Memory_Strains)),
             setdiff(colnames(B),h(x)),drop=F]
      D <- D[sp(Memory_Strains),sp(Memory_Strains)]
      
      SI['S']=SI['S']-1
      N=N-1
      Population <- Population[setdiff(1:nrow(Population),ix)]
    }
    
    
  } else if (event=='Recovery'){
    ################################### RECOVERY #####################################
    Is <- which(Population$Status=='I')
    if (length(Is)>1){
      ix <- sample(Is,1)
    } else {
      ix <- Is
    }
    strain <- Population[ix,strain]
    Active_Strains <- setdiff(Active_Strains,strain)
    Population[ix,Status:='S']
    Population[ix,strain:=0]
    x=Population[ix,Individual]
    histories[[h(x)]] <- setdiff(unique(c(histories[[h(x)]],strain)),0)
    strains <- unique(histories[[h(x)]])
    ## we remove this strain from HistoryMap's rows and add a column for our new susceptible's history
    ## and recompute Betas in light of our new tree=tree+dt
    HistoryMap <- HistoryMap[setdiff(rownames(HistoryMap),sp(strain)),,drop=F]
    B <- B[rownames(HistoryMap),,drop=F]
    HistoryMap <- cbind(HistoryMap,minDists(rownames(HistoryMap),sp(strains),D))
    colnames(HistoryMap)[ncol(HistoryMap)] <- h(x)
    B <- cbind(B,rep(0,nrow(B)))
    colnames(B)[ncol(B)] <- h(x)

    SI=SI+c(1,-1)
  } else {
    
    ################################### INFECTION #####################################
    # ix <- sample(Population[,which(Status=='S')],1)
    # x <- Population[ix,Individual]  ## who gets infected?
    x <- as.numeric(gsub('h','',sample(colnames(B),1,prob=colSums(B))))
    ix <- which(Population$Individual==x)
    splitting.tip <- as.numeric(gsub('s','',sample(rownames(B),1,prob=B[,h(x)])))
    # tree <- speciate(tree,splitting.tip)
    tree_nb <- speciate(tree,nodebank,splitting.tip)
    tree <- tree_nb$tree
    nodebank <- tree_nb$nodebank
    rm('tree_nb')
    new.strain <- as.numeric(gsub('s','',data.table::last(tree$tip.label)))
    Active_Strains <- c(Active_Strains,new.strain)
    Memory_Strains <- c(Memory_Strains,new.strain)         ## Is this unnecessary in light of B?
    histories[[h(x)]] <- c(histories[[h(x)]],new.strain)
    
    HistoryMap <- HistoryMap[,setdiff(colnames(HistoryMap),h(x)),drop=F]
    HistoryMap <- rbind(HistoryMap,HistoryMap[sp(splitting.tip),,drop=F])
    rownames(HistoryMap)[nrow(HistoryMap)] <- sp(new.strain)
    B <- B[,setdiff(colnames(HistoryMap),h(x)),drop=F]
    B <- rbind(B,B[sp(splitting.tip),,drop=F])
    rownames(B)[nrow(B)] <- sp(new.strain)
    if (any(is.na(B))){
      stop()
    }
    
    D <- rbind(D,D[sp(splitting.tip),,drop=F])
    D <- cbind(D,rbind(D[,sp(splitting.tip),drop=F]))
    rownames(D)[nrow(D)] <- sp(new.strain)
    colnames(D)[ncol(D)] <- sp(new.strain)
    Population[ix,Status:='I']
    Population[ix,strain:=new.strain]
    SI=SI+c(-1,1)
    
    ###############################################################################
  }
  
  output$time <- c(output$time,tt)
  output$SI <- rbind(output$SI,SI)
  if (any(HistoryMap!=0)){
    B[HistoryMap!=0] <- beta(computeDists(HistoryMap,D))
  }
  if (dim(B)[1]==0){
    return(list(''))
  }
}

x <- list('SI_dynamics'=output,
          'tree'=tree)
class(x) <- 'phylodt'

plot(x)
# plot.phylodt <- function(x,legend.position='topleft',...){
#   par(mfrow=c(2,1))
#   mx <- max(x$SI_dynamics$SI)
#   mn <- min(x$SI_dynamics$SI)
#   plot(x$SI_dynamics$time,x$SI_dynamics$SI[,'S'],col='black',type = 'l',lwd=2,ylim=c(0.5*mn,1.1*mx),
#        xlab='time',ylab='# Individuals',main='SI dynamics')
#   lines(x$SI_dynamics$time,x$SI_dynamics$SI[,'I'],type = 'l',pch=NA,col='red',lwd=2)
#   legend(x=legend.position,lty=c(1,1),col=c('black','red'),legend = c('S','I'),lwd=c(2,2))
#   
#   plot(x$tree,show.tip.label = F,edge.color = rgb(0,0,0,0.2))
# }



rm(list=ls())
gc()
list.files('R/') %>% sapply(FUN=function(x) source(paste('R/',x,sep='')))
x <- phylodt(N_avg = 200,R0 = 3,I0=2)
plot.phylodt(x,legend.position='top',edge.alpha=0.4)
nodelabels(pch=16)

PD <- pdDynamics(x$tree)

beta <- function(d,R0=1.1){
  d[]=R0*(1/7+1/(50*365))/(90) 
  return(d)
}
x2 <- phylodt(N_avg=200,beta = beta,I0=3)
plot.phylodt(x2,legend.position='right')

PD2 <- pdDynamics(x2$tree)


par(mfrow=c(1,1))
plot(PD2$time,PD2$avg_pd,type='l',log='y',col='blue',lwd=2,xlab='time (days)',ylab='Phylogenetic Diversity')
lines(PD$time,PD$avg_pd,type='l',col='black',lwd=2)
legend('bottomright',legend = c('Cross-Reactivity','Random'),lty=c(1,1),col=c('blue','black'),lwd=c(2,2))


par(mfrow=c(2,1))
d=ape::cophenetic.phylo(x2$tree)
hist(d[upper.tri(d)],breaks=20,main='Random infections',xlab='Phylogentic Distances')
d2=ape::cophenetic.phylo(x$tree)
hist(d2[upper.tri(d2)],breaks=20,col='blue',main='Cross-Reactivity',xlab='Phylogenetic Distances')


x3=phylodt(Tmax=200)
