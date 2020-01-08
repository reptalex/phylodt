phylodt <- function(Tmax=50,N_avg=100,mu=(1/(50*365)),nu=1/7,I0=10,d0=1,R0=1.6,lambda=NULL,bmax=NULL,beta=NULL,tree=NULL){
  # N_avg=100 ## average population size at equilibrium without infections
  # mu=(1/(50*365)) ## life expectancy of 50 years
  if (is.null(lambda)){
    lambda=N_avg*mu
  }
  # lambda=0.01  ## rate of migration of new susceptibles in days - on average 100 days
  # nu <- 1/7 # clearance rate
  # I0=10
  # d0=1
  # R0=1.6
  if (is.null(bmax)){
    bmax=R0*(nu+mu)/(N_avg-I0)
  }
  
  if(is.null(beta)){
    beta <- function(d,bmax.=bmax,r=log(2)/(15*365)) bmax*(1-exp(-r*d)) #default memory half-life of 15 years
  }
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
  
  if (is.null(tree)){
    tree <- ape::rtree(I0,tip.label = sp(1:I0)) # initialize w/ star tree
    tree <- phytools::force.ultrametric(tree)
    ## initialize with average 
    D=ape::cophenetic.phylo(tree)
    mn <- mean(D[upper.tri(D)])
    tree$edge.length <- d0*tree$edge.length/mn
    D=ape::cophenetic.phylo(tree)
    D <- D[sp(1:I0),sp(1:I0)]
  } else {
    I0=ape::Ntip(tree)
    D=ape::cophenetic.phylo(tree)
    D <- D[sp(1:I0),sp(1:I0)]
  }
  nodebank <- makeNodeMap(tree)
  Active_Strains <- 1:I0         #strains still evolving - determines +dt updates of D
  Memory_Strains <- 1:I0         #strains still in serological memory - determines rows/cols of D
  
  B <- matrix(Inf,nrow=I0,ncol=N_avg-I0) %>% beta
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
      strains_in_history <- unique(histories[[h(x)]])
      ## we remove this strain from HistoryMap's rows and add a column for our new susceptible's history
      ## and recompute Betas in light of our new tree=tree+dt
      HistoryMap <- HistoryMap[setdiff(rownames(HistoryMap),sp(strain)),,drop=F]
      B <- B[rownames(HistoryMap),,drop=F]
      HistoryMap <- cbind(HistoryMap,
                          minDists(rownames(HistoryMap),
                                   sp(strains_in_history),
                                   D))
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
      return(list('SI_dynamics'=output,
                  'tree'=tree))
    }
  }
  return(list('SI_dynamics'=output,
              'tree'=tree))
}
