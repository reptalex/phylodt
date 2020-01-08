speciate <- function(tree,nodebank,splitting.tip,just_tree=F){ #split tip
  N=length(tree$tip.label)
  nds=tree$edge>N
  x=max(tree$edge[nds])+1 #number of our new node
  rp <- rootpath(tree,splitting.tip) %>% getOrientation(tree)
  nodeSeq <- function(nd,ng,p){
    if (ng==0){
      return(NULL)
    } else {
      return(seq(nd+p+1,nd+ng+p))
    }
  }
  
  if ('pos' %in% rp$orientation){ #we have some nodes on - side sister clades of splitting.tip
    neg_anc=nodebank[(node %in% rp[orientation=='pos']$node) & neg>0,]
    if (nrow(neg_anc)!=0){
      neg_nds= neg_anc %>% 
                apply(1,FUN=function(x) nodeSeq(x['node'],x['neg'],x['pos'])) %>%
                unlist %>% unique %>% c
      neg_nds <- setdiff(neg_nds,rp$node)
    } else {
      neg_nds <- NULL
    }
  } else {
    neg_nds <- NULL
  }
  
  ### there are three groups of nodes:
  ### 1) nodes on rootpath of our splitting.tip
  ### 2) nodes on + side sister clades of splitting.tip
  ### 3) nodes on -side sister clades of splitting.tip
  
  if (rp$orientation[2]=='pos'){
    new_node <- rp$node[2]+2
  } else {
    new_node <- rp$node[2]+nodebank[node==rp$node[2],pos]+2
  }
  
  ### update tree
  split_tip_edge <- which(tree$edge[,2]==splitting.tip)
  nedges <- nrow(tree$edge)
  new_edge_mat <- matrix(c(rp$node[2]+1,rep(new_node,2),new_node,splitting.tip,N+1),nrow=3,byrow=F)
  
  tree$edge[tree$edge %in% neg_nds] <- tree$edge[tree$edge %in% neg_nds]+1
  tree$edge[tree$edge>N] <- tree$edge[tree$edge>N]+1 ##adding new species ads +1 to all nodes
  
  if (split_tip_edge==1){
    tree$edge <- rbind(new_edge_mat,tree$edge[2:nedges,])
    tree$edge.length <- c(tree$edge.length[1],0,0,tree$edge.length[2:nedges])
  } else if (split_tip_edge==nedges){
    tree$edge <- rbind(tree$edge[1:(split_tip_edge-1),],new_edge_mat)
    tree$edge.length <- c(tree$edge.length,0,0)
  } else {
    tree$edge <- rbind(tree$edge[1:(split_tip_edge-1),],
                       new_edge_mat,
                       tree$edge[(split_tip_edge+1):nedges,])
    tree$edge.length <- c(tree$edge.length[1:(split_tip_edge)],0,0,
                          tree$edge.length[(split_tip_edge+1):nedges])
  }
  new_edges <- 1:2+split_tip_edge
  tree$tip.label <- c(tree$tip.label,sp(N+1))
  tree$Nnode <- tree$Nnode+1
  ### add our new node in the right place: precisely where the splitting tip was found
  
  ### update nodebank
   #account for new node by adding +1 on correct side (pos/neg) for all nodes on rootpath
  if (just_tree){
    return(tree)
  } else {
    for (i in 2:nrow(rp)){
      if (rp$orientation[i]=='pos'){
        nodebank[node==rp$node[i],pos:=pos+1]
      } else {
        nodebank[node==rp$node[i],neg:=neg+1]
      }
    }
    nodebank[node %in% neg_nds,node:=node+1]   #account for new species by adding +1 to all nodes on neg side
    nodebank[,node:=node+1]                 #account for new species by adding +1 to all nodes
    ## add new node
    nodebank <- rbind(nodebank,data.table('node'=new_node,'pos'=0,'neg'=0))
    setkey(nodebank,node,neg)
    return(list('tree'=tree,'nodebank'=nodebank))
  }
}
