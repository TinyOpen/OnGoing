#' Bind two factors
#'
#' Create a new factor from two existing factors, where the new factor's levels
#' are the union of the levels of the input factors.
#'
#' @param a factor
#' @param b factor
#'
#' @return factor
#' @export
#' @examples
#'#' fbind(iris$Species[c(1, 51, 101)], PlantGrowth$group[c(1, 11, 21)])

BossaClust <- function(data = NULL, data.pre = NULL, alpha = 1, 
                       p = c(0.9, 0.8, 0.7), lin = 0.2,
                       is.pca = TRUE, pca.sum.prop = 0.95, fix.pca.comp = FALSE, n.comp = 50,
                       cri = 1, lintype = "ward.D2", perplexity = 30)
{
  
  require("psych")
  
  # Check input data --------------------------
  try(if(is.null(data.pre) & is.null(data)) stop("No data to process"))
  
  if(is.null(data.pre)){
    print("Do boosa transformation and calculate the similarity matrix and disimilarity matrix...")
    data.pre <- BossaSimi(data, is.pca = is.pca, pca.sum.prop = pca.sum.prop, 
                          fix.pca.comp = fix.pca.comp, n.comp = n.comp, alpha = alpha)
    n <- dim(data)[1]
  }
  
  data.simi <- data.pre$bossa.simi
  data.dis <- data.pre$bossa.disimi
  
  if(is.null(data)){
    n <- dim(data.pre$data)[1]
    data <- data.pre$data
  } 
  
  # Do overlap cluster with "SC" method -----------------------------
  print("Do overlap cluster...")
  overlap.pre <- OverlapClust(data.simi, p = p, lin = lin)
  overlap.clu <- overlap.pre$overlap.clu
  clust.center <- overlap.pre$clust.center
  
  # Merge clusters-----------------------------
  print("Merge some subclusters...")
  sum.clu <- dim(overlap.clu)[2] - 2
  colnames(overlap.clu) <- c("first.clu", "belong.layer", 
                             paste("clust.", 1:sum.clu, sep = ""))
  ori.clu <- overlap.clu[,-c(1,2)]
  
  shmat <- clush(overlap.clu[, -c(1, 2)])
  sig.lev <- ifelse(cri < 1, cri,
                    ifelse(cri == 1, 0.05/sum.clu, 0.05/sum.clu/(sum.clu-1)))
  
  
  if(sum.clu < 2) return(list(clust.center = clust.center, 
                              overlap.clu = overlap.clu, shmat = shmat, p = p))
  
  clumatch<-keyfeat(ori.clu, sig.lev)
  scrit0<-clumatch$scrit0
  scrit1<-clumatch$scrit1
  
  clu.dis<-as.dist(clumatch$stat)
  merclu<-clumatch$kfp
  sepclu<-clumatch$kfn
  
  # Take charge of the left cells
  non.core.ind <- (1:n)[apply(overlap.clu[, -c(1,2)], 1, sum) == 0]
  k1<-dim(clust.center)[1]
  for (i in non.core.ind){
    maxci<-rep(0,k1)
    ij<-0
    for (j in 1:k1){			
      ij<-ij+1	
      maxci[ij]<-quantile(data.simi[i,overlap.clu[,1]==j], 0.5)
    }
    max.ind<-which.max(maxci)
    
    overlap.clu[i,1]<-max.ind
    overlap.clu[i,(max.ind+2)]<-1
  }
  
  clu.hc <- hclust(clu.dis,lintype)
  tree.max <-max(cutree(clu.hc, h = scrit0))
  tree.min <-max(cutree(clu.hc, h = scrit1))
  
  clu.merge <- sapply(tree.min:tree.max, ClustMerge)  
  cell.hc.clust <- sapply(tree.min:tree.max, function(x){
    hc.clust <- hclust(as.dist(data.pre$bossa.disimi), lintype)
    hc.tree <- cutree(hc.clust, k = x)
    hc.tree
  })
  
  # Do tsne for visualization--------------------------
  print("Do tsne....")
  my.tsne <- Rtsne(data, perplexity = perplexity)
  tsne.y <- transform(my.tsne$Y, cell = 1:301)
  print("tsne done.")
  
  
  return(list(overlap.clu = clu.merge, non.overlap.clu = cell.hc.clust,
              ori.overlap = overlap.clu, clust.center = clust.center,
              clu.dis = clu.dis, tree.max = tree.max, tree.min = tree.min,
              cell.simi = data.simi, tsne.y = tsne.y, data.pre = data.pre))
  
}