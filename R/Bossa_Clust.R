#' Bossa Clustering
#'
#' With the previous calculated similarity matrix or the original categorical
#' dataframe, the results of both overlap clustering and hierarchical clustering
#' are obtained with several recommended cluster numbers(k) aftering processing the merge
#' cluster step.
#'
#' @param data an original categorical data with n observations and p variables.
#' @param data.pre an list obtained by \code{\link{BoosaSimi}} including original
#' categorical data, similarity matrix, disimilarity matrix and transformed data,
#' Bossa scores. It is recommended to calculate the data.pre first and then do
#' \code{\link{BoosaClust}} in order to save time when trying to change parameters
#' of this function.
#' @param alpha A power scaling for Bossa scores, representing the weight of
#' variable sigma value.
#' @param p A set of quantiles(90%, 75% and median) of the positive values of
#' similarity matrix to form clusters at different levels of within-cluster similarity.
#' @param lin A tuning parameter to control the size of each overlap cluster before
#' merging, smaller lin leads to larger cluster size.
#' @param is.pca A logical variable indicating if the Bossa scores should transformed
#' to principle components and then calculate the similarity matrix. It is recommended
#' when processing the ultra-dimention data.
#' @param pca.sum.prop A numeric indicating how many components should be reserved
#' in order to make this propotion of variance. The default is \code{pca.sum.prop =  0.95}.
#' @param n.comp The number of components of PCA. The default is \code{n.comp = 50}.
#' @param fix.pca.comp A numeric variable indicating whether choosing the fixed
#' number of components or the fixed porpotion of variance and the default is to
#' choose fixed porpotion.
#' @param cri A tuning parameter, if p value smaller than cri, then reject
#' the NULL hypothesis and merge overlap subclusters. And cri can be any numeric less
#' than \code{1}, if \code{cri = 1} then the criteria will be reset to \code{0.05/N}
#' (N is the numer of all overlap subcluster), and if \code{cri = 2} then the
#' criteria \code{0.05/N(N-1)}.
#' @param lintype The agglomeration method to be used in \code{\link[stats]{hclust}}.
#' This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2",
#' "single", "complete", "average" and so on. The default is "ward.D2".
#'
#' @import  Rtsne
#'
#' @return An object including overlap clusters after merging and non-overlap 
#' clusters, which can be showed by function \code{\link{bossa_interactive}}
#' @export


BossaClust <- function(data, data.pre = NULL, alpha = 1, p = c(0.9, 0.75, 0.5),
                       lin = 0.25, is.pca = TRUE, pca.sum.prop = 0.95, n.comp = 50,
                       fix.pca.comp = FALSE, cri = 1, lintype = "ward.D2", perplexity = 30)
{
  
  # Check input data --------------------------
  try(if (is.null(data.pre) & is.null(data))
    stop("No data to process"))

  
  if (is.null(data.pre)) {
    data.pre <- BossaSimi(data, is.pca = is.pca, pca.sum.prop = pca.sum.prop,
                          fix.pca.comp = fix.pca.comp, n.comp = n.comp, alpha = alpha)
  }
  
  data.simi <- data.pre$bossa.simi
  
  n <- dim(data.simi)[1]
  cell.name <- row.names(data.simi)
  
  try(if (dim(data.simi)[1] != dim(data)[1])
    stop("there is conflict between data.pre and data."))
  
  # Do overlap cluster with 'SC' method -----------------------------
  cat("Do overlap cluster...\n")
  overlap.pre <- OverlapClust(data.simi)
  overlap.clu <- overlap.pre$overlap.clu
  clust.center <- overlap.pre$clust.center
  
  # Merge clusters-----------------------------
  cat("\nMerge some subclusters...\n")
  sum.clu <- dim(overlap.clu)[2] - 2
  
  # Cluster of the left cells--------------
  overlap.clu <- LeftClust(overlap.clu, sum.clu, data.simi, n = n)

  # Set some output variables

  if (sum.clu < 2){
    cat("/n there just one subcluster.")
    clu.merge <- overlap.clu[,-c(1,2)]
    clu.merge <- data.frame(cell = rownames(data.simi), clust = overlap.clu[,3])
    clu.merge <- as.list(clu.merge)
    
    tree.max = 1
    tree.min = 1
    
    scrit0 = NULL
    scrit1 = NULL
    share.mat = NULL
    clu.match = NULL
  }
  if (sum.clu >= 2){
    share.mat <- ClustShare(overlap.clu[, -c(1,2)]) # calculate the overlap of each clusters
    clu.match <- KeyFeature(overlap.clu, cri, sum.clu = sum.clu)  # determine key features of each overlap clusters
    scrit0 <- clu.match$scrit0
    scrit1 <- clu.match$scrit1
    mer.clu <- clu.match$mer.clu
    clu.dis <- as.dist(clu.match$clu.math.stat)
    
    clu.hc <- hclust(clu.dis, lintype) # calculate the recommended number of clusters
    tree.max <- max(cutree(clu.hc, h = scrit0))
    tree.min <- max(cutree(clu.hc, h = scrit1))
    
    # Merge overlap clusters with recommended k 
    clu.merge <- lapply(tree.min:tree.max, function(x){ # merge the overlap subclusters
      ClustMerge(x, clu.hc = clu.hc, n = n, overlap.clu = overlap.clu, data.simi = data.simi)}
    )
  }

  # Do non-overlap clustering with recommended k and do test find different variables(genes)
  cell.hc.clust <- sapply(tree.min:tree.max, function(x) {
    hc.clust <- hclust(as.dist(data.pre$bossa.disimi), lintype)
    order.clust <- OrderClust(hc.clust, k = x)
    order.clust
  }) # it's a data.frame
  

  # HC: find different variables(genes) and prepare the data for heatmap
  U.score.non.pca <- data.pre$U.score.non.pca 
  hc.de.plot.all <- FindHcDe(cell.hc.clust = cell.hc.clust, 
                         U.score.non.pca = U.score.non.pca, cell.name = cell.name)
  hc.de.plot <- hc.de.plot.all$hc.de.plot
  # Overlap: find different variables(genes) and prepare the data for heatmap
  
  overlap.de.plot.all <- FindOverlapDe(clu.merge, U.score.non.pca, cell.name = cell.name, n = n)
  overlap.de.plot <- overlap.de.plot.all$overlap.de.plot
  clu.merge <- overlap.de.plot.all$clu.merge.1
  
  clu.merge.share <- lapply(clu.merge, FUN = function(x){
    ClustShare(x[,-1])
  })
  
  # Before Merge: find different variables(genes) and prepare the data for heatmap
  bef.de.plot.all <- FindBefDe(mer.clu, overlap.clu, U.score.non.pca, cell.name)
  bef.de.plot <- bef.de.plot.all$bef.de.plot
  mer.clu <- bef.de.plot.all$mer.clu
  
  # Do tsne for visualization--------------------------
  cat("\n\nDo tsne....\n")
  my.tsne <- Rtsne(data, perplexity = perplexity)
  cell = matrix(1:n, n)
  tsne.y <- cbind(my.tsne$Y, cell) 
  cat("tsne done.\n")
  
  K.level <- tree.max - tree.min + 1
  
  # prepare the data for Overlap tsne visualization 
  overlap.melt.data <- lapply(1:K.level, function(x){
    OverlapMelt(K = x, overlap.clu = clu.merge, data.tsne = tsne.y)
  })
  
  overlap.note <- OverlapNote(clu.merge)

  return(list(overlap.clu = clu.merge, non.overlap.clu = cell.hc.clust,
              ori.overlap = overlap.clu, clu.dis = clu.dis, clu.share = share.mat,
              clu.merge.share = clu.merge.share, mer.clu = mer.clu, 
              tree.max = tree.max, tree.min = tree.min, tsne.y = tsne.y,
              hc.de.plot = hc.de.plot, overlap.de.plot = overlap.de.plot, bef.de.plot = bef.de.plot,
              overlap.melt.data = overlap.melt.data, overlap.clu.res = overlap.note))
  
}

