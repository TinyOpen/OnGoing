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
#' @export

BossaClust <- function(data = NULL, data.pre = NULL, alpha = 1, p = c(0.9, 0.75, 0.5),
                       lin = 0.25, is.pca = TRUE, pca.sum.prop = 0.95, n.comp = 50,
                       fix.pca.comp = FALSE, cri = 1, lintype = c, perplexity = 30)
{
  
  # Check input data --------------------------
  try(if (is.null(data.pre) & is.null(data))
    stop("No data to process"))
  
  if (is.null(data.pre)) {
    print("Do boosa transformation and calculate the similarity matrix and disimilarity matrix...")
    data.pre <- BossaSimi(data, is.pca = is.pca, pca.sum.prop = pca.sum.prop,
                          fix.pca.comp = fix.pca.comp, n.comp = n.comp, alpha = alpha)
    n <- dim(data)[1]
  }
  
  data.simi <- data.pre$bossa.simi
  data.dis <- data.pre$bossa.disimi
  
  if (is.null(data)) {
    n <- dim(data.pre$data)[1]
    data <- data.pre$data
  }
  
  # Do overlap cluster with 'SC' method -----------------------------
  print("Do overlap cluster...")
  overlap.pre <- OverlapClust(data.simi, p = p, lin = lin)
  overlap.clu <- overlap.pre$overlap.clu
  clust.center <- overlap.pre$clust.center
  
  # Merge clusters-----------------------------
  print("Merge some subclusters...")
  sum.clu <- dim(overlap.clu)[2] - 2
  colnames(overlap.clu) <- c("first.clu", "belong.layer", paste("clust.",
                                                                1:sum.clu, sep = ""))
  ori.clu <- overlap.clu[, -c(1, 2)]
  
  shmat <- clush(overlap.clu[, -c(1, 2)])
  sig.lev <- ifelse(cri < 1, cri, ifelse(cri == 1, 0.05/sum.clu, 0.05/sum.clu/(sum.clu -
                                                                                 1)))
  
  
  if (sum.clu < 2)
    return(list(clust.center = clust.center, overlap.clu = overlap.clu,
                shmat = shmat, p = p))
  
  clumatch <- keyfeat(ori.clu, sig.lev)
  scrit0 <- clumatch$scrit0
  scrit1 <- clumatch$scrit1
  
  clu.dis <- as.dist(clumatch$stat)
  merclu <- clumatch$kfp
  sepclu <- clumatch$kfn
  
  # Take charge of the left cells--------------
  non.core.ind <- (1:n)[apply(overlap.clu[, -c(1, 2)], 1, sum) == 0]
  k1 <- dim(clust.center)[1]
  for (i in non.core.ind) {
    maxci <- rep(0, k1)
    ij <- 0
    for (j in 1:k1) {
      ij <- ij + 1
      maxci[ij] <- quantile(data.simi[i, overlap.clu[, 1] == j],
                            0.5)
    }
    max.ind <- which.max(maxci)
    
    overlap.clu[i, 1] <- max.ind
    overlap.clu[i, (max.ind + 2)] <- 1
  }
  
  clu.hc <- hclust(clu.dis, lintype)
  tree.max <- max(cutree(clu.hc, h = scrit0))
  tree.min <- max(cutree(clu.hc, h = scrit1))
  
  clu.merge <- sapply(tree.min:tree.max, ClustMerge)
  cell.hc.clust <- sapply(tree.min:tree.max, function(x) {
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
              ori.overlap = overlap.clu, clust.center = clust.center, clu.dis = clu.dis,
              tree.max = tree.max, tree.min = tree.min, cell.simi = data.simi,
              tsne.y = tsne.y, data.pre = data.pre))
  
}
