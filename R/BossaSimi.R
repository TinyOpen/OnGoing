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


BossaSimi <- function(data, is.pca = TRUE, pca.sum.prop = 0.95, fix.pca.comp = FALSE, 
                      n.comp = 50, alpha = 1)
{
  if (!is.data.frame(data))
  {
    if (is.matrix(data))
    {
      data <- as.data.frame(data)
      warning("The input data should be a data frame.", call. = FALSE)
    }
  }
  
  n <- dim(data)[1]
  p <- dim(data)[2]
  
  if (!all(apply(data, 2, sum) != 0)) warning("The input data contains 'All Zero' gene which should be deleted.")
  
  if (!all(complete.cases(data))) warning("The input data contains missing value.")
  data <- na.omit(data)
  
  #if (!all(data == as.integer(data)))
  #warning("the input data should be binary or odinary")
  
  U.score <- apply(data, 2, FUN = function(x)
  {
    expres <- rle(sort(x))
    expres.level <- expres$values
    prop.level <- expres$lengths/n
    n.level <- length(expres.level)
    P <- prop.level[-n.level]
    q <- qnorm(cumsum(P))
    e <- -diff(c(0, exp(-q^2/2), 0))/sqrt(2 * pi)/prop.level
    U.score.x <- vector(length = n)
    for (j in 1:n.level)
    {
      U.score.x[x == expres.level[j]] <- e[j]
    }
    return(U.score.x)
  })
  if (is.pca)
  {
    data.pca <- prcomp(U.score)
    data.pca.sdev <- data.pca$sdev
    var.prop <- data.pca.sdev^2/sum(data.pca.sdev^2)
    var.prop.sum <- cumsum(var.prop)
    sum.prop.index <- which(var.prop.sum > pca.sum.prop)[1]
    if (!fix.pca.comp)
    {
      print(paste("need", sum.prop.index, "components to get", pca.sum.prop, 
                  "% variance"))
      data.pca.x <- data.pca$x[, 1:sum.prop.index]
    }
    else
    {
      if (n.comp > n | n.comp > p) 
        warning("The proposed number of components is more than the rank of data.", 
                call. = FALSE)
      data.pca.x <- data.pca$x[, 1:min(n.comp, n, p)]
    }
    U.score <- data.pca.x
  }
  U.score.scale <- apply(U.score, 2, FUN = function(x) x * var(x)^(alpha/2))
  bossa.simi <- U.score.scale %*% cor(U.score) %*% t(U.score.scale)
  max <- max(bossa.simi)
  min <- min(bossa.simi)
  diff <- max - min
  bossa.disimi <- (max - bossa.simi)/(max - min)
  return(list(bossa.disimi = bossa.disimi, U.score = U.score, bossa.simi = bossa.simi))
}

clush <- function(obs){
  m<-dim(obs)[2]
  shmat<-matrix(0,m,m)
  for(i in 1:(m-1)){
    for(j in i:m){
      shmat[i,j]<-shmat[j,i]<-sum(obs[,i]==1 & obs[,j]==1)
    }
  }
  diag(shmat)<-colSums(obs)
  return(shmat)
}

BossaClust <- function(data.simi, data = NULL, data.dis = NULL, method = "SC", alpha = 1, 
                       p = c(0.9,  0.8, 0.7,0.5), lin = 0.2,
                       is.pca = TRUE, pca.sum.prop = 0.95, fix.pca.comp = FALSE, n.comp = 50)
{
  
  
  require("psych")
  stopifnot(method %in% c("SC", "HC"))
  if (!is.null(data.dis)) {
    stopifnot(is.dist(data.dis)) 
    stopifnot(dim(data.dis)[1] != dim(data.dis)[2])
  }
  
  # Check input data for "SC" method --------------------------
  
  if (method == "SC")
  {
    if (is.null(data.simi))
    {
      data.simi <- BossaSimi(data, is.pca = is.pca, pca.sum.prop = pca.sum.prop, 
                             fix.pca.comp = fix.pca.comp, n.comp = n.comp, alpha = alpha)$bossa.simi
      n <- dim(data)[1]
    } else {
      stopifnot(all(apply(data.simi, 2, is.numeric)))
      stopifnot(dim(data.simi)[1] == dim(data.simi)[2])
      
      if (!is.null(data)) 
      {
        warning("Original data is useless since 'data.simi' isn't NULL")
        if((dim(data)[1] != dim(data.simi)[1])) warning("The row of original binary data should equal to that of data.simi.")
      }
      n <- dim(data.simi)[1]
    }
  }
  
  # Do overlap cluster with "SC" method -----------------------------
  require("psych")
  overlap.clu <- cbind(first.clu = rep(0, n), belong.layer = rep(0, n))
  n.clu <- 1
  clust <- list()
  
  
  simi.tri <- data.simi[lower.tri(data.simi)]
  core.thresh <- quantile(simi.tri[simi.tri > 0], p)
  self <- diag(data.simi)
  clust.center <- data.frame()
  
  for (l in 1:length(p))
  {
    thresh <- core.thresh[l]
    jump0 <- 0
    
    # Look for the critical centers ------------------------
    
    while (jump0 == 0)
    {
      max.var.idx <- which.max(self)[1]
      candidt <- sort(data.simi[max.var.idx, ], decreasing = T, index.return = T)
      candidt.idx <- (1:n)[candidt$ix[candidt$x > thresh]]
      
      clust[[n.clu]] <- max.var.idx
      
      for(i in candidt.idx){
        if(i == max.var.idx) next
        candidt.simi <- data.simi[i, c(unlist(clust[[n.clu]]), i)]
        candidt.simi.check <- candidt.simi[1] # distance between i and the center
        
        if (lin < 0 || lin >= 1) warning("The parameter 'lin' is not appropriate.")
        
        candidt.simi.check <- ifelse(lin > 0 & lin < 1, quantile(candidt.simi, lin),
                                     ifelse(lin == 0, min(candidt.simi), candidt.simi.check))
        
        
        if(candidt.simi.check > thresh) clust[[n.clu]] <- append(clust[[n.clu]], i)
      }
      
      sum.new.idx <- sum(overlap.clu[unlist(clust[[n.clu]]), 1] == 0)
      print(paste(sum.new.idx, " new points are assigned to the current cluster, sum to ",
                  length(clust[[n.clu]]), " cells.", sep = ""))
      
      if(sum.new.idx > 1 & n.clu <= 50){
        
        self[unlist(clust[[n.clu]])] <- 0
        
        clust.new <- ifelse(1:n %in% unlist(clust[[n.clu]]), 1, 0)
        overlap.clu <- cbind(overlap.clu, clust.new)
        
        change.position.1 <- as.logical(clust.new) & overlap.clu[, 1] == 0
        overlap.clu[change.position.1, 1] <- n.clu
        change.position.2 <- as.logical(clust.new) & overlap.clu[, 2] == 0
        overlap.clu[change.position.2, 2] <- l
        
        clust.center <- rbind(clust.center, c(l, max.var.idx))
        n.clu <- n.clu + 1
        
      } else {
        jump0 <- 1
        clust <- clust[-n.clu]
        if(n.clu == 1) n.clu <- 0
        
      }
    }
    
    # Consider the non-neighbours of max variance points ------------------
    
    if(n.clu == 0) next
    sum.layer.0 <- sum(overlap.clu[, 2] == l)
    
    candidt.idx.2 <- which(overlap.clu[, 2] != l)
    candidt.clu <- which(clust.center[,1] == l)
    
    for(i in candidt.idx.2)
    {
      clu.simi <- sapply(candidt.clu, FUN = function(j){
        candidt.simi <- data.simi[i, unlist(clust[[j]])]
        candidt.simi.check <- max(candidt.simi)
        
        candidt.simi.check <- ifelse(lin > 0 & lin < 1, quantile(candidt.simi, lin),
                                     ifelse(lin == 0, min(candidt.simi), candidt.simi.check))
        return(candidt.simi.check)
      })
      
      max.simi <- max(clu.simi)
      
      if(max.simi > thresh) {
        max.clu.idx <- candidt.clu[which.max(clu.simi)]
        clust[[max.clu.idx]] <- append(clust[[max.clu.idx]], i)
        if(overlap.clu[i, 1] == 0) overlap.clu[i, 1] <- max.clu.idx
        if(overlap.clu[i, 2] == 0) overlap.clu[i, 2] <- l
        overlap.clu[i, 2 + max.clu.idx] <- 1
      }
    }
    
    sum.layer.1 <- sum(overlap.clu[, 2] == l)
    
    if(l > 1) print(paste("After second search, there are", sum.layer.1 - sum.layer.0, 
                          "new cells are assigned.", sep = " "))
  }
  
  shmat <- clush(overlap.clu[, -c(1, 2)])
  
  return(list(clust.center = clust.center, overlap.clu = overlap.clu, shmat = shmat))
  
  
  
  if (method == "HC")
  {
    if (is.null(data.dis))
    {
      if (is.null(data.simi)) 
        data.disimi <- BossaSimi(data, is.pca = is.pca, pca.sum.prop = pca.sum.prop, 
                                 fix.pca.comp = fix.pca.comp, n.comp = n.comp, alpha = alpha)$bossa.disimi
      else
      {
        max <- max(data.simi)
        min <- min(data.simi)
        diff <- max - min
        data.disimi <- (max - data.simi)/(max - min)
      }
    }
  }
  
}

