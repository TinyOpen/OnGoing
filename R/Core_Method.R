#' Boosa Similarity
#'
#' Calculate the similarity matrix of Bossa scores which are obtained by boosting
#' on single attribute.
#'
#' @param is.pca A logical variable indicating if the Bossa scores should transformed
#' to principle components and then calculate the similarity matrix. It is recommended
#' when processing the ultra-dimention data.
#' @param pca.sum.prop A numeric indicating how many components should be reserved
#' in order to make this propotion of variance. The default is \code{pca.sum.prop =  0.95}.
#' @param n.comp The number of components of PCA. The default is \code{n.comp = 50}.
#' @param fix.pca.comp A numeric variable indicating whether choosing the fixed
#' number of components or the fixed porpotion of variance and the default is to
#' choose fixed porpotion.
#' @param alpha A power scaling for Bossa scores, representing the weight of
#' variable sigma value.
#'
#' @return An object including Bossa scores, original data, disimilarity matrix and
#' similarity matrix.
#'
#' @export
BossaSimi <- function(data, is.pca = TRUE, pca.sum.prop = 0.95, fix.pca.comp = FALSE,
                      n.comp = 50, alpha = 1) {
  if (!is.data.frame(data)) {
    if (is.matrix(data)) {
      data <- as.data.frame(data)
      warning("The input data should be a data frame.", call. = FALSE)
    }
  }
  
  n <- dim(data)[1]
  p <- dim(data)[2]
  
  if (!all(apply(data, 2, sum) != 0))
    warning("The input data contains 'All Zero' gene which should be deleted.")
  
  if (!all(complete.cases(data)))
    warning("The input data contains missing value.")
  data <- na.omit(data)
  
  # if (!all(data == as.integer(data))) warning('the input data should be
  # binary or odinary')
  print("Do transformation...")
  U.score <- apply(data, 2, FUN = function(x) {
    expres <- rle(sort(x))
    expres.level <- expres$values
    prop.level <- expres$lengths/n
    n.level <- length(expres.level)
    P <- prop.level[-n.level]
    q <- qnorm(cumsum(P))
    e <- -diff(c(0, exp(-q^2/2), 0))/sqrt(2 * pi)/prop.level
    U.score.x <- vector(length = n)
    for (j in 1:n.level) {
      U.score.x[x == expres.level[j]] <- e[j]
    }
    
    return(U.score.x)
  })
  if (is.pca) {
    data.pca <- prcomp(U.score)
    data.pca.sdev <- data.pca$sdev
    var.prop <- data.pca.sdev^2/sum(data.pca.sdev^2)
    var.prop.sum <- cumsum(var.prop)
    sum.prop.index <- which(var.prop.sum > pca.sum.prop)[1]
    if (!fix.pca.comp) {
      print(paste("need", sum.prop.index, "components to get", pca.sum.prop,
                  "% variance"))
      data.pca.x <- data.pca$x[, 1:sum.prop.index]
    } else {
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
  return(list(data = data, U.score = U.score, bossa.disimi = bossa.disimi,
              bossa.simi = bossa.simi))
}



#' Overlap Clustering
#'
#' Do overlap clustering with Bossa similarity in different levels of \code{p}.
#'
#' @param data.simi The similarity matrix of Bossa scores.
#'
OverlapClust <- function(data.simi) {
  n <- dim(data.simi)[1]
  overlap.clu <- cbind(first.clu = rep(0, n), belong.layer = rep(0, n))
  n.clu <- 1
  clust <- list()
  
  simi.tri <- data.simi[lower.tri(data.simi)]
  core.thresh <- quantile(simi.tri[simi.tri > 0], p)
  self <- diag(data.simi)
  clust.center <- data.frame()
  
  for (l in 1:length(p)) {
    thresh <- core.thresh[l]
    jump0 <- 0
    
    # Look for the critical centers ------------------------
    
    while (jump0 == 0) {
      max.var.idx <- which.max(self)[1]
      candidt <- sort(data.simi[max.var.idx, ], decreasing = T, index.return = T)
      candidt.idx <- (1:n)[candidt$ix[candidt$x > thresh]]
      
      clust[[n.clu]] <- max.var.idx
      
      for (i in candidt.idx) {
        if (i == max.var.idx)
          next
        candidt.simi <- data.simi[i, c(unlist(clust[[n.clu]]),
                                       i)]
        candidt.simi.check <- candidt.simi[1]  # distance between i and the center
        
        if (lin < 0 || lin >= 1)
          warning("The parameter 'lin' is not appropriate.")
        
        candidt.simi.check <- ifelse(lin > 0 & lin < 1, quantile(candidt.simi,
                                                                 lin), ifelse(lin == 0, min(candidt.simi), candidt.simi.check))
        
        
        if (candidt.simi.check > thresh)
          clust[[n.clu]] <- append(clust[[n.clu]], i)
      }
      
      sum.new.idx <- sum(overlap.clu[unlist(clust[[n.clu]]), 1] ==
                           0)
      print(paste(sum.new.idx, " new points are assigned to the current cluster, sum to ",
                  length(clust[[n.clu]]), " cells.", sep = ""))
      
      if (sum.new.idx > 1 & n.clu <= 50) {
        
        self[unlist(clust[[n.clu]])] <- 0
        
        clust.new <- ifelse(1:n %in% unlist(clust[[n.clu]]), 1,
                            0)
        overlap.clu <- cbind(overlap.clu, clust.new)
        
        change.position.1 <- as.logical(clust.new) & overlap.clu[,
                                                                 1] == 0
        overlap.clu[change.position.1, 1] <- n.clu
        change.position.2 <- as.logical(clust.new) & overlap.clu[,
                                                                 2] == 0
        overlap.clu[change.position.2, 2] <- l
        
        clust.center <- rbind(clust.center, c(l, max.var.idx))
        n.clu <- n.clu + 1
        
      } else {
        jump0 <- 1
        clust <- clust[-n.clu]
        if (n.clu == 1)
          n.clu <- 0
        
      }
    }
    
    # Consider the non-neighbours of max variance points ------------------
    
    if (n.clu == 0)
      next
    sum.layer.0 <- sum(overlap.clu[, 2] == l)
    
    candidt.idx.2 <- which(overlap.clu[, 2] != l)
    candidt.clu <- which(clust.center[, 1] == l)
    
    for (i in candidt.idx.2) {
      clu.simi <- sapply(candidt.clu, FUN = function(j) {
        candidt.simi <- data.simi[i, unlist(clust[[j]])]
        candidt.simi.check <- max(candidt.simi)
        
        candidt.simi.check <- ifelse(lin > 0 & lin < 1, quantile(candidt.simi,
                                                                 lin), ifelse(lin == 0, min(candidt.simi), candidt.simi.check))
        return(candidt.simi.check)
      })
      
      max.simi <- max(clu.simi)
      
      if (max.simi > thresh) {
        max.clu.idx <- candidt.clu[which.max(clu.simi)]
        clust[[max.clu.idx]] <- append(clust[[max.clu.idx]], i)
        if (overlap.clu[i, 1] == 0)
          overlap.clu[i, 1] <- max.clu.idx
        if (overlap.clu[i, 2] == 0)
          overlap.clu[i, 2] <- l
        overlap.clu[i, 2 + max.clu.idx] <- 1
      }
    }
    
    sum.layer.1 <- sum(overlap.clu[, 2] == l)
    
    if (l > 1)
      print(paste("After second search, there are", sum.layer.1 -
                    sum.layer.0, "new cells are assigned.", sep = " "))
  }
  
  return(list(overlap.clu = overlap.clu, clust.center = clust.center))
  
}



#' Clust Share
#'
#' Calculate the share matrix among overlap clusters.
#'
#' @param obs The result of overlap clustering.
#'
clush <- function(obs) {
  m <- dim(obs)[2]
  shmat <- matrix(0, m, m)
  for (i in 1:(m - 1)) {
    for (j in i:m) {
      shmat[i, j] <- shmat[j, i] <- sum(obs[, i] == 1 & obs[, j] ==
                                          1)
    }
  }
  diag(shmat) <- colSums(obs)
  return(shmat)
}



#' Merge test
#'
#' Do a test when try to merge two overlap clusters.
#'
#' @param data A data.frame with two columns with dichotomous or polytomous.
#' @param r A numeric indicating tetrachoric or polychoric correlation coefficient.
#'
testlmc <- function(data, r) {
  n <- dim(data)[1]
  P <- list()
  for (i in 1:2) {
    dat <- data[, i]
    rlei <- rle(sort(dat))
    labeli <- rlei$values
    pri <- rlei$lengths/n
    leni <- length(labeli)
    P[[i]] <- pri[-leni]
  }
  p1 <- P[[1]]
  p2 <- P[[2]]
  K1 <- length(p1) + 1
  K2 <- length(p2) + 1
  q1 <- qnorm(cumsum(p1))
  q2 <- qnorm(cumsum(p2))
  e1 <- -diff(c(0, exp(-q1^2/2), 0))/sqrt(2 * pi)/c(p1, 1 - sum(p1))
  e2 <- -diff(c(0, exp(-q2^2/2), 0))/sqrt(2 * pi)/c(p2, 1 - sum(p2))
  U1 <- e1[data[, 1] + 1]
  U2 <- e2[data[, 2] + 1]
  stat <- sum(U1 * U2)
  Low <- c(-Inf, q1)
  Up <- c(q1, Inf)
  q2l <- c(-Inf, q2)
  q2u <- c(q2, Inf)
  cvs <- 0
  cv <- 0
  
  for (i in 1:K1) {
    
    f1 <- function(x) {
      dnorm(x) * pnorm((q2u[1] - r * x)/sqrt(1 - r^2))
    }
    If1 <- integrate(f1, lower = Low[i], upper = Up[i])$value
    cv <- cv + If1 * e1[i] * e2[1]
    cvs <- cvs + If1 * e1[i]^2 * e2[1]^2
    
    fk2 <- function(x) {
      dnorm(x) * pnorm((r * x - q2l[K2])/sqrt(1 - r^2))
    }
    Ifk2 <- integrate(fk2, lower = Low[i], upper = Up[i])$value
    cv <- cv + Ifk2 * e1[i] * e2[K2]
    cvs <- cvs + Ifk2 * e1[i]^2 * e2[K2]^2
    
    if (K2 > 2) {
      for (j in 2:(K2 - 1)) {
        f <- function(x) {
          dnorm(x) * (pnorm((q2u[j] - r * x)/sqrt(1 - r^2)) - pnorm((q2l[j] -
                                                                       r * x)/sqrt(1 - r^2)))
        }
        If <- integrate(f, lower = Low[i], upper = Up[i])$value
        cv <- cv + If * e1[i] * e2[j]
        cvs <- cvs + If * e1[i]^2 * e2[j]^2
      }
    }
  }
  
  vu <- n * (cvs - cv^2)
  sstat <- stat/sqrt(vu)
  pval <- 2 * (1 - pnorm(abs(sstat)))
  return(list(stat = sstat, pval = pval))
}



#' Key feature
#'
#' Determine key features of each overlap clusters.
#'
#' @param cdata Overlap clustering result.
#' @param criter A numeric indicating tetrachoric or polychoric correlation coefficient.
#'
keyfeat <- function(cdata, criter) {
  n <- dim(cdata)[1]
  m <- dim(cdata)[2]
  cmax <- max(apply(cdata, 2, max))
  if (cmax == 1) {
    r <- tetrachoric(cdata)$rho
  } else {
    r <- polychoric(cdata)$rho
  }
  N <- m * (m - 1)/2
  Stats <- matrix(0, m, m)
  Pval <- matrix(1, m, m)
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      # coding change to 0:k-------------
      c1 <- cdata[, i]
      uc1 <- unique(c1)
      sc1 <- sort(uc1)
      for (c1i in sc1) {
        c1[c1 == c1i] <- (1:length(uc1))[sc1 == c1i] - 1
      }
      c2 <- cdata[, j]
      uc2 <- unique(c2)
      sc2 <- sort(uc2)
      for (c2i in sc2) {
        c2[c2 == c2i] <- (1:length(uc2))[sc2 == c2i] - 1
      }
      tij <- testlmc(cbind(c1, c2), r[i, j])
      Stats[i, j] <- tij$stat
      Pval[i, j] <- tij$pval
    }
  }
  Sp <- Stats
  Sp[Pval > criter] <- 0
  Pp <- Pval
  Pp[Pval > criter] <- 1
  Pp[Sp < 0] <- 1
  psortp <- sort(Pp, index.return = T)
  featindp <- psortp$ix[psortp$x < 1]
  flenp <- length(featindp)
  if (flenp == 0) {
    kfeatp <- 0
  } else {
    kfeatp <- rep(0, flenp)
    for (ii in 1:flenp) {
      fii <- floor(featindp[ii]/m)
      fij <- featindp[ii] - fii * m
      kfeatp[ii] <- paste(colnames(cdata)[fij], colnames(cdata)[fii +
                                                                  1], sep = ",")
    }
  }
  Sn <- Stats
  Sn[Pval > criter] <- 0
  Pn <- Pval
  Pn[Pval > criter] <- 1
  Pn[Sn > 0] <- 1
  psortn <- sort(Pn, index.return = T)
  featindn <- psortn$ix[psortn$x < 1]
  flenn <- length(featindn)
  if (flenn == 0) {
    kfeatn <- 0
  } else {
    kfeatn <- rep(0, flenn)
    for (ii in 1:flenn) {
      fii <- floor(featindn[ii]/m)
      fij <- featindn[ii] - fii * m
      kfeatn[ii] <- paste(colnames(cdata)[fij], colnames(cdata)[fii +
                                                                  1], sep = ",")
    }
  }
  stat <- Stats
  stat <- stat + t(stat) - diag(diag(stat))
  maxs <- max(stat)
  mins <- min(stat)
  crit <- qnorm(criter)
  stat <- (maxs - stat)/(maxs - mins)
  scrit0 <- (maxs + crit)/(maxs - mins)
  scrit1 <- (maxs - crit)/(maxs - mins)
  return(list(kfp = kfeatp, kfn = kfeatn, stat = stat, scrit0 = scrit0,
              scrit1 = scrit1))
}


#' Clusr Merge
#'
#' With every recommended k, merge the overlap clusters depending on the distance
#' between those.
#'
#' @param tree.size The label of each overlap clusters.
#'
ClustMerge <- function(tree.size) {
  k <- max(tree.size)
  tree.merge.index <- cutree(clu.hc, k = k)
  new.overlap.clu <- data.frame(cell = c(1:n))
  
  for (clu in 1:k) {
    clu.ind <- which(tree.merge.index == clu)
    new.overlap <- overlap.clu[, (clu.ind[1] + 2)]
    
    for (j in clu.ind) {
      new.overlap <- (overlap.clu[, j + 2] | new.overlap)
    }
    
    new.overlap.clu <- cbind(new.overlap.clu, new.overlap)
  }
  colnames(new.overlap.clu) <- c("cell", paste("merge.clu", 1:k, sep = ""))
  for (i in 1:dim(new.overlap.clu)[1]) {
    for (j in 1:dim(new.overlap.clu)[2]) {
      new.overlap.clu[i, j] <- as.numeric(new.overlap.clu[i, j])
    }
  }
  return(new.overlap.clu)
}
