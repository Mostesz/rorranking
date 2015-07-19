generatePerformances <- function(crits.nr, alts.nr, distribution) {
  if (crits.nr < 2) {
    stop('Criteria number must be greater or equal 2')
  }
  if (!(distribution %in% c('UNIFORM', 'SKEW_NORMAL'))) {
    stop('Incorrect distribution type')
  }
  
  if (distribution == 'UNIFORM') {
    genVector <- genRandomVectorNormalDistribution
  } else if (distribution == 'SKEW_NORMAL') {
    genVector <- genRandomVectorSkewNormalDistribution
  }
  perfs <- matrix(genVector(crits.nr * alts.nr), ncol=crits.nr, nrow=alts.nr)
  while(TRUE) {
    had.domination <- FALSE
    for(i in 1:(alts.nr-1)) {
      for(j in (i+1):alts.nr) {
        while(isDominated(perfs, i, j) || isDominated(perfs, j, i)) {
          had.domination <- TRUE
          perfs[i,] <- genVector(crits.nr)
          perfs[j,] <- genVector(crits.nr)
        }
      }
    }
    has.duplications = FALSE
    for(i in 1:crits.nr) {
      if (length(unique(perfs[,i])) != alts.nr) {
        has.duplications = TRUE
        break
      }
    }
    if (!had.domination && !has.duplications) {
      break
    }
  }
  
  return(perfs)
}

genRandomVectorNormalDistribution <- function(n) {
  return(runif(n, 0, 1))
}

genRandomVectorSkewNormalDistribution <- function(n) {
  return(rsn(n, 5, 2, 5))
}

generatePreferences <- function(perfs, preferences.nr) {
  result <- matrix(0, ncol=2, nrow=preferences.nr)
  alts.nr <- nrow(perfs)  
  edges.max <- alts.nr * (alts.nr - 1) / 2
  triangular.mat.indexes <- sample(edges.max, preferences.nr)
  alts.indexes <- sample(alts.nr)
  for (pref.idx in 1:preferences.nr) {
    row.idx <- getTriangularMatrixRowIndex(triangular.mat.indexes[pref.idx], alts.nr)
    col.idx <- getTriangularMatrixColumnIndex(triangular.mat.indexes[pref.idx], alts.nr)
    result[pref.idx, 1] <- alts.indexes[row.idx]
    result[pref.idx, 2] <- alts.indexes[col.idx]
  }
  return(result)
}

generateRankingPreferences <- function(perfs) {
  result <- t(combn(sample(nrow(perfs)), 2))
  return(result)
}

generatePreferencesFromLinearFunc <- function(perfs, preferences.nr) {
  ranking <- generateRankingPreferencesFromLinearFunc(perfs)
  return(ranking[sample(nrow(ranking), size=preferences.nr),])
}

generateRankingPreferencesFromLinearFunc <- function(perfs) {
  crits.nr <- ncol(perfs)
  
  max.utilities <- rand.fixed.vector(crits.nr, 1, 0.05, 1)
  max.grades <- apply(perfs, 2, max)
  slopes <- max.utilities / max.grades
  
  perfs.utilities <- t(t(perfs) * slopes)
  ranking <- t(combn(sort(rowSums(perfs.utilities), index.return=TRUE, decreasing=TRUE)$ix, 2))
  
  return(ranking)
}

isDominated <- function(perfs, dominant.idx, dominated.idx) {
  for(i in 1:ncol(perfs)) {
    if (perfs[dominated.idx, i] > perfs[dominant.idx, i]) {
      return(FALSE)
    }
  }
  return(TRUE)
}

getTriangularMatrixRowIndex <- function(value.idx, matrix.size) {
  i <- (matrix.size-1) * (matrix.size) / 2 - value.idx
  K <- floor((sqrt(8 * i+1)-1)/2)
  return(matrix.size - 1 - K)
}

getTriangularMatrixColumnIndexUsingRowIndex <- function(value.idx, matrix.size, row.idx) {
  return (value.idx - (matrix.size - 1) * (row.idx - 1) + (row.idx - 1) * row.idx / 2 + 1);
}

getTriangularMatrixColumnIndex <- function(value.idx, matrix.size) {
  row <- getTriangularMatrixRowIndex(value.idx, matrix.size)
  return(getTriangularMatrixColumnIndexUsingRowIndex(value.idx, matrix.size, row))
}

rand.fixed.vector <- function(n, s, a, b) {
  return(rand.fixed.sum.matrix(n, 1, s, a, b)[,1])
}

rand.fixed.sum.matrix <- function(n, m, s, a, b) {
  # Roger Stafford - Jan. 19, 2006
  # http://www.mathworks.com/matlabcentral/fileexchange/9700-random-vectors-with-fixed-sum/content/randfixedsum.m
  
  if ((m != round(m)) || (n != round(n)) || (m < 0) || (n < 1)) {
    stop('n must be a whole number and m a non-negative integer.')
  } else if ((s < n * a) || (s > n * b) || (a >= b)) {
    stop('Inequalities n*a <<- s <<- n*b and a < b must hold.')
  }
  
  s <- (s - n * a) / (b - a)
  
  k <- max(min(floor(s), n - 1), 0) # Must have 0 <<- k <<- n-1
  s <- max(min(s, k+1), k) # Must have k <<- s <<- k+1
  s1 <- s - (k:(k-n+1)) # s1 & s2 will never be negative
  s2 <- (k+n):(k+1) - s
  w <- matrix(0,n,n+1)
  w[1,2] <- .Machine$double.xmax # Scale for full 'double' range
  t <- matrix(0,n-1,n)
  tiny <- .Machine$double.xmin # The smallest positive matlab 'double' no.
  for (i in 2:n) {
    tmp1 <- w[i-1, 2:(i+1)] * s1[1:i] / i
    tmp2 <- w[i-1, 1:i] * s2[(n-i+1):n] / i
    w[i, 2:(i+1)] <- tmp1 + tmp2
    tmp3 <- w[i, 2:(i+1)] + tiny # In case tmp1 & tmp2 are both 0,
    tmp4 <- (s2[(n-i+1):n] > s1[1:i]) # then t is 0 on left & 1 on right
    t[i-1,1:i] <- (tmp2 / tmp3) * tmp4 + (1 - tmp1 / tmp3) * (!tmp4)
  }
  
  v <- n^(3/2) * (w[n,k+2] / .Machine$double.xmax) * (b - a)^(n - 1)
  
  x <- matrix(0,n,m)
  if (m == 0) {
    return(c())
  } # If m is zero, quit with x <- []
  rt <- matrix(runif((n-1) * m), n-1, m) # For random selection of simplex type
  rs <- matrix(runif((n-1) * m), n-1, m) # For random location within a simplex
  s <- matrix(s, 1, m)
  j <- matrix(k+1, 1, m) # For indexing in the t table
  sm <- matrix(0, 1, m)
  pr <- matrix(1, 1, m) # Start with sum zero & product 1
  for (i in (n-1):1) { # Work backwards in the t table ) {
    e <- (rt[n-i,] <= t[i,j]) # Use rt to choose a transition
    sx <- rs[n-i,]^(1 / i) # Use rs to compute next simplex coord.
    sm <- sm + (1 - sx) * pr * s / (i + 1) # Update sum
    pr <- sx * pr # Update product
    x[n-i,] <- sm + pr * e # Calculate x using simplex coords.
    s <- s - e
    j <- j - e # Transition adjustment
  }
  x[n,] <- sm + pr * s # Compute the last x
  
  rp <- matrix(runif(n*m), n, m)# Use rp to carry out a matrix 'randperm'
  p <- matrix(0, n, m)
  for (i in 1:ncol(rp)) {
    p[,i] <- sort(rp[,i], index.return=TRUE)$ix
  }
  x.idxs = p + repmat(matrix(seq(0, n*(m-1), n), 1), n, 1)
  for (i in 1:ncol(x)) {
    x.idxs[,i] <- x.idxs[,i] - (i-1)*nrow(x)
    x[,i] <- x[,i]
  }
  x <- (b - a) * x + a # Permute & rescale x
  
  return(x)
}

repmat <- function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx <- dim(X)[1]
  nx <- dim(X)[2]
  return(matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T))
}
