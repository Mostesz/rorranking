generatePerformances <- function(crits.nr, alts.nr, min.value = 0, max.value = 1) {
  if (crits.nr < 2) {
    stop('Criteria number must be greater or equal 2')
  }
  
  perfs <- matrix(runif(crits.nr * alts.nr, min.value, max.value),
                  ncol=crits.nr, nrow=alts.nr)
  while(TRUE) {
    had.domination <- FALSE
    for(i in 1:(alts.nr-1)) {
      for(j in (i+1):alts.nr) {
        while(isDominated(perfs, i, j) || isDominated(perfs, j, i)) {
          had.domination <- TRUE
          perfs[i,] <- runif(crits.nr, min.value, max.value)
          perfs[j,] <- runif(crits.nr, min.value, max.value)
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