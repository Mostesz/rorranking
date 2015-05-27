generatePerformances <- function(crits.nr, alts.nr, min.value = 0, max.value = 1) {
  if (crits.nr < 2) {
    stop('Criteria number must be greater or equal 2')
  }
  
  perfs = matrix(runif(crits.nr * alts.nr, min.value, max.value),
                  ncol=crits.nr, nrow=alts.nr)
  while(TRUE) {
    has.domination = FALSE
    for(i in 1:(alts.nr-1)) {
      for(j in (i+1):alts.nr) {
        while(isDominated(perfs, i, j) || isDominated(perfs, j, i)) {
          has.domination = TRUE
          perfs[i,] = runif(crits.nr, min.value, max.value)
          perfs[j,] = runif(crits.nr, min.value, max.value)
        }
      }
    }
    if (!has.domination) {
      break
    }
  }
  
  return(perfs)
}

isDominated <- function(perfs, dominant.idx, dominated.idx) {
  for(i in 1:ncol(perfs)) {
    if (perfs[dominated.idx, i] > perfs[dominant.idx, i]) {
      return(FALSE)
    }
  }
  return(TRUE)
}