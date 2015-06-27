findSolution <- function(perfs, preferences,
                         func.type, charact.points.number=NULL, discretization.method=NULL,
                         check.consistency.only=FALSE) {
  
  validateFuncType(func.type, charact.points.number, discretization.method)
  
  crits.nr = ncol(perfs)
  params <- getFuncParameters(func.type, charact.points.number, discretization.method, crits.nr)
  
  solution <- findPreferenceRelations(perfs, 
                                     TRUE, 
                                     strong.prefs = preferences,
                                     nums.of.characteristic.points = params$nums.of.characteristic.points,
                                     discretization.method = params$discretization.method,
                                     check.consistency.only = check.consistency.only)
  return(solution)
}

runExperiment <- function(crits.nr, alts.nr, generated.perfs.number,
                          preferences.numbers.list, repetitions.number.per.pref,
                          examined.chact.points.numbers,
                          examine.robustness=TRUE) {
  for(problem.idx in 1:generated.perfs.number) {
    perfs <- generatePerformances(crits.nr, alts.nr)
    for (preferences.number in 1:len(preferences.numbers.list)) {
      runSeries(perfs, preferences.number, examined.chact.points.numbers, examine.robustness)
    }
  }
}

runSeries <- function(perfs, preferences.number, examined.chact.points.numbers, examine.robustness) {
  preferences.models.list = getAllPreferencesModels(examined.chact.points.numbers)
  for (pref.model in preferences.models.list) {
    eps.values <- c()
    found.solutions <- 0
    for (pref.repetition.idx in 1:repetitions.number.per.pref) {
      preferences <- getPreferences(perfs, preferences.number)
      solution <- generatePreferencesAndFindSolution(perfs, preferences,
                                                     func.type=pref.model$func.type,
                                                     charact.points.number=pref.model$charact.points.number,
                                                     discretization.method=pref.model$discretization.method,
                                                     check.consistency.only=!examine.robustness)
      if (solution$found.solution) {
        eps.values <- c(eps.values, solution$eps)
        found.solutions <- found.solutions + 1
        if (examine.robustness) {
          relations.number <- sum(solution$nec.relations)
          
          pref.adj.matrix <- convertAdjacencyListToAdjacencyMatrix(preferences)
          pref.adj.matrix <- findTransitiveClosure(pref.adj.matrix)
          pref.independent.relations.number <- 0
          for (i in nrow(pref.adj.matrix)) {
            for (j in ncol(pref.adj.matrix)) {
              if (i != j && solution$nec.relations[i,j] && !pref.adj.matrix[i,j]) {
                pref.independent.relations.number <- pref.independent.relations.number + 1
              }
            }
          }
          #TODO save robustness results
          #TODO exctract to new method
        }
      }
    }
    #TODO save and return all results
  }
}

convertAdjacencyListToAdjacencyMatrix <- function(alts.number, adjacency.list) {
  result <- matrix(FALSE, nrow=alts.number, ncol=alts.number)
  for (i in 1:nrow(adjacency.list)) {
    result[adjacency.list[i, 1],adjacency.list[i, 2]] <- TRUE
  }
  return(result)
}

findTransitiveClosure <- function(matrix) {
  if (nrow(matrix) != ncol(matrix)) {
    stop('Incorrect matrix argument. Should be square matrix')
  }
  n = nrow(matrix)
  for (k in 1:n) {
    for (i in 1:n) {
      for (j in 1:n) {
        matrix[i,j] <- matrix[i,j] || (matrix[i,k] && matrix[k,j])
      }
    }
  }
  return(matrix)
}

getAllPreferencesModels <- function(examined.charact.points.numbers) {
  result <- list()
  pref.idx = 1
  for (func.type in c('LINEAR', 'GENERAL', 'SEGMENTED')) {
    if (func.type == 'SEGMENTED') {
      for (discretization.method in getDiscretizationAlgorithmsTypes()) {
        for (charact.points in examined.charact.points.numbers) {
          result[[pref.idx]] <- list(func.type=func.type,
                                     charact.points.number=charact.points,
                                     discretization.method=discretization.method)
          pref.idx = pref.idx + 1
        }
      }
    } else {
      result[[pref.idx]] <- list(func.type=func.type)
      pref.idx = pref.idx + 1
    }
  }
  return(result)
}

getPreferences <- function(perfs, preferences.number) {
  if (is.null(preferences.number)) {
    return(generateRankingPreferences(perfs))
  } else {
    return(generatePreferences(perfs, preferences.number))
  }
}

getFuncParameters <- function(func.type, charact.points.number, discretization.method, crits.nr) {
  if (func.type == 'LINEAR') {
    return(list(nums.of.characteristic.points=rep(2, crits.nr), discretization.method=NULL))
  } else if (func.type == 'GENERAL') {
    return(list(nums.of.characteristic.points=NULL, discretization.method=NULL))
  } else if (func.type == 'SEGMENTED') {
    return(list(nums.of.characteristic.points=charact.points.number, discretization.method=discretization.method))
  }
}

validateFuncType <- function(func.type, charact.points.number, discretization.method) {
  if (func.type == 'LINEAR' || func.type == 'GENERAL') {
    if (any(!is.null(charact.points.number), !is.null(discretization.method))) {
      stop(paste('Characteristic points number and discretization methods should not be undefined for', func.type,'function'))
    }
  } else if (func.type == 'SEGMENTED') {
    if (any(is.null(charact.points.number), is.null(discretization.method))) {
      stop('Characteristic points number and discretization methods should be defined for SEGMENTED function')
    }
  } else {
    stop('Invalid function type')
  }
}