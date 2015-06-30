runExperiment <- function(results.base.dir,
                          matrix.sizes, generating.perfs.number,
                          preferences.numbers.list, repetitions.number.per.pref,
                          examined.chact.points.numbers,
                          examine.robustness=TRUE) {
  if (!file.exists(results.base.dir)) {
    stop('Given base directory does not exist')
  } else {
    for (dir.to.rm in file.path(results.base.dir, list.files(results.base.dir))) {
      unlink(dir.to.rm, recursive=TRUE)
    }
  }
  
  start.time <- Sys.time()
  
  for (m.size.idx in 1:nrow(matrix.sizes)) {
    crits.nr <- matrix.sizes[m.size.idx,1]
    alts.nr <- matrix.sizes[m.size.idx,2]
    matrix.size.dir.name <- paste(crits.nr, alts.nr, sep='x')
    
    for (preferences.number in preferences.numbers.list) {
      preferences.type.dir.name <- if (is.null(preferences.number)) 'ranking' else paste('pref', preferences.number, sep='-')
      
      preferences.models.list = getAllPreferencesModels(examined.chact.points.numbers)
      for (pref.model in preferences.models.list) {
        pref.model.dir.name <- pref.model$func.type
        if (pref.model$func.type == 'SEGMENTED') {
          pref.model.dir.name <- paste(pref.model$func.type, pref.model$charact.points.number, pref.model$discretization.method, sep='-')
        }
        pref.model.dir.path <- file.path(results.base.dir, matrix.size.dir.name, preferences.type.dir.name, pref.model.dir.name)
        dir.create(pref.model.dir.path, recursive=TRUE)
        
        series.result <- runSeries(crits.nr, alts.nr, generating.perfs.number,
                                   preferences.number, repetitions.number.per.pref,
                                   pref.model,
                                   examine.robustness)
        write.csv(series.result$eps.values, file=file.path(pref.model.dir.path, 'eps.values'))
        write.csv(series.result$found.solutions, file=file.path(pref.model.dir.path, 'found.solutions'))
        write.csv(series.result$relations.numbers, file=file.path(pref.model.dir.path, 'relations.numbers'))
        write.csv(series.result$pref.ind.relations.numbers, file=file.path(pref.model.dir.path, 'pref.ind.relations.numbers'))
      }
    }
  }
  
  print(Sys.time() - start.time)
}

runSeries <- function(crits.nr, alts.nr, generating.perfs.number,
                      preferences.number, repetitions.number.per.pref,
                      pref.model,
                      examine.robustness) {
  perfs.names <- paste('m', 1:generating.perfs.number, sep='')
  
  eps.values.df <- data.frame(matrix(ncol=length(perfs.names), nrow=repetitions.number.per.pref))
  colnames(eps.values.df) <- perfs.names
  found.solutions.df <- data.frame(matrix(0, ncol=length(perfs.names), nrow=1))
  colnames(found.solutions.df) <- perfs.names
  relations.numbers.df <- data.frame(matrix(ncol=length(perfs.names), nrow=repetitions.number.per.pref))
  colnames(relations.numbers.df) <- perfs.names
  pref.ind.relations.numbers.df <- data.frame(matrix(ncol=length(perfs.names), nrow=repetitions.number.per.pref))
  colnames(pref.ind.relations.numbers.df) <- perfs.names
  for(perfs.name in perfs.names) {
    perfs <- generatePerformances(crits.nr, alts.nr)
    
    for (pref.repetition.idx in 1:repetitions.number.per.pref) {
      preferences <- getPreferences(perfs, preferences.number)
      solution <- findSolution(perfs, preferences,
                               func.type=pref.model$func.type,
                               charact.points.number=rep(pref.model$charact.points.number, crits.nr),
                               discretization.method=pref.model$discretization.method,
                               check.consistency.only=!examine.robustness)
      if (solution$found.solution) {
        eps.values.df[pref.repetition.idx, perfs.name] <- solution$eps
        found.solutions.df[1, perfs.name] <- found.solutions.df[1, perfs.name] + 1
        if (examine.robustness) {
          relations.numbers.df[pref.repetition.idx, perfs.name] <- sum(solution$nec.relations)
          
          pref.adj.matrix <- convertAdjacencyListToAdjacencyMatrix(preferences, alts.nr)
          pref.adj.matrix <- findTransitiveClosure(pref.adj.matrix)
          pref.independent.relations.number <- 0
          for (i in nrow(pref.adj.matrix)) {
            for (j in ncol(pref.adj.matrix)) {
              if (i != j && solution$nec.relations[i,j] && !pref.adj.matrix[i,j]) {
                pref.independent.relations.number <- pref.independent.relations.number + 1
              }
            }
          }
          pref.ind.relations.numbers.df[pref.repetition.idx, perfs.name] <- pref.independent.relations.number
        }
      }
    }
  }
  return(list(eps.values=eps.values.df,
              found.solutions=found.solutions.df,
              relations.numbers=relations.numbers.df,
              pref.ind.relations.numbers=pref.ind.relations.numbers.df))
}

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

convertAdjacencyListToAdjacencyMatrix <- function(adjacency.list, values.nr) {
  result <- matrix(FALSE, nrow=values.nr, ncol=values.nr)
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