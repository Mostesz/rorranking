launchAllExperiments <- function(results.base.dir,
                                 matrix.sizes, generating.perfs.number,
                                 preferences.numbers.list, pref.repetitions.number.expr, pref.repetitions.number.rbst,
                                 examined.chact.points.numbers) {
  for(perf.idx in 1:generating.perfs.number) {
    print(paste('/////////////// PROBLEM REPETITION NO:', perf.idx))
    perfs <- NULL
    
    perfs.processing.time <- Sys.time()
    results <- list()
    for (distribution in c('UNIFORM', 'SKEW_NORMAL')){
      matrix.sizes.results <- list()
      for (m.size.idx in 1:nrow(matrix.sizes)) {
        single.problem.processing.time <- Sys.time()
        
        crits.nr <- matrix.sizes[m.size.idx,1]
        alts.nr <- matrix.sizes[m.size.idx,2]
        if(is.null(perfs)) {
          perfs <- generatePerformances(crits.nr, alts.nr, distribution)
        }
        matrix.size.dir.name <- paste(crits.nr, alts.nr, sep='x')
        
        print(paste(matrix.size.dir.name, distribution))
        
        series.path <- file.path(results.base.dir, matrix.size.dir.name, distribution)
        
        res <- launchAllExperimentsForPerfs(perfs,
                                            preferences.numbers.list, pref.repetitions.number.expr, pref.repetitions.number.rbst,
                                            examined.chact.points.numbers)
        matrix.sizes.results[[matrix.size.dir.name]] <- res
        
        print('Matrix size execution time:')
        print(Sys.time() - single.problem.processing.time)
      }
      results[[distribution]] <- matrix.sizes.results
    }
    print('Results saving... Please do not terminate the script.')
    saveResults(results.base.dir, results)
    print('Saving finished.')
    print('Repetition execution time:')
    print(Sys.time() - perfs.processing.time)
  }
}

saveResultsForPerfs <- function(results.base.dir, results) {
  for(distributions.dir in names(results)) {
    distributions.results <- results[[distributions.dir]]
    for(m.sizes.dir in names(distributions.results)) {
      m.sizes.results <- distributions.results[[m.sizes.dir]]
      for(pref.type.dir in names(m.sizes.results)) {
        pref.type.results <- m.sizes.results[[pref.type.dir]]
        for(pref.model.dir in names(pref.type.results)) {
          results <- pref.type.results[pref.model.dir]
          path <- file.path(results.base.dir, distributions.dir, m.sizes.dir, pref.type.dir, pref.model.dir)
          if (!file.exists(path)) {
            dir.create(path, recursive=TRUE)
          }
          write.table(results$found.solutions.number, file=file.path(path, 'foundsolutionsnumber.csv'), append=T, row.names=F, col.names=F, sep=',')
          write.table(results$eps.values, file=file.path(path, 'epsvalues.csv'), append=T, row.names=F, col.names=F, sep=',')
          if (!is.null(results$relations.numbers)) {
            write.table(results$relations.numbers, file=file.path(path, 'relationsnumbers.csv'), append=T, row.names=F, col.names=F, sep=',')
          }
          if (!is.null(results$pref.ind.relations.numbers)) {
            write.table(results$pref.ind.relations.numbers, file=file.path(path, 'prefindrelationsnumbers.csv'), append=T, row.names=F, col.names=F, sep=',')
          }
        }
      }
    }
  }
}

launchAllExperimentsForPerfs <- function(perfs,
                                         preferences.numbers.list, pref.repetitions.number.expr, pref.repetitions.number.rbst,
                                         examined.chact.points.numbers) {
  preferences.types.results <- list()
  for (preferences.number in preferences.numbers.list) {
    if (!is.null(preferences.number) && nrow(perfs)*(nrow(perfs)-1)/2 < preferences.number) {
      break
    }
    
    preferences.type.dir.name <- if (is.null(preferences.number)) 'ranking' else paste('pref', preferences.number, sep='-')
    preferences.models.list = getAllPreferencesModels(examined.chact.points.numbers, ncol(perfs))
    
    models.results <- list()
    for (pref.model in preferences.models.list) {
      pref.model.dir.name <- pref.model$func.type
      if (pref.model$func.type == 'SEGMENTED') {
        pref.model.dir.name <- paste(pref.model$func.type, pref.model$charact.points.number, pref.model$discretization.method, sep='-')
      }
      print(paste('>>', paste(pref.model.dir.name, preferences.type.dir.name, sep=', ')))
      
      exp.res <- launchExpressivenessExperiment(perfs,
                                                preferences.number, pref.repetitions.number.expr,
                                                pref.model)
      rob.res <- launchRobustnessExperiment(perfs,
                                            preferences.number, pref.repetitions.number.expr, pref.repetitions.number.rbst,
                                            pref.model)
      res <- rob.res
      res$found.solutions.number <- exp.res
      
      models.results[[pref.model.dir.name]] <- res
    }
    preferences.types.results[[preferences.type.dir.name]] <- models.results
  }
  return(preferences.types.results)
}

launchExpressivenessExperiment <- function(perfs,
                                           preferences.number, pref.repetitions.number.expr,
                                           pref.model) {
  found.solutions.number <- 0
  for (pref.repetition.idx in 1:pref.repetitions.number.expr) {
    preferences <- getPreferences(perfs, preferences.number)
    solution <- findSolution(perfs, preferences,
                             func.type=pref.model$func.type,
                             charact.points.number=rep(pref.model$charact.points.number, ncol(perfs)),
                             discretization.method=pref.model$discretization.method,
                             check.consistency.only=TRUE)
    if (solution$found.solution) {
      found.solutions.number <- found.solutions.number + 1
    }
  }
  return(found.solutions.number)
}

launchRobustnessExperiment <- function(perfs,
                                       preferences.number, pref.repetitions.number.expr, pref.repetitions.number.rbst,
                                       pref.model) {
  eps.values <- matrix(ncol=pref.repetitions.number.expr, nrow=1)
  relations.numbers <- NULL
  pref.ind.relations.numbers <- NULL
  
  for (pref.repetition.idx in 1:pref.repetitions.number.expr) {
    examine.robustness <- pref.repetition.idx <= pref.repetitions.number.rbst && !is.null(preferences.number)
    
    preferences <- getPreferences(perfs, preferences.number, using.linear.func=TRUE)
    solution <- findSolution(perfs, preferences,
                             func.type=pref.model$func.type,
                             charact.points.number=rep(pref.model$charact.points.number, ncol(perfs)),
                             discretization.method=pref.model$discretization.method,
                             check.consistency.only=!examine.robustness)
    if (solution$found.solution) {
      eps.values[1, pref.repetition.idx] <- solution$eps
      
      if (examine.robustness) {
        relations.numbers <- matrix(ncol=pref.repetitions.number.rbst, nrow=1)
        pref.ind.relations.numbers <- matrix(ncol=pref.repetitions.number.rbst, nrow=1)
        
        rel.no <- computeRelationsNumbers(solution$nec.relations, preferences, nrow(perfs))
        relations.numbers[1, pref.repetition.idx] <- rel.no$relations.numbers
        pref.ind.relations.numbers[1, pref.repetition.idx] <- rel.no$pref.independent.relations.number
      }
    }
  }
  
  return(list(
    eps.values=eps.values,
    relations.numbers=relations.numbers,
    pref.ind.relations.numbers=pref.ind.relations.numbers
    ))
}

runExperiment <- function(results.base.dir,
                          matrix.sizes, generating.perfs.number,
                          preferences.numbers.list, repetitions.number.per.pref,
                          examined.chact.points.numbers,
                          distribution = 'UNIFORM',
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
    print(matrix.size.dir.name)
    
    for (preferences.number in preferences.numbers.list) {
      preferences.type.dir.name <- if (is.null(preferences.number)) 'ranking' else paste('pref', preferences.number, sep='-')
      
      preferences.models.list = getAllPreferencesModels(examined.chact.points.numbers, crits.nr)
      for (pref.model in preferences.models.list) {
        series.start.time <- Sys.time()
        
        pref.model.dir.name <- pref.model$func.type
        if (pref.model$func.type == 'SEGMENTED') {
          pref.model.dir.name <- paste(pref.model$func.type, pref.model$charact.points.number, pref.model$discretization.method, sep='-')
        }
        pref.model.dir.path <- file.path(results.base.dir, matrix.size.dir.name, preferences.type.dir.name, pref.model.dir.name)
        dir.create(pref.model.dir.path, recursive=TRUE)
        
        series.result <- runSeries(crits.nr, alts.nr, generating.perfs.number,
                                   preferences.number, repetitions.number.per.pref,
                                   pref.model,
                                   distribution,
                                   examine.robustness)
        write.csv(series.result$eps.values, file=file.path(pref.model.dir.path, 'eps.values'))
        write.csv(series.result$found.solutions, file=file.path(pref.model.dir.path, 'found.solutions'))
        if (examine.robustness) {
          write.csv(series.result$relations.numbers, file=file.path(pref.model.dir.path, 'relations.numbers'))
          write.csv(series.result$pref.ind.relations.numbers, file=file.path(pref.model.dir.path, 'pref.ind.relations.numbers'))
        }
        
        print('Series finished')
        print(Sys.time() - series.start.time)
      }
    }
  }
  
  print(Sys.time() - start.time)
}

runSeries <- function(crits.nr, alts.nr, generating.perfs.number,
                      preferences.number, repetitions.number.per.pref,
                      pref.model,
                      distribution,
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
    perfs <- generatePerformances(crits.nr, alts.nr, distribution)
    
    found.sol.idx <- 1
    for (pref.repetition.idx in 1:repetitions.number.per.pref) {
      preferences <- getPreferences(perfs, preferences.number)
      solution <- findSolution(perfs, preferences,
                               func.type=pref.model$func.type,
                               charact.points.number=rep(pref.model$charact.points.number, crits.nr),
                               discretization.method=pref.model$discretization.method,
                               check.consistency.only=!examine.robustness)
      if (solution$found.solution) {
        eps.values.df[found.sol.idx, perfs.name] <- solution$eps
        found.solutions.df[1, perfs.name] <- found.solutions.df[1, perfs.name] + 1
        if (examine.robustness) {
          robustness.relations.numbers <- computeRelationsNumbers(solution$nec.relations, preferences, alts.nr)
          relations.numbers.df[found.sol.idx, perfs.name] <- robustness.relations.numbers$relations.numbers
          pref.ind.relations.numbers.df[found.sol.idx, perfs.name] <- robustness.relations.numbers$pref.independent.relations.number
        }
        found.sol.idx <- found.sol.idx + 1
      }
    }
  }
  return(list(eps.values=eps.values.df,
              found.solutions=found.solutions.df,
              relations.numbers=relations.numbers.df,
              pref.ind.relations.numbers=pref.ind.relations.numbers.df))
}

computeRelationsNumbers <- function(nec.relations, preferences, alts.nr) {
  relations.numbers <- sum(nec.relations)
  
  pref.adj.matrix <- convertAdjacencyListToAdjacencyMatrix(preferences, alts.nr)
  pref.adj.matrix <- findTransitiveClosure(pref.adj.matrix)
  pref.independent.relations.number <- 0
  for (i in nrow(pref.adj.matrix)) {
    for (j in ncol(pref.adj.matrix)) {
      if (i != j && nec.relations[i,j] && !pref.adj.matrix[i,j]) {
        pref.independent.relations.number <- pref.independent.relations.number + 1
      }
    }
  }
  return(list(
    relations.numbers=relations.numbers,
    pref.independent.relations.number=pref.independent.relations.number))
}

findSolution <- function(perfs, preferences,
                         func.type, charact.points.number=NULL, discretization.method=NULL,
                         check.consistency.only=FALSE) {
  
  validateFuncType(func.type, charact.points.number, discretization.method)
  
  crits.nr = ncol(perfs)
  params <- getFuncParameters(func.type, charact.points.number, discretization.method, crits.nr)
  
  solution <- findPreferenceRelations(perfs, 
                                      FALSE, 
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

getAllPreferencesModels <- function(examined.charact.points.numbers, crits.nr) {
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

getPreferences <- function(perfs, preferences.number, using.linear.func=FALSE) {
  if (is.null(preferences.number)) {
    return(generateRankingPreferences(perfs))
  } else if (using.linear.func) {
    return(generatePreferencesFromLinearFunc(perfs, preferences.number))
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