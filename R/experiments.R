generateProblemAndFindSolution <- function(crits.nr, alts.nr, preferences.number=NULL, func.type,
                                           charact.points.number=NULL, discretization.method=NULL,
                                           check.consistency.only=FALSE) {
  validateFuncType(func.type, charact.points.number, discretization.method)
  
  perfs <- generatePerformances(crits.nr, alts.nr)
  preferences <- getPreferences(perfs, preferences.number)
  
  params <- getFuncParameters(func.type, charact.points.number, discretization.method, crits.nr)
  
  solution <- findPreferenceRelations(perfs, 
                                     TRUE, 
                                     strong.prefs = preferences,
                                     nums.of.characteristic.points = params$nums.of.characteristic.points,
                                     discretization.method = params$discretization.method,
                                     check.consistency.only = check.consistency.only)
  return(solution)
}

runExpressivenessExperiment <- function(crits.nr, alts.nr, preference.type, func.type,
                                        charact.points.number=NULL, discretization.method=NULL, preferences.number=NULL) {
  
}

runResistanceExperiment <- function(crits.nr, alts.nr, preference.type, func.type,
                                    charact.points.number=NULL, discretization.method=NULL, preferences.number=NULL) {
  
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
  } else if (func.type == 'SERIALIZED') {
    return(list(nums.of.characteristic.points=charact.points.number, discretization.method=discretization.method))
  }
}

validateFuncType <- function(func.type, charact.points.number, discretization.method) {
  if (func.type == 'LINEAR' || func.type == 'GENERAL') {
    if (any(!is.null(charact.points.number), !is.null(discretization.method))) {
      stop(paste('Characteristic points number and discretization methods should not be undefined for', func.type,'function'))
    }
  } else if (func.type == 'SERIALIZED') {
    if (any(is.null(charact.points.number), is.null(discretization.method))) {
      stop('Characteristic points number and discretization methods should be defined for SERIALIZED function')
    }
  } else {
    stop('Invalid function type')
  }
}