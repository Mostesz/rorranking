getCharacteristicPoints <- function(perfs, nums.of.characteristic.points,
                                    discretization.method = 'EQUAL_WIDTH_INTERVAL',
                                    strong.prefs = NULL,
                                    weak.prefs = NULL,
                                    indif.prefs = NULL) {
  if (is.null(discretization.method)) {
    discretization.method <- 'EQUAL_WIDTH_INTERVAL'
  }
  if (!discretization.method %in% getDiscretizationAlgorithmsTypes()) {
    stop(paste("Argument 'discretizationAlgorithmType' must be a vector of either ",
               paste(getDiscretizationAlgorithmsTypes(), collapse=", "), collapse=""))
  }
  
  switch(discretization.method,
         EQUAL_WIDTH_INTERVAL = {
           getCharacteristicPointsEqualWidthIntervalARules(perfs, nums.of.characteristic.points)
         },
         EQUAL_FREQ_INTERVAL = {
           getCharacteristicPointsEqualFreqInterval(perfs, nums.of.characteristic.points)
         },
         K_MEANS = {
           getCharacteristicPointsKMeans(perfs, nums.of.characteristic.points)
         },
         KERNEL_DENSITY_ESTIMATION = {
           getCharacteristicPointsKernelDensityEstimation(perfs, nums.of.characteristic.points)
         },
         GHADERI_DISCRETIZATION = {
           getCharacteristicPointsGhaderi(perfs, nums.of.characteristic.points, 
                                          strong.prefs, weak.prefs, indif.prefs)
         }
  )
}

getDiscretizationAlgorithmsTypes <- function() {
  return(c(
    'EQUAL_WIDTH_INTERVAL', # równa odległośc punktów charakterystycznych
    'EQUAL_FREQ_INTERVAL', # równa liczba ocen wariantów w ramach odcinków
    'K_MEANS', # k-means uruchomiony na ocenach z euklidesową metryką odległości
    'KERNEL_DENSITY_ESTIMATION', # dyskretyzacja oparta na tzw. funkcji kernelizującej
    'GHADERI_DISCRETIZATION' # metoda Ghaderi
  ));
}

getCharacteristicPointsKernelDensityEstimation <- function(perfs, nums.of.characteristic.points) {
  intervals.numbers = nums.of.characteristic.points - 1
  list.of.characteristic.points <- list()
  
  nr.crit <- ncol(perfs)
  for (crit.idx in 1:nr.crit) {
    instances <- sort(perfs[,crit.idx])
    unique.instances <- unique(instances)
    
    candidates <- list()
    for (cand.idx in 1:(length(unique.instances)-1)) {
      candidates[[cand.idx]] <- (unique.instances[cand.idx] + unique.instances[cand.idx + 1]) / 2
    }
    landmarks <- c(unique.instances[1], tail(unique.instances, n=1))
    cut.point.idx <- 0
    
    while (cut.point.idx < intervals.numbers[crit.idx] - 1) {
      Score <- c()
      
      for (cand.idx in 1:length(candidates)) {
        neighbors <- c(landmarks[sum(landmarks < candidates[[cand.idx]])],
                       landmarks[length(landmarks) - sum(landmarks > candidates[[cand.idx]]) + 1])
        
        left.inst.idx <- NULL
        n.interval <- 0
        cut.inst.idx <- NULL
        for (inst.idx in 1:length(instances)) {
          if (instances[inst.idx] > neighbors[2]) {
            break
          }
          if (instances[inst.idx] >= neighbors[1]) {
            if (is.null(left.inst.idx)) {
              left.inst.idx <- inst.idx
            }
            if (instances[inst.idx] < candidates[[cand.idx]]) {
              cut.inst.idx <- inst.idx
            }
            right.inst.idx <- inst.idx
            
            n.interval <- n.interval + 1
          }
        }
        
        f.value.left <- sum(instances >= neighbors[1] && instances < candidates[[cand.idx]]) / ((candidates[[cand.idx]] - neighbors[1]) * n.interval)
        f.value.right <- sum(instances >= candidates[[cand.idx]] && instances <= neighbors[2]) / ((neighbors[2] - candidates[[cand.idx]]) * n.interval)
        left.score <- 0
        for (i in left.inst.idx:cut.inst.idx) {
          kernel.values.sum <- 0
          for (j in left.inst.idx:right.inst.idx) {
            k = kernelFunction(instances[i], instances[j], candidates[[cand.idx]] - neighbors[1])
            kernel.values.sum <- kernel.values.sum + k
          }
          p.value = kernel.values.sum / ((candidates[[cand.idx]] - neighbors[1]) * n.interval)
          left.score <- left.score + (p.value - f.value.left)
        }
        right.score <- 0
        for (i in (cut.inst.idx+1):right.inst.idx) {
          kernel.values.sum <- 0
          for (j in left.inst.idx:right.inst.idx) {
            k = kernelFunction(instances[i], instances[j], neighbors[2] - candidates[[cand.idx]])
            kernel.values.sum <- kernel.values.sum + k
          }
          p.value = kernel.values.sum / ((neighbors[2] - candidates[[cand.idx]]) * n.interval)
          right.score <- right.score + (p.value - f.value.right)
        }
        Score[cand.idx] <- left.score + right.score
      }
      BestCandidate <- candidates[[which.max(Score)]]
      landmarks <- c(landmarks, BestCandidate)
      landmarks <- sort(landmarks)
      candidates[[which.max(Score)]] <- NULL
      
      cut.point.idx <- cut.point.idx + 1
    }
    list.of.characteristic.points[[crit.idx]] <- landmarks
  }
  return(list.of.characteristic.points)
}

kernelFunction <- function(x1,x2,h) {
  u = (x1-x2)/h
  return(kernel.function(u, kernel="gaussian"))
}

getCharacteristicPointsKMeans <- function(perfs, nums.of.characteristic.points) {
  intervals.numbers = nums.of.characteristic.points - 1
  nr.crit <- ncol(perfs)
  list.of.characteristic.points <- list()  
  for (i in 1:nr.crit) {
    crit.values <- sort(perfs[,i])
    clusters.assignments <- Kmeans(crit.values, intervals.numbers[i], iter.max=10, method = "euclidean")$cluster
    
    charac.points = c()
    last.assignment <- NULL
    for (j in 1:length(clusters.assignments)) {
      if (is.null(last.assignment)) {
        charac.points = c(charac.points, crit.values[j])
      } else if (last.assignment != clusters.assignments[j]) {
        charac.points = c(charac.points, (crit.values[j] + crit.values[j - 1]) / 2)
      }
      last.assignment <- clusters.assignments[j]
    }
    charac.points = c(charac.points, crit.values[length(clusters.assignments)])
    
    list.of.characteristic.points[[i]] = charac.points
  }
  return(list.of.characteristic.points)
}

getCharacteristicPointsGhaderi <- function(perfs, nums.of.characteristic.points,
                                           strong.prefs, weak.prefs, indif.prefs) {
  intervals.numbers = nums.of.characteristic.points - 1
  preferences = rbind(strong.prefs, weak.prefs, indif.prefs)
  
  nr.crit <- ncol(perfs)
  nr.alts <- nrow(perfs)
  levels.list <- getLevels(perfs);
  list.of.characteristic.points <- list()  
  for (i in 1:nr.crit) {
    att.values = sort(perfs[,i])
    levels <- levels.list[[i]]
    candidates <- vector(mode='numeric', length=length(levels)-1)
    for (j in 1:(length(levels)-1)) {
      candidates[j] <- (levels[j] + levels[j+1])/2
    }
    
    candidates.obj <- vector(mode='numeric', length=length(candidates))
    if (length(preferences) > 0) {
      for (preference.idx in 1:nrow(preferences)) {
        left.perfs <- perfs[preferences[preference.idx, 1], i]
        right.perfs <- perfs[preferences[preference.idx, 2], i]
        lower.bound <- min(left.perfs, right.perfs)
        upper.bound <- max(left.perfs, right.perfs)
        for (cand.idx in 1:length(candidates)) {
          if (candidates[cand.idx] >= lower.bound && candidates[cand.idx] <= upper.bound) {
            candidates.obj[cand.idx] <- candidates.obj[cand.idx] + 1
          }
        }
      }
    }
    candidates.number.constraint <- rep.int(1, length(candidates))
    
    obj <- candidates.obj
    mat <- matrix(candidates.number.constraint, nrow = 1)
    dir <- c('==')
    rhs <- c(intervals.numbers[i] - 1)
    types <- c('B')
    max <- TRUE
    lp.result <- Rglpk_solve_LP(obj, mat, dir, rhs, types = types, max = max)
    
    list.of.characteristic.points[[i]] <- c(att.values[1])
    if(lp.result$status == 0) {
      for(sol.idx in 1:length(lp.result$solution)) {
        if (lp.result$solution[sol.idx] > 0) {
          list.of.characteristic.points[[i]] <- c(list.of.characteristic.points[[i]], candidates[sol.idx])
        }
      }
    }
    list.of.characteristic.points[[i]] <- c(list.of.characteristic.points[[i]], att.values[nr.alts])
  }
  return(list.of.characteristic.points)
}

getGeneralCharacteristicPoints <- function(perfs, nums.of.characteristic.points, method.name) {
  intervals.numbers = nums.of.characteristic.points - 1
  nr.crit <- ncol(perfs)
  list.of.characteristic.points <- list()  
  for (i in 1:nr.crit) {
    list.of.characteristic.points[[i]] <- discretize(perfs[,i],
                                                     method=method.name,
                                                     categories=intervals.numbers[i],
                                                     onlycuts=TRUE)
  }
  
  return(list.of.characteristic.points)
}

getCharacteristicPointsEqualFreqInterval <- function(perfs, nums.of.characteristic.points) {
  intervals.numbers = nums.of.characteristic.points - 1
  list.of.characteristic.points <- list()
  
  nr.crit <- ncol(perfs)
  nr.alts <- nrow(perfs)
  for (i in 1:nr.crit) {
    if (length(unique(perfs[,i])) != nr.alts) {
      stop('The values of the criteria should be unique')
    }
    instances <- sort(perfs[,i])
    nrepl <- floor(nr.alts / intervals.numbers[i])
    rest <- nr.alts - nrepl * intervals.numbers[i]
    
    characteristic.points <- c(instances[1])
    instances.idx <- 1
    for(interval.idx in 1:intervals.numbers[i]) {
      instances.idx <- instances.idx + nrepl
      if (interval.idx <= rest) {
        instances.idx <- instances.idx + 1
      }
      point <- (instances[min(instances.idx, nr.alts)] + instances[instances.idx-1])/2
      characteristic.points <- c(characteristic.points, point)
    }
    
    list.of.characteristic.points[[i]] <- characteristic.points
  }
  return(list.of.characteristic.points)
}

getCharacteristicPointsEqualWidthIntervalARules <- function (perfs, nums.of.characteristic.points) {
  crit.intervals.numbers = nums.of.characteristic.points - 1
  list.of.characteristic.points = list()
  nr.crit <- ncol(perfs)
  for (i in 1:nr.crit){
    intervals.no = crit.intervals.numbers[i]
    list.of.characteristic.points[[i]] = discretize(perfs[,i], categories=intervals.no, onlycuts=TRUE)
  }
  return(list.of.characteristic.points)
}

getCharacteristicPointsEqualWidthInterval <- function (perfs, nums.of.characteristic.points) {
  # Function return list of lists of characteristic points for each criterion
  # 
  # perfs: the performance matrix
  # nums.of.characteristic.points: list of nums of characteristic points for each criterion 
  levels.list <- getLevels(perfs);
  list.of.characteristic.points = list()
  nr.crit <- ncol(perfs)
  for (i in c(1:nr.crit)){  #dla każdego kryterium
    list.of.characteristic.points[[i]] = vector(mode="numeric", length=0)
    num.of.characteristic.points = nums.of.characteristic.points[i]
    if (!is.numeric(num.of.characteristic.points)) {
      num.of.characteristic.points = 0
    }
    if (num.of.characteristic.points > 1) {
      for (j in c(1:num.of.characteristic.points)) {
        last.level = getLastElement(levels.list[[i]])      
        range.of.performance = last.level - levels.list[[i]][1]
        characteristic.point = levels.list[[i]][1] + (range.of.performance)*(j-1)/(num.of.characteristic.points-1) 
        list.of.characteristic.points[[i]] = c(list.of.characteristic.points[[i]], characteristic.point)
      }  
    } else {
      list.of.characteristic.points[[i]] = list()
    }
  }
  return(list.of.characteristic.points);
}