getCharacteristicPoints <- function(perf, nums.of.characteristic.points,
                                    discretization.algorithm.type = 'EQUAL_WIDTH_INTERVAL',
                                    strong.prefs = NULL,
                                    weak.prefs = NULL,
                                    indif.prefs = NULL) {
  if (!discretization.algorithm.type %in% getDiscretizationAlgorithmsTypes()) {
    stop(paste("Argument 'discretizationAlgorithmType' must be a vector of either ",
               paste(getDiscretizationAlgorithmsTypes(), collapse=", "), collapse=""))
  }
  
  switch(discretization.algorithm.type,
         EQUAL_WIDTH_INTERVAL = {
           getCharacteristicPointsEqualWidthInterval(perf, nums.of.characteristic.points)
         },
         EQUAL_FREQ_INTERVAL = {
           getCharacteristicPoints(perf, nums.of.characteristic.points, 'frequency')
         },
         K_MEANS = {
           getCharacteristicPoints(perf, nums.of.characteristic.points, 'cluster')
         },
         KERNEL_DENSITY_ESTIMATION = {
           getCharacteristicPointsKernelDensityEstimation(perf, nums.of.characteristic.points)
         },
         GHADERI_DISCRETIZATION = {
           getCharacteristicPointsGhaderi(perf, nums.of.characteristic.points, 
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

getCharacteristicPointsKernelDensityEstimation <- function(perf, nums.of.characteristic.points) {
  intervals.numbers = nums.of.characteristic.points - 1
  list.of.characteristic.points <- list()
  
  nr.crit <- ncol(perf)
  for (i in 1:nr.crit) {
    instances <- sort(perf[,i])
    unique.instances <- unique(instances)
    
    candidates <- list()
    for (j in 1:(length(unique.instances)-1)) {
      candidates[[j]] <- (unique.instances[j] + unique.instances[j+1])/2
    }
    landmarks <- c(unique.instances[1], tail(unique.instances, n=1))
    cut.point.idx <- 0
    
    while (cut.point.idx < intervals.numbers[i] - 1) {
      pl <- matrix(0, nrow=length(instances), ncol=length(instances))
      pr <- matrix(0, nrow=length(instances), ncol=length(instances))
      Score <- c()
      
      for (cand.idx in 1:length(candidates)) {
        neighbors <- c(landmarks[sum(landmarks < candidates[[cand.idx]])],
                     landmarks[length(landmarks) - sum(landmarks > candidates[[cand.idx]]) + 1])
        
        left.inst.idx <- NULL
        n.interval <- 0
        for (inst.idx in 1:length(instances)) {
          if (instances[inst.idx] > neighbors[2]) {
            break
          }
          if (instances[inst.idx] >= neighbors[1]) {
            if (is.null(left.inst.idx)) {
              left.inst.idx <- inst.idx
            }
            if (instances[inst.idx] < candidates[cand.idx]) {
              cut.inst.idx <- inst.idx
            }
            right.inst.idx <- inst.idx
            
            n.interval <- n.interval + 1
          }
        }

        fleft <- sum(instances >= neighbors[1] & instances < candidates[[cand.idx]]) / ((candidates[[cand.idx]] - neighbors[1]) * n.interval)
        fright <- sum(instances >= candidates[[cand.idx]] & instances <= neighbors[2]) / ((neighbors[2] - candidates[[cand.idx]]) * n.interval)
        
        for (m in 1:length(instances)) {
          if (instances[m] >= neighbors[1] && instances[m] < candidates[[cand.idx]]) {
            for (n in 1:length(instances)) {
              pl[n,m] <- kernelFunction(instances[m], instances[n], candidates[[cand.idx]] - neighbors[1])
              pr[n,m] <- 0
            }
          } else if (instances[m] >= candidates[[cand.idx]] && instances[m] <= neighbors[2]) {
            for (n in 1:length(instances)) {
              pr[n,m] <- kernelFunction(instances[m], instances[n], neighbors[2] - candidates[[cand.idx]])
              pl[n,m] <- 0
            }
          }
        }
        
        PLeft <- colSums(pl) / ((candidates[[cand.idx]] - neighbors[1]) * n.interval)
        PLeft <- PLeft[left.inst.idx:cut.inst.idx]
        PRight <- colSums(pr) / ((neighbors[2] - candidates[[cand.idx]]) * n.interval)
        PRight <- PRight[(cut.inst.idx+1):right.inst.idx]
        Score[cand.idx] <- abs(sum(PLeft - fleft) + sum(PRight - fright))
      }
      
      BestCandidate <- candidates[[which.max(Score)]]
      landmarks <- c(landmarks, BestCandidate)
      landmarks <- sort(landmarks)
      candidates[which.max(Score)] <- NULL
      
      cut.point.idx <- cut.point.idx + 1
    }
    list.of.characteristic.points[[i]] <- landmarks
  }
  return(list.of.characteristic.points)
}

kernelFunction <- function(x1,x2,h) {
  return(( 1 / ( (sqrt(2*pi)) * (h/6)) ) * exp( (-(x1-x2)^2) / (2*(h/6)^2) ));
}

getCharacteristicPointsKMeans <- function(perf, nums.of.characteristic.points) {
  intervals.numbers = nums.of.characteristic.points - 1
  nr.crit <- ncol(perf)
  list.of.characteristic.points <- list()  
  for (i in 1:nr.crit) {
    crit.values <- sort(perf[,i])
    clusters.assignments <- Kmeans(crit.values, intervals.numbers[i], iter.max=100, method = "euclidean")$cluster
    
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

getCharacteristicPointsGhaderi <- function(perf, nums.of.characteristic.points,
                                           strong.prefs, weak.prefs, indif.prefs) {
  intervals.numbers = nums.of.characteristic.points - 1
  preferences = rbind(strong.prefs, weak.prefs, indif.prefs)
  
  nr.crit <- ncol(perf)
  nr.alts <- nrow(perf)
  levels.list <- getLevels(perf);
  list.of.characteristic.points <- list()  
  for (i in 1:nr.crit) {
    att.values = sort(perf[,i])
    levels <- levels.list[[i]]
    candidates <- vector(mode='numeric', length=length(levels)-1)
    for (j in 1:(length(levels)-1)) {
      candidates[j] <- (levels[j] + levels[j+1])/2
    }
    
    candidates.obj <- vector(mode='numeric', length=length(candidates))
    if (length(preferences) > 0) {
      for (preference.idx in 1:nrow(preferences)) {
        left.perf <- perf[preferences[preference.idx, 1], i]
        right.perf <- perf[preferences[preference.idx, 2], i]
        for (cand.idx in 1:length(candidates)) {
          if (candidates[cand.idx] >= left.perf && candidates[cand.idx] <= right.perf) {
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

getCharacteristicPoints <- function(perf, nums.of.characteristic.points, method.name) {
  intervals.numbers = nums.of.characteristic.points - 1
  nr.crit <- ncol(perf)
  list.of.characteristic.points <- list()  
  for (i in 1:nr.crit) {
    list.of.characteristic.points[[i]] <- discretize(perf[,i],
                                                     method=method.name,
                                                     categories=intervals.numbers[i],
                                                     onlycuts=TRUE)
  }
  
  return(list.of.characteristic.points)
}

getCharacteristicPointsEqualWidthInterval <- function (perf, nums.of.characteristic.points) {
  # Function return list of lists of characteristic points for each criterion
  # 
  # perf: the performance matrix
  # nums.of.characteristic.points: list of nums of characteristic points for each criterion 
  levels.list <- getLevels(perf);
  list.of.characteristic.points = list()
  nr.crit <- ncol(perf)
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