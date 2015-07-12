shouldAlwaysFindSolutionForLinearFunctionAndRankingPref <- function() {
  discretization.methods <- getDiscretizationAlgorithmsTypes()
  for (d.method in discretization.methods) {
    alts.nr <- 10
    ranking <- sample(alts.nr)
    x <- c()
    for (i in 1:alts.nr) {
      x[ranking[i]] <- alts.nr + 1 - i
    }
    y <- ranking * 3 + 2
    perfs <- matrix(c(x, y), nrow=alts.nr)
    preferences <- t(combn(ranking, 2))
    
    solution = findSolution(perfs, preferences,
                       'SEGMENTED', charact.points.number=c(4,4), discretization.method=d.method,
                       check.consistency.only=TRUE)
    expect_that(solution$found.solution, is_equivalent_to(TRUE))
  } 
}