#' @export
isPreferred <- function(model, alternative, referenceAlternative, necessarily) {
  if (is.null(model$epsilonIndex)) {
    stop("Function isPreffered requries model with epsilon as variable. Use function buildModel with includeEpsilonAsVariable = TRUE.")
  }
  
  if (necessarily) {
    model$constraints <- combineConstraints(model$constraints,
                                            buildPairwiseComparisonConstraint(referenceAlternative, alternative, model, "strongPreference"))
  } else {
    model$constraints <- combineConstraints(model$constraints,
                                          buildPairwiseComparisonConstraint(alternative, referenceAlternative, model, "weakPreference"))
  }
  
  optimizedEpsilon <- maximizeEpsilon(model)
  
  if (necessarily) {
    return (optimizedEpsilon$status != 0 || optimizedEpsilon$optimum < model$minEpsilon)
  } else {
    return (optimizedEpsilon$status == 0 && optimizedEpsilon$optimum >= model$minEpsilon)
  }
}


#' @export
maxEpsilonSolution <- function(model, allowInconsistency = FALSE) {
  if (is.null(model$epsilonIndex)) {
    stop("Function isPreffered requries model with epsilon as variable. Use function buildModel with includeEpsilonAsVariable = TRUE.")
  }
  
  solution <- maximizeEpsilon(model)
  
  if ((solution$status == 0 && solution$optimum >= model$minEpsilon) || allowInconsistency) {
    return (toSolution(model, solution$solution))
  }
  
  return (NULL)
}

