#### HELPERS

buildPairwiseComparisonConstraint <- function(alternative, referenceAlternative, model, type) {
  stopifnot(type %in% c("weakPreference", "strongPreference", "indifference"))
  
  lhs <- ua(referenceAlternative, ncol(model$constraints$lhs), model$perfToModelVariables) - ua(alternative, ncol(model$constraints$lhs), model$perfToModelVariables)
  dir <- "<="
  rhs <- 0
  
  if (type == "strongPreference") {
    if (is.null(model$epsilonIndex)) {
      rhs <- -model$minEpsilon
    } else {
      lhs[model$epsilonIndex] <- 1
    }
  } else if (type == "indifference") {
    dir <- "=="
  }
  
  return (list(lhs = lhs, dir = dir, rhs = rhs))
}


addVarialbesToModel <- function(constraints, variables) {
  for (var in variables)
    constraints$lhs <- cbind(constraints$lhs, 0)
  constraints$types <- c(constraints$types, variables)
  return (constraints)
}

combineConstraints <- function(...) {
  allConst <- list(...)
  
  lhs <- c()
  dir <- c()
  rhs <- c()
  types <- c()
  
  for (const in allConst) {
    if (!is.null(const)) {      
      lhs <- rbind(lhs, const$lhs)
      dir <- c(dir, const$dir)
      rhs <- c(rhs, const$rhs)
      types <- c(types, const$types)
    }
  }
  
  return (list(lhs = lhs, dir = dir, rhs = rhs, types = types))
}

removeConstraints <- function(allConst, constraintsToRemoveIndices) {
  return (list(lhs = allConst$lhs[-c(constraintsToRemoveIndices), ],
               dir = allConst$dir[-c(constraintsToRemoveIndices)],
               rhs = allConst$rhs[-c(constraintsToRemoveIndices)],
               types = allConst$types))
}

#### BUILDING MODEL

#' @export
buildModel <- function(problem, minEpsilon = 1e-4) { # includeEpsilonAsVariable,
  nrAlternatives <- nrow(problem$perf)
  nrCriteria <- ncol(problem$perf)
  
  # criterion value to alternative indices
  
  criterionValues <- replicate(nrCriteria, list())
  
  for (j in seq_len(nrCriteria)) {
    for (i in seq_len(nrAlternatives)) {
      value <- problem$perf[i, j]
      
      found <- FALSE
      
      for (k in seq_len(length(criterionValues[[j]]))) {
        if (criterionValues[[j]][[k]]$value == value) { # todo: consider epsilon
          found <- TRUE
          criterionValues[[j]][[k]]$alternatives <- c(criterionValues[[j]][[k]]$alternatives, i)
        }
      }
      
      if (!found) {
        criterionValues[[j]][[length(criterionValues[[j]]) + 1]] <- list(
          value=value,
          alternatives=c(i)
        )
      }
    }
    
    if (length(criterionValues[[j]]) < 2) {
      stop(paste("Criterion ", j, " is superfluous!"))
    }
    
    # sort criterion values
    criterionValues[[j]] <- criterionValues[[j]][order(
      sapply(criterionValues[[j]],
             function(x) x$value, simplify=TRUE
      ), decreasing=FALSE)]
  }
  
  perfToModelVariables <- replicate(nrCriteria, replicate(nrAlternatives, list()))
  firstChPointVariableIndex <- c(1)
  chPoints <- c()
  
  for (j in seq_len(nrCriteria)) {
    numberOfCharacteristicPoints <- problem$characteristicPoints[j]
    
    if (numberOfCharacteristicPoints == 0) {
      numberOfCharacteristicPoints <- length(criterionValues[[j]])
    }
    
    if (j != nrCriteria) {
      firstChPointVariableIndex[j + 1] <- firstChPointVariableIndex[j] + numberOfCharacteristicPoints - 1
    }
    
    chPoints[j] <- numberOfCharacteristicPoints
  }
  
  numberOfVariables <- firstChPointVariableIndex[length(firstChPointVariableIndex)] + chPoints[nrCriteria] - 2
  
  for (j in seq_len(nrCriteria)) {
    firstValue <- criterionValues[[j]][[1]]$value
    lastValue <- criterionValues[[j]][[length(criterionValues[[j]])]]$value
    direction <- problem$criteria[j]
    
    if (problem$characteristicPoints[j] == 0) {      
      for (i in seq_len(nrAlternatives)) {
        value <- problem$perf[i, j]
        criterionValueIndex <- which(sapply(criterionValues[[j]], function(x){x$value == value}))
        
        if (direction == "g" && criterionValueIndex > 1) {
          perfToModelVariables[[i, j]][[1]] = c(firstChPointVariableIndex[j] + criterionValueIndex - 2, 1.0)
        } else if (direction == "c" && criterionValueIndex < length(criterionValues[[j]])) {
          perfToModelVariables[[i, j]][[1]] = c(firstChPointVariableIndex[j] + criterionValueIndex - 1, 1.0)
        }
      }
    } else {
      numberOfCharacteristicPoints <- problem$characteristicPoints[j]
      intervalLength <- (lastValue - firstValue) / (numberOfCharacteristicPoints - 1);
      coeff <- 1.0 / intervalLength;
      
      for (i in seq_len(nrAlternatives)) {
        value <- problem$perf[i, j]
        
        if (direction == "g") {
          if (value == lastValue) {
            perfToModelVariables[[i, j]][[1]] <- c(firstChPointVariableIndex[j] + numberOfCharacteristicPoints - 2, 1.0)
          } else if (value > firstValue) {
            lowerChPointIndex <- floor((value - firstValue) * coeff)
            
            if (lowerChPointIndex >= numberOfCharacteristicPoints - 1) {
              stop("InternalError?: lowerChPointIndex >= numberOfCharacteristicPoints - 1: This should never happen.");
            }
            
            lowerValue = firstValue + intervalLength * lowerChPointIndex
            upperValue = firstValue + intervalLength * (lowerChPointIndex + 1)
            
            lowerCoeff <- 0.0
            upperCoeff <- 0.0
            
            if (value <= lowerValue) {
              # comp accuracy
              lowerCoeff = 1.0
              upperCoeff = 0.0
            } else if (value >= upperValue) {
              # comp accuracy
              lowerCoeff = 0.0
              upperCoeff = 1.0
            } else {
              lowerCoeff = (lowerValue - value) / (upperValue - lowerValue) + 1.0
              upperCoeff = (value - lowerValue) / (upperValue - lowerValue)
            }
            
            if (lowerChPointIndex > 0) {
              perfToModelVariables[[i, j]][[1]] = c(firstChPointVariableIndex[j] + lowerChPointIndex - 1, lowerCoeff)
              perfToModelVariables[[i, j]][[2]] = c(firstChPointVariableIndex[j] + lowerChPointIndex, upperCoeff)
            } else {
              perfToModelVariables[[i, j]][[1]] = c(firstChPointVariableIndex[j] + lowerChPointIndex, upperCoeff)
            }
          }
        } else {
          if (value == firstValue) {
            perfToModelVariables[[i, j]][[1]] = c(firstChPointVariableIndex[j], 1.0)
          } else if (value < lastValue) {
            lowerChPointIndex <- floor((value - firstValue) * coeff)
            
            if (lowerChPointIndex >= numberOfCharacteristicPoints - 1) {
              stop("InternalError?: lowerChPointIndex >= numberOfCharacteristicPoints - 1: This should never happen.");
            }
            
            lowerValue = firstValue + intervalLength * lowerChPointIndex
            upperValue = firstValue + intervalLength * (lowerChPointIndex + 1)
            
            lowerCoeff <- 0.0
            upperCoeff <- 0.0
            
            if (value <= lowerValue) {
              # comp accuracy
              lowerCoeff = 1.0
              upperCoeff = 0.0
            } else if (value >= upperValue) {
              # comp accuracy
              lowerCoeff = 0.0
              upperCoeff = 1.0
            } else {
              lowerCoeff = (upperValue - value) / (upperValue - lowerValue)
              upperCoeff = (value - upperValue) / (upperValue - lowerValue) + 1.0
            }
            
            if (lowerChPointIndex < numberOfCharacteristicPoints - 2) {
              perfToModelVariables[[i, j]][[1]] = c(firstChPointVariableIndex[j] + lowerChPointIndex, lowerCoeff)
              perfToModelVariables[[i, j]][[2]] = c(firstChPointVariableIndex[j] + lowerChPointIndex + 1, upperCoeff)
            } else {
              perfToModelVariables[[i, j]][[1]] = c(firstChPointVariableIndex[j] + lowerChPointIndex, lowerCoeff)
            }
          }
        }
      }
    }
  }
  
  # epsilon index
  
  #epsilonIndex <- NULL
  #if (includeEpsilonAsVariable) {
  numberOfVariables <- numberOfVariables + 1
  epsilonIndex <- numberOfVariables
  #}
  
  # constraints
  
  # sum to 1
  
  lhs <- rep(0, numberOfVariables)
  
  for (j in seq_len(nrCriteria)) {
    if (problem$criteria[j] == 'g')
      lhs[firstChPointVariableIndex[j] + chPoints[j] - 2] <- 1
    else
      lhs[firstChPointVariableIndex[j]] <- 1
  }
  
  constraints <- list(lhs = lhs, dir = "==", rhs = 1)
  
  # monotonicity of vf
  
  for (j in seq_len(nrCriteria)) {
    for (k in seq_len(chPoints[j] - 2)) {
      lhs <- rep(0, numberOfVariables)
      rhs <- 0
      
      if (problem$criteria[j] == "g") {
        lhs[firstChPointVariableIndex[j] + k - 1] <- 1
        lhs[firstChPointVariableIndex[j] + k] <- -1
      } else {
        lhs[firstChPointVariableIndex[j] + k - 1] <- -1
        lhs[firstChPointVariableIndex[j] + k] <- 1
      }
      
      if (problem$strictVF) {
        #if (includeEpsilonAsVariable) {
        lhs[epsilonIndex] <- 1
        #} else {
        #  rhs <- -minEpsilon
        #}
      }
      
      constraints <- combineConstraints(constraints,
                                        list(lhs = lhs, dir = "<=", rhs = rhs))
    }
    
    lhs <- rep(0, numberOfVariables)
    rhs <- 0
    
    if (problem$criteria[j] == 'g')
      lhs[firstChPointVariableIndex[j]] <- -1
    else
      lhs[firstChPointVariableIndex[j] + chPoints[j] - 2] <- -1
    
    if (problem$strictVF) {
      #if (includeEpsilonAsVariable) {
      lhs[epsilonIndex] <- 1
      #} else {
      #  rhs <- -minEpsilon
      #}
    }
    
    constraints <- combineConstraints(constraints,
                                      list(lhs = lhs, dir = "<=", rhs = rhs))
  }
  
  constraints$types <- rep("C", numberOfVariables)
  
  # building model
  
  model <- list(
    constraints = constraints,
    firstChPointVariableIndex = firstChPointVariableIndex,
    epsilonIndex = epsilonIndex,
    chPoints = chPoints,
    perfToModelVariables = perfToModelVariables,
    criterionValues = criterionValues,
    criterionPreferenceDirection = problem$criteria,
    prefInfoToConstraints = list(),
    generalVF = problem$characteristicPoints == 0,
    minEpsilon = minEpsilon
  )
  
  # preference information
  
  # assignment examples
  
  prefInfoIndex <- 1
  
  if (is.matrix(problem$strongPreference)) {
    for (k in seq_len(nrow(problem$strongPreference))) {
      alternative <- problem$strongPreference[k, 1]
      referenceAlternative <- problem$strongPreference[k, 2]
      
      model$constraints <- combineConstraints(model$constraints,
                                              buildPairwiseComparisonConstraint(alternative, referenceAlternative,
                                                                                model, type = "strongPreference"))
      
      model$prefInfoToConstraints[[prefInfoIndex]] <- nrow(model$constraints$lhs)
      prefInfoIndex <- prefInfoIndex + 1
    }
  }
  
  if (is.matrix(problem$weakPreference)) {
    for (k in seq_len(nrow(problem$weakPreference))) {
      alternative <- problem$weakPreference[k, 1]
      referenceAlternative <- problem$weakPreference[k, 2]
      
      model$constraints <- combineConstraints(model$constraints,
                                              buildPairwiseComparisonConstraint(alternative, referenceAlternative,
                                                                                model, type = "weakPreference"))
      
      model$prefInfoToConstraints[[prefInfoIndex]] <- nrow(model$constraints$lhs)
      prefInfoIndex <- prefInfoIndex + 1
    }
  }
  
  if (is.matrix(problem$indifference)) {
    for (k in seq_len(nrow(problem$indifference))) {
      alternative <- problem$indifference[k, 1]
      referenceAlternative <- problem$indifference[k, 2]
      
      model$constraints <- combineConstraints(model$constraints,
                                              buildPairwiseComparisonConstraint(alternative, referenceAlternative,
                                                                                model, type = "indifference"))
      
      model$prefInfoToConstraints[[prefInfoIndex]] <- nrow(model$constraints$lhs)
      prefInfoIndex <- prefInfoIndex + 1
    }
  }
  
  return (model)
}

ua <- function(alternative, nrVariables, perfToModelVariables) {
  res <- rep(0, nrVariables)
  
  for (j in seq_len(ncol(perfToModelVariables))) {
    for (k in seq_len(length(perfToModelVariables[[alternative, j]]))) {
      res[perfToModelVariables[[alternative, j]][[k]][1]] <- perfToModelVariables[[alternative, j]][[k]][2]
    }
  }
  
  return (res)
}

eliminateEpsilon <- function(model) {
  stopifnot(!is.null(model$epsilonIndex))
  
  model$constraints$rhs <- model$constraints$rhs - model$constraints$lhs[, model$epsilonIndex] * model$minEpsilon
  model$constraints$lhs <- model$constraints$lhs[, -c(model$epsilonIndex)]
  model$constraints$types <- model$constraints$types[-c(model$epsilonIndex)]
  model$epsilonIndex <- NULL
  
  return (model)
}
