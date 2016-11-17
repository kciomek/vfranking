#' @import Rglpk
extremizeVariable <- function(constraints, variableIndex, maximize) {
  obj <- rep(0, ncol(constraints$lhs))
  obj[variableIndex] <- 1  
  Rglpk_solve_LP(obj, constraints$lhs, constraints$dir, constraints$rhs, max = maximize,
                 types = constraints$types)
}

maximizeEpsilon <- function(model) {
  stopifnot(!is.null(model$epsilonIndex))
  
  return (extremizeVariable(model$constraints, model$epsilonIndex, TRUE))
}

isModelConsistent <- function(model) {
  if (all(model$constraints$lhs[, model$epsilonIndex] == 0)) {
    # epsilon is not used, but still the model can be feasible
    ret <- extremizeVariable(model$constraints, model$epsilonIndex - 1, TRUE)
    
    return (ret$status == 0)
  } else {
    ret <- maximizeEpsilon(model)
  
    return (ret$status == 0 && ret$optimum >= model$minEpsilon)
  }
}

getRanksFromF <- function(model, values, accuracy = 1e-16) {
  nrVariables <- ncol(model$constraints$lhs)
  nrAlternatives <- nrow(model$perfToModelVariables)
  
  comprehensiveValue <- sapply(seq_len(nrAlternatives), function(i) {return (sum(ua(i, nrVariables, model$perfToModelVariables) * values) )} )
  
  ranks <- sapply(seq_len(nrAlternatives), function(i) {
    rank <- 1
    for (j in seq_len(nrAlternatives)) {
      if (i != j && comprehensiveValue[j] - comprehensiveValue[i] > accuracy) {
        rank <- rank + 1
      }
    }
    
    return (rank)
  } )
  
  return (ranks)
}

#' @export
toSolution <- function(model, values) {
  nrVariables <- ncol(model$constraints$lhs)
  nrAlternatives <- nrow(model$perfToModelVariables)
  nrCriteria <- ncol(model$perfToModelVariables)
  
  if (!is.null(model$epsilonIndex) && length(values) == nrVariables - 1) {
    model <- eliminateEpsilon(model)
    nrVariables <- ncol(model$constraints$lhs)
  }
  
  stopifnot(length(values) == nrVariables)
  
  #ranks
  
  ranks <- getRanksFromF(model, values)
  
  # epsilon
  
  epsilon <- NULL
  
  if (!is.null(model$epsilonIndex)) {
    epsilon <- values[model$epsilonIndex]
  } else {
    epsilon <- model$minEpsilon
  }
  
  # vf
  
  vf <- list()
  
  for (j in seq_len(nrCriteria)) {
    nrValues <- length(model$criterionValues[[j]])
    
    if (model$generalVF[j]) {
      x <- sapply(model$criterionValues[[j]], function(w) { w$value })
    } else {
      firstValue <- model$criterionValues[[j]][[1]]$value
      lastValue <- model$criterionValues[[j]][[length(model$criterionValues[[j]])]]$value
      intervalLength <- (lastValue - firstValue) / (model$chPoints[j] - 1)
      
      x <- c(firstValue,
             unlist(sapply(seq_len(model$chPoints[j] - 2), function(w) { firstValue + intervalLength * w })),
             lastValue)
    }
    
    y <- values[model$firstChPointVariableIndex[j] : (model$firstChPointVariableIndex[j] + model$chPoints[j] - 2)]
        
    if (model$criterionPreferenceDirection[j] == "g") {
      y <- c(0, y)
    } else {
      y <- c(y, 0)
    }
    
    vf[[j]] <- cbind(x, y)
  }
  
  # alternative values
  alternativeValues <- matrix(nrow=nrAlternatives, ncol=nrCriteria)
  
  for (i in seq_len(nrAlternatives)) {
    for (j in seq_len(nrCriteria)) {
      alternativeValues[i, j] <- 0
      
      for (k in seq_len(length(model$perfToModelVariables[[i, j]]))) {
        alternativeValues[i, j] <- alternativeValues[i, j] + values[model$perfToModelVariables[[i, j]][[k]][1]] * model$perfToModelVariables[[i, j]][[k]][2]
      }
    }
  }
  
  return (list(
    vf = vf,
    ranks = ranks,
    alternativeValues = alternativeValues,
    solution = values,
    epsilon = epsilon,
    generalVF = model$generalVF
    ))
}

#' @export
ranksToRanking <- function(ranks) {
  result <- c()
  altOrder <- order(ranks)
  prevRank <- -1
  
  for (alt in altOrder) {
    if (length(result) > 0) {
      if (prevRank == ranks[alt]) {
        result <- c(result, ",")
      } else {
        result <- c(result, "-")
      }
    }
    
    prevRank <- ranks[alt]
    
    result <- c(result, alt)
  }
  
  return (paste(result, collapse = ""))
}
