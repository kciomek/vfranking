#' @import hitandrun
#' @export
sampleParameters <- function(model, numberOfSamples = 1000) {
  stopifnot(numberOfSamples > 0)
  stopifnot(all(model$constraints$types == "C"))
  
  if (!isModelConsistent(model)) {
    stop("Model infeasible.")
  }
  
  nrAlternatives <- nrow(model$perfToModelVariables)
  nrCriteria <- ncol(model$perfToModelVariables)
  
  model <- eliminateEpsilon(model)
  
  constraints <- model$constraints
  constraints$dir[which(constraints$dir == "==")] <- "="
  geq <- which(constraints$dir == ">=")
  
  for (i in geq) {
    constraints$rhs[i] <- -1 * constraints$rhs[i]
    constraints$lhs[i, ] <- -1 * constraints$lhs[i, ]
  }
  
  constraints$dir[geq] <- "<="
  names(constraints)[1] <- "constr"
  constraints[[4]] <- NULL
  
  return (hitandrun(constraints, n.samples = numberOfSamples, thin.fn = function(n) { n^3 }))
}


#' @export
pwi <- function(model, samples) {
  stopifnot(nrow(samples) > 0)
  
  nrAlternatives <- nrow(model$perfToModelVariables)
  nrCriteria <- ncol(model$perfToModelVariables)
  
  model <- eliminateEpsilon(model)
  
  result <- matrix(data = 0, nrow = nrAlternatives, ncol = nrAlternatives)
  
  for (i in seq_len(nrow(samples))) {
    ranks <- getRanksFromF(model, samples[i, ])
    
    for (i in seq_len(nrAlternatives)) {
      for (j in seq_len(nrAlternatives)) {
        if (ranks[i] < ranks[j]) {
          result[i, j] <- result[i, j] + 1
        }
      }
    }
  }
  
  result <- result / nrow(samples)
  
  return (result)
}


#' @export
rai <- function(model, samples) {
  stopifnot(nrow(samples) > 0)
  
  nrAlternatives <- nrow(model$perfToModelVariables)
  nrCriteria <- ncol(model$perfToModelVariables)
  
  model <- eliminateEpsilon(model)
  
  result <- matrix(data = 0, nrow = nrAlternatives, ncol = nrAlternatives)
  
  for (i in seq_len(nrow(samples))) {
    ranks <- getRanksFromF(model, samples[i, ])
    
    for (i in seq_len(nrAlternatives)) {
      result[i, ranks[i]] <- result[i, ranks[i]] + 1
    }
  }
  
  result <- result / nrow(samples)
  
  return (result)
}

