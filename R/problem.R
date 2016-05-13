#' @export
buildProblem <- function(perf, criteria, characteristicPoints, strictVF = TRUE,
                         strongPreference = NULL, weakPreference = NULL, indifference = NULL) {
  stopifnot(is.matrix(perf))
  stopifnot(is.vector(characteristicPoints))
  stopifnot(ncol(perf) == length(criteria))
  stopifnot(all(criteria %in% c('g', 'c')))
  stopifnot(ncol(perf) == length(characteristicPoints))
  stopifnot(is.logical(strictVF))
  stopifnot(all(characteristicPoints == 0 | characteristicPoints >= 2))
  
  if (is.null(weakPreference)) {
    weakPreference <- matrix(nrow=0, ncol=2)
  }
  
  stopifnot(ncol(weakPreference) == 2)
  stopifnot(all(weakPreference[,1] >= 1))
  stopifnot(all(weakPreference[,1] <= nrow(perf)))
  stopifnot(all(weakPreference[,2] >= 1))
  stopifnot(all(weakPreference[,2] <= nrow(perf)))
  
  if (is.null(strongPreference)) {
    strongPreference <- matrix(nrow=0, ncol=2)
  }
  
  stopifnot(ncol(strongPreference) == 2)
  stopifnot(all(strongPreference[,1] >= 1))
  stopifnot(all(strongPreference[,1] <= nrow(perf)))
  stopifnot(all(strongPreference[,2] >= 1))
  stopifnot(all(strongPreference[,2] <= nrow(perf)))
  
  if (is.null(indifference)) {
    indifference <- matrix(nrow=0, ncol=2)
  }
  
  stopifnot(ncol(indifference) == 2)
  stopifnot(all(indifference[,1] >= 1))
  stopifnot(all(indifference[,1] <= nrow(perf)))
  stopifnot(all(indifference[,2] >= 1))
  stopifnot(all(indifference[,2] <= nrow(perf)))
  
  return (list(perf = perf,
               strictVF = strictVF,
               criteria = criteria,
               characteristicPoints = characteristicPoints,
               strongPreference = strongPreference,
               weakPreference = weakPreference,
               indifference = indifference))
}

