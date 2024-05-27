#'@rdname FHN_Splitting_Cpp
#'@title FHN_Splitting_Cpp
#'@description Simulate a trajectory from the stochastic FHN model
#' using Strang splitting
#'@return path of the FHN model
#'@export
FHN_Splitting_Cpp <- function(grid, h, startv, dm, cm, eps, beta){
  return(SplittingFHN::FHN_Splitting_Cpp_(grid, h, startv, dm, cm, eps, beta))
}
