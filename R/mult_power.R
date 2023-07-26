#' Simple power calculator that can easily applied to constrained regression problems.
#' This function will determine the power of an alternative vs. a null hypothesis given constraints, a covariance estimate, and alpha level.

#' DISCALIMER: Like all functions in this package, this power calculation assumes normality in both the null and alternative distributions and uses the estimated covariance to construct both.
#' In the context of regression, these assumptions will hold a large proportion of the time, but it is important to check prior to using this function and the package in general.

#' @param G The hypothesized quadratic inequality constraint matrix such that: GB >= h
#' @param h The hypothesized quadratic inequality constraint vector such that: GB >= h
#' @param A The hypothesized quadratic equality constraint matrix such that: AB = b
#' @param b The hypothesized quadratic equality constraint vector such that: AB = b
#' @param cov: The estimated covariance matrix
#' @param true_mu: The mean you wish to test against the null constraints. (Test your estimated mean - power is another metric to see how badly the constraints are violated)
#' @param type: To do the IU test, input "IU". To do the EI test, input "EI"
#' @param alpha: The parameter that dictates the level of signifigance you wish to test the power at. It will default to .05 if nothing else is entered.
#' @return The power of true_mu vs. mu_hat (least favorable mean given the constraints)
#' @export
mult_power <- function (cov,true_mu,G,h,A,b,type,alpha=.05){
  draws1 <- 10000
  draws2 <- 1000
  R <- rbind(A,G)
  r <- c(b,h)
  ineq <- length(h)
  eq <- length(b)
  RXXR <- R %*% cov %*% t(R)
  mean <- matrix(rep(0,ineq+eq))
  X <- mvrnorm(draws1, mean, RXXR)
  iucdf <- c()
  eicdf <- c()
  for (j in 1:draws1) {
    x <- X[j,]
    P <- solve(RXXR)
    q <- (x %*% P)
    G <- diag(ineq+eq)
    h <- matrix(0, nrow = ineq+eq, ncol = 1)
    solution <- solve.QP(Dmat = P, dvec = q, Amat = t(G), bvec = h, meq = eq)
    optimal <- solution$solution
    x_lagrange <- as.vector(solution$Lagrangian)
    indices <- which(x_lagrange != 0)
    optimal[indices] <- 0
    IU_TS <- t(x - optimal) %*% P %*% (x - optimal)
    iucdf <- c(iucdf, IU_TS)
    EI_TS <- t(optimal) %*% P %*% (optimal)
    eicdf <- c(eicdf, EI_TS)
  }
  sorted_iucdf <- sort(iucdf)
  sorted_eicdf <- sort(eicdf)
  p <- 1-alpha
  iu_upper_percentile <- sorted_iucdf[as.integer(p * draws1)]
  ei_upper_percentile <- sorted_eicdf[as.integer(p * draws1)]
  if (type=="IU"){
    upper_percentile <- iu_upper_percentile
  }
  else{
    upper_percentile <- ei_upper_percentile
  }
  true_mean <- matrix((R %*% true_mu) - r)
  Y <- mvrnorm(draws2,true_mean,RXXR)
  frequency <- 0
  for (i in 1:draws2) {
    y <- Y[i,]
    P <- solve(RXXR)
    q <- (y %*% P)
    G <- diag(ineq+eq)
    h <- matrix(0, nrow = ineq+eq, ncol = 1)
    solution <- solve.QP(Dmat = P, dvec = q, Amat = t(G), bvec = h, meq = eq)
    optimal <- solution$solution
    x_lagrange <- as.vector(solution$Lagrangian)
    indices <- which(x_lagrange != 0)
    optimal[indices] <- 0
    if (type=="IU"){
      TS <- t(y - optimal) %*% P %*% (y - optimal)
    }
    else{
      TS <- t(optimal) %*% P %*% (optimal)
    }
    if (TS >= upper_percentile) {
      frequency <- frequency + 1
    }
  }
  result <- frequency/draws2
  if (type=="IU"){
    print(paste("IU Power:",result))
  }
  else{
    print(paste("EI Power:",result))
  }
  return(result)
}



