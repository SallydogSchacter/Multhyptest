#' Performs a hypothesis test on your regression coefficients (Can also do hypothesis on non-regression related problems).
#' This will tell inform you if placing constraints on your Beta coefficients is reasonable given the covariance structure.
#' The complex_hyptest function can perform the basics equality vs. inequality (EI) and inequality vs. unrestricted (IU) tests
#' like hyptest but adds a level of complexity to the problem. Unlike hyptest, complex_hyptest tacks on
#' another set of equality restrictions (Ax=b) into the null hypothesis. In both the EI and
#' IU cases, the alternative associated with these new restrictions is unrestricted.
#' A small p-value will indicate that the constraints you want to impose on your regression are unreasonable and do not respect the covariance structure of your data
#' If all your constraints are equality or inequality, don't use this function; the hyptest function in this package will work better and faster.

#' DISCLAIMER: Like all functions in this package, this test assumes normality in both the null and alternative distributions and uses the estimated covariance to construct both.
#' In the context of regression, these assumptions will hold a large proportion of the time, but it is important to check prior to using this function and the package in general.

#' @param cov: estimated covariance matrix
#' @param mean_hat: estimated mean
#' @param G The hypothesized quadratic inequality constraint matrix such that: GB >= h
#' @param h The hypothesized quadratic inequality constraint vector such that: GB >= h
#' @param A The hypothesized quadratic equality constraint matrix such that: AB = b
#' @param b The hypothesized quadratic equality constraint vector such that: AB = b
#' @param type: To do the IU test, input "IU". To do the EI test, input "EI"
#' @return Graphs the distribution of the test statistic and shows you how likely your unrestricted Beta estimate is given the restrictions you entered. REMEMBER: THIS IS A ONE SIDED TEST!!!
#' @export

complex_hyptest <- function(cov,mean_hat,G,h,A,b,type,input_draws=10000) {
  R <- rbind(A,G)
  r <- c(b,h)
  ineq <- length(h)
  eq <- length(b)
  design <- cov
  RXXR <- R %*% design %*% t(R)
  lambda_hat_bar <- matrix((R %*% mean_hat) - r)
  quad <- function(l_bar){
    P <- (solve(RXXR))
    q <- matrix((t(l_bar) %*% P))
    G <- diag(ineq+eq)
    h <- matrix(0, nrow = ineq+eq, ncol = 1)
    solution <- solve.QP(P, q, t(G), h, meq = eq)
    x_optimal <- as.vector(solution$solution)
    x_lagrange <- as.vector(solution$Lagrangian)
    sol <- (list(x_optimal,x_lagrange))
    tilda_hat <- (sol[[1]])
    indices <- which(sol[[2]] != 0)
    tilda_hat[indices] <- 0
    return(tilda_hat)
  }
  lambda_hat_optimal <- quad(lambda_hat_bar)
  iu_stat <- t(lambda_hat_bar-lambda_hat_optimal) %*% solve(RXXR) %*% (lambda_hat_bar-lambda_hat_optimal)
  ei_stat <- t(lambda_hat_optimal) %*% solve(RXXR) %*% (lambda_hat_optimal)
  if (type=="EI"){
    print(paste("EI stat:",(ei_stat)))
  }
  else{
    print(paste("IU stat:",(iu_stat)))
  }
  null_mean <- matrix(0, nrow = ineq+eq, ncol = 1)
  empirical_cdf <- function(input_draws){
    ei_stats <- c()
    iu_stats <- c()
    tildas <- list()
    for (j in 1:input_draws){
      lambda_bar <- mvrnorm(n = 1, mu = null_mean, Sigma = (RXXR))
      tilda <- quad(lambda_bar)
      tildas <- append(tildas,list(tilda))
      iu_term <- t(lambda_bar-tilda) %*% solve(RXXR) %*% (lambda_bar-tilda)
      iu_stats <- c(iu_stats,iu_term)
      ei_term <- t(tilda) %*% solve(RXXR) %*% (tilda)
      ei_stats <- c(ei_stats,ei_term)
    }
    return(list(iu_stats,ei_stats,tildas))
  }

  op <- empirical_cdf(input_draws)
  c_iucdf <- op[[1]]
  c_eicdf <- op[[2]]
  ts <- op[[3]]

  write.csv(c_iucdf, "c_iucdf.csv", row.names = FALSE)
  write.csv(c_eicdf, "c_eicdf.csv", row.names = FALSE)
  dual_dist_cdf <- function(c,t){
    emp_dual <- unlist(t)
    sorted_ecdf <- sort(emp_dual)
    upper_percentile <- sorted_ecdf[0.950 * input_draws]
    print(paste("95% Rejection Boundary", upper_percentile))
    edual <- ecdf(emp_dual)
    probability <- edual(c)
    return(probability)
  }
  if (type=="EI"){
    ei_p <- (dual_dist_cdf(ei_stat,c_eicdf))
    print(paste("Non-Parametrically Calculated EI P-Val:", 1-ei_p))
    ed <- read.csv("c_eicdf.csv")
    st <- ei_stat
  }
  else{
    iu_p <- (dual_dist_cdf(iu_stat,c_iucdf))
    print(paste("Non-Parametrically Calculated IU P-Val:", 1-iu_p))
    ed <- read.csv("c_iucdf.csv")
    st <- iu_stat
  }


  # Define the number of bins
  num_bins <- 200
  eds = ed[[1]]
  # Calculate the bin edges
  bin_edges <- seq(min(eds), max(eds), length.out = num_bins + 1)
  # Compute the histogram
  hist_result <- hist(eds, breaks = bin_edges, plot = FALSE)
  counts <- hist_result$counts
  breaks <- hist_result$breaks
  plot(breaks, c(counts, 0), type = "n", xlab = "Values", ylab = "Frequency", main = "Empirical Distribution", xlim = range(breaks), ylim = c(0, max(counts)), xaxs = "i", yaxs = "i")
  rect(breaks[-length(breaks)], 0, breaks[-1], counts, col = "turquoise", border = NA)
  box()
  abline(v = st, col = "red", lty = "dashed")
  png(file = "histogram.png", width = 800, height = 600, units = "px", res = 300)
  dev.off()

  monte_carlo_weights <- function(input_draws){
    hats <- list()
    for (z in 1:input_draws){
      tilda <- ts[[z]]
      selected_tilda <- tilda[(length(tilda) - ineq + 1):length(tilda)]
      hats <- append(hats, sum(selected_tilda > 0))
    }
    return(hats)
  }
  counts <- monte_carlo_weights(input_draws)

  iu_dist_monte <- function(c) {
    total <- 0
    for (k in 0:ineq) {
      df <- k+eq
      if (df == 0) {
        probability <- 0
      }
      else {
        probability <- 1 - pchisq(c, df)
      }
      count <- sum(counts == ineq-k) / input_draws
      total <- total + count * probability
    }
    print(paste("Parametrically Calculated IU P-val:", total))
  }
  ei_dist_monte <- function(c) {
    total <- 0
    for (k in 0:ineq) {
      df <- k
      if (df == 0) {
        probability <- 0
      }
      else {
        probability <- 1 - pchisq(c, df)
      }
      count <- sum(counts == k) / input_draws
      total <- total + count * probability
    }
    print(paste("Parametrically Calculated EI P-val:", total))
  }
  if (type ==  "IU"){
    iu_dist_monte(iu_stat)
  }
  else{
    ei_dist_monte(ei_stat)
  }
}
