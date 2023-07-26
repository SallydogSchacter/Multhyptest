#' Performs the IU and EI hypothesis tests which are especially useful in the context of constrained regression.
#' IU: Shorthand for null hypothesis of inequality constraints vs. an unrestricted alternative.
#' EI: Shorthand for equality constraints null hypothesis vs. inequality constraints alternative (greater than the null).
#' The IU hypothesis test in this function assumes that all null constraints that you enter are inequality constraints.
#' The EI hypothesis test in this function assumes that all null constraints that you enter are equality constraints.
#' For a more complex analysis that can take non-diagonal constraint matrices and employ inequality and equality constraints at the same time, use complex_hyptest.

#' DISCLAIMER: Like all functions in this package, both tests assumes normality in both the null and alternative distributions and uses the estimated covariance to construct both.
#' In the context of regression, these assumptions will hold a large proportion of the time, but it is important to check prior to using this function and the package in general.

#' @param cov_matrix: estimated covariance matrix
#' @param mean_hat: estimated mean
#' @param G: Must be DIAGONAL The hypothesized quadratic constraint matrix such that: Gx >= h or Gx = h under the null (both tests are done here)
#' @param h: The hypothesized quadratic constraint vector such that: Gx >= h or Gx = h under the null (both tests are done here)
#' @param type: To do the IU test, input "IU". To do the EI test, input "EI"
#' @return Graphs the CDF of the EI and IU distributions and the p-value of the respective hypothesis test (remember this is always a one-sided test!)
#' @export
hyptest <- function(cov,mean_hat,R,r,type,input_draws=10000) {
  design <- cov
  size <- length(r)
  RXXR <- R %*% design %*% t(R)
  mu_hat_bar <- matrix((R %*% mean_hat) - r)
  quad <- function(l_bar){
    P <- (solve(RXXR))
    q <- matrix((t(l_bar) %*% P))
    G <- diag(size)
    h <- matrix(0, nrow = size, ncol = 1)
    if (type == "EI"){
      mq = 0
    }
    else{
      mq = size
    }
    solution <- solve.QP(P, q, t(G), h, meq = mq)
    x_optimal <- as.vector(solution$solution)
    x_lagrange <- as.vector(solution$Lagrangian)
    sol <- (list(x_optimal,x_lagrange))
    tilda_hat <- (sol[[1]])
    indices <- which(sol[[2]] != 0)
    tilda_hat[indices] <- 0
    return(matrix(tilda_hat))
  }
  tilda_hat <- quad(mu_hat_bar)
  iu_stat <- t(mu_hat_bar - tilda_hat) %*% (solve(RXXR)) %*% (mu_hat_bar - tilda_hat)
  ei_stat <- t(tilda_hat) %*% solve(RXXR) %*% tilda_hat
  if (type=="EI"){
    print(paste("EI stat:",(ei_stat)))
  }
  else{
    print(paste("IU stat:",(iu_stat)))
  }

  n_mean <- matrix(0, nrow = size, ncol = 1)
  empirical_cdf <- function(input_draws){
    ei_stats <- c()
    iu_stats <- c()
    tildas <- list()
    for (j in 1:input_draws){
      mean_hatted <- mvrnorm(1, mu = n_mean, Sigma = (RXXR))
      tilda <- quad(mean_hatted)
      tildas <- append(tildas, list(tilda))
      iu_stats <- c(iu_stats, t(mean_hatted - tilda) %*% solve(RXXR) %*% (mean_hatted - tilda))
      ei_stats <- c(ei_stats, t(tilda) %*% solve(RXXR) %*% tilda)
    }
    return(list(iu_stats, ei_stats, tildas))
  }

  op <- empirical_cdf(input_draws)
  iucdf <- op[1]
  eicdf <- op[2]
  ts <- op[3]

  write.csv(iucdf, "iucdf.csv", row.names = FALSE)
  write.csv(eicdf, "eicdf.csv", row.names = FALSE)

  ei_dist_cdf <- function(c){
    emp_ei <- unlist(eicdf)
    sorted_ecdf <- sort(emp_ei)
    upper_percentile <- sorted_ecdf[0.950 * input_draws]
    print(paste("EI 95% Rejection Boundary:", upper_percentile))
    eie <- ecdf(emp_ei)
    probability <- eie(c)
    return(probability)
  }
  iu_dist_cdf <- function(c){
    emp_iu <- unlist(iucdf)
    sorted_ecdf <- sort(emp_iu)
    upper_percentile <- sorted_ecdf[0.950 * input_draws]
    print(paste("IU 95% Rejection Boundary:", upper_percentile))
    iue <- ecdf(emp_iu)
    probability <- iue(c)
    return(probability)
  }
  if (type=="EI"){
    ei_p <- (ei_dist_cdf(ei_stat))
    print(paste("EI Nonparametrically Calculated P-Val:", 1-ei_p))
    ed <- read.csv("eicdf.csv")
    # Define the number of bins
    num_bins <- 200
    eds = ed[[1]]
    # Calculate the bin edges
    bin_edges <- seq(min(eds), max(eds), length.out = num_bins + 1)
    # Compute the histogram
    hist_result <- hist(eds, breaks = bin_edges, plot = FALSE)
    counts <- hist_result$counts
    breaks <- hist_result$breaks
    plot(breaks, c(counts, 0), type = "n", xlab = "Values", ylab = "Frequency", main = "EI Distribution", xlim = range(breaks), ylim = c(0, max(counts)), xaxs = "i", yaxs = "i")
    rect(breaks[-length(breaks)], 0, breaks[-1], counts, col = "turquoise", border = NA)
    box()
    abline(v = ei_stat, col = "red", lty = "dashed")
    png(file = "histogram.png", width = 800, height = 600, units = "px", res = 300)
    dev.off()
  }
  else{
    iu_p <- (iu_dist_cdf(iu_stat))
    print(paste("IU Nonparametrically Calculated P-Val:", 1-iu_p))
    IU <- read.csv("iucdf.csv")
    # Define the number of bins
    num_bins <- 200
    ius = IU[[1]]
    # Calculate the bin edges
    bin_edges <- seq(min(ius), max(ius), length.out = num_bins + 1)
    # Compute the histogram
    hist_result <- hist(ius, breaks = bin_edges, plot = FALSE)
    counts <- hist_result$counts
    breaks <- hist_result$breaks
    plot(breaks, c(counts, 0), type = "n", xlab = "Values", ylab = "Frequency", main = "IU Distribution", xlim = range(breaks), ylim = c(0, max(counts)), xaxs = "i", yaxs = "i")
    rect(breaks[-length(breaks)], 0, breaks[-1], counts, col = "tan", border = NA)
    box()
    abline(v = iu_stat, col = "red", lty = "dashed")
    png(file = "histogram.png", width = 800, height = 600, units = "px", res = 300)
    dev.off()
  }

  monte_carlo_weights <- function(input_draws){
    hats <- list()
    for (z in 1:input_draws){
      tilda <- ts[[1]][[z]]
      hats <- append(hats, sum(tilda > 0))
    }
    return(hats)
  }
  weights <- monte_carlo_weights(input_draws)
  ei_dist_monte <- function(c){
    total <- 0
    for (k in 0:(size)){
      df <- k
      if (df == 0){
        probability <- 0
      } else {
        probability <- 1 - pchisq(c, df)
      }
      count <- sum(unlist(weights) == k)/input_draws
      total <- total + count * probability
    }
    print(paste("EI Parametrically Calculated P-Val:", total))
  }

  iu_dist_monte <- function(c){
    total <- 0
    for (k in 0:(size)){
      df <- k
      if (df == 0){
        probability <- 0
      } else {
        probability <- 1 - pchisq(c, df)
      }
      count <- sum(unlist(weights) == (size-k))/input_draws
      total <- total + count * probability
    }
    print(paste("IU Parametrically Calculated P-Val:", total))
  }
  if (type=="EI"){
    ei_dist_monte(ei_stat)
  }
  else{
    iu_dist_monte(iu_stat)
  }
}
