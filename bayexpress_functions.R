# assume a flat prior
BE.default_u_1 <- 1
BE.default_u_2 <- 1

#' Calculate the logarithm of the Bayes' factor between 2 experiments
#'
#' @param N_1 Total number of reads for genes in experiment 1
#' @param n_1 Number of reads for a gene in experiment 1
#' @param N_2 Total number of reads for genes in experiment 2
#' @param n_2 Number of reads for a gene in experiment 2
#' @param u_1 Optional hyper-parameter (default is 1 for a flat prior)
#' @param u_2 Optional hyper-parameter (default is 1 for a flat prior)
#'
#' @return log10 Bayes' factor
#'
calc_log10_bayes_factor <- function(
    N_1,
    n_1,
    N_2,
    n_2,
    u_1 = BE.default_u_1,
    u_2 = BE.default_u_2) {

  # 'lbeta' returns the natural log of the beta function
  result <- lbeta(u_1 + n_1, u_2 + N_1 - n_1) +
    lbeta(u_1 + n_2, u_2 + N_2 - n_2) -
    lbeta(u_1 + n_1 + n_2, u_2 + N_1 - n_1 + N_2 - n_2)

  # convert natural log to log_10
  return(result / log(10))
}

#' Calculate the logarithm of the gene expression fold change between 2 experiments
#'
#' @param N_1 Total number of reads for genes in experiment 1
#' @param n_1 Number of reads for a gene in experiment 1
#' @param N_2 Total number of reads for genes in experiment 2
#' @param n_2 Number of reads for a gene in experiment 2
#' @param u_1 Optional hyper-parameter (default is 1 for a flat prior)
#' @param u_2 Optional hyper-parameter (default is 1 for a flat prior)
#'
#' @return log2 of the gene expression fold change
#'
calc_log2_fold_change <- function(
    N_1,
    n_1,
    N_2,
    n_2,
    u_1 = BE.default_u_1,
    u_2 = BE.default_u_2) {

  rate_1 <- exp(log(u_1 + n_1) - log(u_2 + N_1 - n_1))
  rate_2 <- exp(log(u_1 + n_2) - log(u_2 + N_2 - n_2))

  result <- log2(rate_2) - log2(rate_1)

  return(result)
}
