#' Summary function.
#'
#' Summary function for the \code{\link{ciee}} and \code{\link{ciee_loop}}
#' functions.
#'
#' @param object \code{ciee} object (output of the \code{\link{ciee}}
#'               or \code{\link{ciee_loop}} function).
#' @param ... Additional arguments affecting the summary produced.
#'
#' @return Formatted data frames of the results of all computed methods.
#'
#' @examples
#'
#' maf <- 0.2
#' n <- 1000
#' dat <- generate_data(n = n, maf = maf)
#' datX <- data.frame(X = dat$X)
#' names(datX)[1] <- "X1"
#' for (i in 2:10){
#'   X <- stats::rbinom(n, size = 2, prob = maf)
#'   datX$X <- X
#'   names(datX)[i] <- paste("X", i, sep="")
#' }
#'
#' results1 <- ciee(Y = dat$Y, X = datX$X1, K = dat$K, L = dat$L)
#' summary(results1)
#'
#' results2 <- ciee_loop(Y = dat$Y, X = datX, K = dat$K, L = dat$L)
#' summary(results2)
#'
#' @export
#'

summary.ciee <- function(object = NULL, ...) {
    if (is.null(object)) {
        stop("ciee output has to be supplied.")
    }
    res_out <- NULL
    if ("results_ee" %in% names(object)) {
        res_ee_out <- data.frame(point_estimates = object$results_ee$point_estimates,
                                 SE_estimates = object$results_ee$SE_estimates,
                                 wald_test_stat = object$results_ee$wald_test_stat,
                                 pvalues = object$results_ee$pvalues)
        rownames(res_ee_out) <- paste("CIEE", rownames(res_ee_out), sep = "_")
        print(paste("Results based on estimating equations."))
        print(res_ee_out)
        res_out <- res_ee_out[,c(1,2,4)]
    }
    if ("results_mult_reg" %in% names(object)) {
        res_mr_out <- data.frame(point_estimates = object$results_mult_reg$point_estimates,
                                 SE_estimates = object$results_mult_reg$SE_estimates,
                                 pvalues = object$results_mult_reg$pvalues)
        rownames(res_mr_out) <- paste("MR", rownames(res_mr_out), sep = "_")
        print(paste("Results based on traditional multiple regression."))
        print(res_mr_out)
        res_out <- rbind(res_out, res_mr_out)
    }
    if ("results_res_reg" %in% names(object)) {
        res_rr_out <- data.frame(point_estimates = object$results_res_reg$point_estimates,
                                 SE_estimates = object$results_res_reg$SE_estimates,
                                 pvalues = object$results_res_reg$pvalues)
        rownames(res_rr_out) <- paste("RR", rownames(res_rr_out), sep = "_")
        print(paste("Results based on traditional regression of residuals."))
        print(res_rr_out)
        res_out <- rbind(res_out, res_rr_out)
    }
    if ("results_sem" %in% names(object)) {
        res_sem_out <- data.frame(point_estimates = object$results_sem$point_estimates,
                                  SE_estimates = object$results_sem$SE_estimates,
                                  pvalues = object$results_sem$pvalues)
        rownames(res_sem_out) <- paste("SEM", rownames(res_sem_out), sep = "_")
        print(paste("Results based on structural equation modeling."))
        print(res_sem_out)
        res_out <- rbind(res_out, res_sem_out)
    }
    invisible(res_out)
}
