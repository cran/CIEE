#' CIEE: Causal inference based on estimating equations
#'
#' Functions to perform CIEE under the GLM or AFT setting:
#' \code{\link{ciee}} obtains point and standard error estimates of all parameter estimates,
#' and p-values for testing the absence of effects; \code{\link{ciee_loop}} performs
#' \code{\link{ciee}} in separate analyses of multiple exposure variables with the same
#' outcome measures and factors ond only returns point estimates, standard error
#' estimates and p-values for the exposure variables. Both functions can also compute
#' estimates and p-values from the two traditional regression methods and from the
#' structural equation modeling method.
#'
#' For the computation of CIEE, point estimates of the parameters are obtained
#' using the \code{\link{get_estimates}} function. Robust sandwich (recommended),
#' bootstrap, or naive standard error estimates of the parameter estimates are
#' obtained using the \code{\link{sandwich_se}}, \code{\link{bootstrap_se}}
#' or \code{\link{naive_se}} function. Large-sample Wald-type tests are performed
#' for testing the absence of effects, using either the robust sandwich or
#' bootstrap standard errors.
#'
#' Regarding the traditional regression methods, the multiple regression or
#' regression of residual approaches can be computed using the
#' \code{\link{mult_reg}} and \code{\link{res_reg}} functions. Finally, the
#' structural equation modeling approachcan be performed using the
#' \code{\link{sem_appl}} function.
#'
#' @param setting String with value \code{"GLM"} or \code{"AFT"} indicating
#'                whether a normally-distributed (\code{"GLM"}) or censored
#'                time-to-event (\code{"AFT"}) primary outcome \code{Y} is
#'                analyzed.
#' @param estimates String vector with possible values \code{"ee"}, \code{"mult_reg"},
#'                  \code{"res_reg"}, \code{"sem"} indicating which methods
#'                  are computed. \code{"ee"} computes CIEE, \code{"mult_reg"}
#'                  the traditional multiple regression method, \code{"res_reg"}
#'                  the traditional regression of residuals method, and
#'                  \code{"sem"} the structural equation modeling approach.
#'                  Multiple methods can be specified.
#' @param ee_se String with possible values \code{"sandwich"}, \code{"bootstrap"},
#'              \code{"naive"} indicating how the standard error estimates
#'              of the parameter estiamtes are computed for CIEE approach.
#'              \code{"sandwich"} computes the robust sandwich estimates (default,
#'              recommended), \code{"bootstrap"} the bootstrap estimates and
#'              \code{"naive"} the naive unadjusted standard error estimates
#'              (not recommended, only for comparison).
#'              One method has to be specified.
#' @param BS_rep Integer indicating the number of bootstrap samples that are
#'               drawn (recommended 1000) if bootstrap standard errors are computed.
#' @param Y Numeric input vector for the primary outcome.
#' @param X Numeric input vector for the exposure variable if the \code{\link{ciee}}
#'          function is used; or numeric input dataframe containing the exposure
#'          variables as columns if the \code{\link{ciee_loop}} function is used.
#' @param K Numeric input vector for the intermediate outcome.
#' @param L Numeric input vector for the observed confounding factor.
#' @param C Numeric input vector for the censoring indicator under the AFT setting
#'          (must be coded 0 = censored, 1 = uncensored).
#'
#' @return Object of class \code{ciee}, for which the summary function
#'         \code{\link{summary.ciee}} is implemented.
#'         \code{\link{ciee}} returns a list containing the point and standard error
#'         estimates of all parameters as well as p-values from hypothesis tests
#'         of the absence of effects, for each specified approach.
#'         \code{\link{ciee_loop}} returns a list containing the point and standard
#'         error estimates only of the exposure variables as well as p-values from
#'         hypothesis tests of the absence of effects, for each specified approach.
#'
#' @examples
#'
#' # Generate data under the GLM setting with default values
#' maf <- 0.2
#' n <- 100
#' dat <- generate_data(n = n, maf = maf)
#' datX <- data.frame(X = dat$X)
#' names(datX)[1] <- "X1"
#' # Add 9 more exposure variables names X2, ..., X10 to X
#' for (i in 2:10){
#'   X <- stats::rbinom(n, size = 2, prob = maf)
#'   datX$X <- X
#'   names(datX)[i] <- paste("X", i, sep="")
#' }
#'
#' # Perform analysis of one exposure variable using all four methods
#' ciee(Y = dat$Y, X = datX$X1, K = dat$K, L = dat$L)
#'
#' # Perform analysis of all exposure variables only for CIEE
#' ciee_loop(estimates = "ee", Y = dat$Y, X = datX, K = dat$K, L = dat$L)
#'
#' @export
#'

ciee <- function(setting = "GLM", estimates = c("ee", "mult_reg", "res_reg", "sem"),
                 ee_se = c("sandwich"), BS_rep = NULL, Y = NULL, X = NULL, K = NULL,
                 L = NULL, C = NULL) {
    if (is.null(setting)) {
        stop("setting has to be supplied.")
    }
    if (is.null(estimates)) {
        stop("At least one method has to be computed.")
    }
    if ((("ee" %in% estimates) & is.null(ee_se)) | (("ee" %in% estimates) &
          length(ee_se) > 1)) {
        stop("If the estimating equations approach is chosen, one approach has to be chosen for the computation of standard errors.")
    }
    if (("bootstrap" %in% estimates) & is.null(BS_rep)) {
        stop("For the computation of bootstrap standard errors, the number of bootstrap samples has to be chosen.")
    }
    if (is.null(Y) | is.null(X) | is.null(K) | is.null(L)) {
        stop("Data of one or more variables is not supplied.")
    }
    if (setting == "AFT" & is.null(C)) {
        stop("C has to be supplied for the AFT setting.")
    }
    if (setting == "AFT" & ("sem" %in% estimates)) {
        stop("The structural equations modeling approach is only implemented for the GLM setting.")
    }
    if (setting == "AFT" & ("res_reg" %in% estimates)) {
        stop("The regression of residuals approach is only implemented for the GLM setting.")
    }
    if (setting == "GLM") {
        data_help <- data.frame(Y = Y, X = X, K = K, L = L)
        data_help <- data_help[stats::complete.cases(data_help), ]
        if ("sem" %in% estimates) {
            results_sem <- sem_appl(Y = data_help$Y, X = data_help$X,
                                    K = data_help$K, L = data_help$L)
        }
        if ("mult_reg" %in% estimates) {
            results_mult_reg <- mult_reg(setting = setting, Y = data_help$Y,
                                         X = data_help$X, K = data_help$K,
                                         L = data_help$L)
        }
        if ("res_reg" %in% estimates) {
            results_res_reg <- res_reg(Y = data_help$Y, X = data_help$X,
                                       K = data_help$K, L = data_help$L)
        }
        if ("ee" %in% estimates) {
            point_estimates_ee <- get_estimates(setting = setting, Y = data_help$Y,
                                                X = data_help$X, K = data_help$K,
                                                L = data_help$L)
            if (ee_se == "sandwich") {
                # Obtain estimating functions expressions
                estfunct <- est_funct_expr(setting = "GLM")
                # Obtain matrices with all first and second derivatives
                derivobj <- deriv_obj(setting = setting, logL1 = estfunct$logL1,
                                      logL2 = estfunct$logL2, Y = data_help$Y,
                                      X = data_help$X, K = data_help$K,
                                      L = data_help$L,
                                      estimates = point_estimates_ee)
                # Obtain score and hessian matrices
                results_scores <- scores(derivobj)
                results_hessian <- hessian(derivobj)
                # Obtain sandwich standard error estimates of the parameters
                se_estimates_ee <- sandwich_se(setting = setting,
                                               scores = results_scores,
                                               hessian = results_hessian)
            }
            if (ee_se == "bootstrap") {
                se_estimates_ee <- bootstrap_se(setting = setting, BS_rep = BS_rep,
                                                Y = data_help$Y, X = data_help$X,
                                                K = data_help$K, L = data_help$L)
            }
            if (ee_se == "naive") {
                se_estimates_ee <- naive_se(setting = setting, Y = data_help$Y,
                                            X = data_help$X, K = data_help$K,
                                            L = data_help$L)
            }
            wald_test_stat_ee <- point_estimates_ee[1:8]/se_estimates_ee
            pvalues_ee <- 2 * stats::pnorm(-abs(wald_test_stat_ee))
            results_ee <- list(point_estimates = point_estimates_ee[1:8],
                               SE_estimates = se_estimates_ee,
                               wald_test_stat = wald_test_stat_ee,
                               pvalues = pvalues_ee)
        }
    }
    if (setting == "AFT") {
        data_help <- data.frame(Y = Y, X = X, K = K, L = L, C = C)
        data_help <- data_help[stats::complete.cases(data_help), ]
        if ("mult_reg" %in% estimates) {
            results_mult_reg <- mult_reg(setting = setting, Y = data_help$Y,
                                         X = data_help$X, K = data_help$K,
                                         L = data_help$L, C = data_help$C)
        }
        if ("ee" %in% estimates) {
            point_estimates_ee <- get_estimates(setting = setting, Y = data_help$Y,
                                                X = data_help$X, K = data_help$K,
                                                L = data_help$L, C = data_help$C)
            if (ee_se == "sandwich") {
                # Obtain estimating functions expressions
                estfunct <- est_funct_expr(setting = setting)
                # Obtain matrices with all first and second derivatives
                derivobj <- deriv_obj(setting = setting, logL1 = estfunct$logL1,
                                      logL2 = estfunct$logL2, Y = data_help$Y,
                                      X = data_help$X, K = data_help$K,
                                      L = data_help$L, C = data_help$C,
                                      estimates = point_estimates_ee)
                # Obtain score and hessian matrices
                results_scores <- scores(derivobj)
                results_hessian <- hessian(derivobj)
                # Obtain sandwich standard error estimates of the parameters
                se_estimates_ee <- sandwich_se(setting = setting,
                                               scores = results_scores,
                                               hessian = results_hessian)
            }
            if (ee_se == "bootstrap") {
                se_estimates_ee <- bootstrap_se(setting = setting, BS_rep = BS_rep,
                                                Y = data_help$Y, X = data_help$X,
                                                K = data_help$K, L = data_help$L,
                                                C = data_help$C)
            }
            if (ee_se == "naive") {
                se_estimates_ee <- naive_se(setting = setting, Y = data_help$Y,
                                            X = data_help$X, K = data_help$K,
                                            L = data_help$L, C = data_help$C)
            }
            wald_test_stat_ee <- point_estimates_ee[1:8]/se_estimates_ee
            pvalues_ee <- 2 * stats::pnorm(-abs(wald_test_stat_ee))
            results_ee <- list(point_estimates = point_estimates_ee[1:8],
                               SE_estimates = se_estimates_ee,
                               wald_test_stat = wald_test_stat_ee,
                               pvalues = pvalues_ee)
        }
    }
    output <- list()
    if ("ee" %in% estimates) {
        output$results_ee <- results_ee
    }
    if ("mult_reg" %in% estimates) {
        output$results_mult_reg <- results_mult_reg
    }
    if ("res_reg" %in% estimates) {
        output$results_res_reg <- results_res_reg
    }
    if ("sem" %in% estimates) {
        output$results_sem <- results_sem
    }
    class(output) <- "ciee"
    return(output)
}

#' @rdname ciee
#' @export

ciee_loop <- function(setting = "GLM", estimates = c("ee", "mult_reg", "res_reg",
                      "sem"), ee_se = c("sandwich"), BS_rep = NULL, Y = NULL,
                      X = NULL, K = NULL, L = NULL, C = NULL) {
    if (is.null(setting)) {
        stop("setting has to be supplied.")
    }
    if (is.null(estimates)) {
        stop("At least one method has to be computed.")
    }
    if ((("ee" %in% estimates) & is.null(ee_se)) | (("ee" %in% estimates)
         & length(ee_se) > 1)) {
        stop("If the estimating equations approach is chosen, one approach has to be chosen for the computation of standard errors.")
    }
    if (("bootstrap" %in% estimates) & is.null(BS_rep)) {
        stop("For the computation of bootstrap standard errors, the number of bootstrap samples has to be chosen.")
    }
    if (is.null(Y) | is.null(X) | is.null(K) | is.null(L)) {
        stop("Data of one or more variables is not supplied.")
    }
    if (setting == "AFT" & is.null(C)) {
        stop("C has to be supplied for the AFT setting.")
    }
    if (setting == "AFT" & ("sem" %in% estimates)) {
        stop("The structural equations modeling approach is only implemented for the GLM setting.")
    }
    if (setting == "AFT" & ("res_reg" %in% estimates)) {
        stop("The regression of residuals approach is only implemented for the GLM setting.")
    }
    if ("sem" %in% estimates) {
        results_sem <- list(point_estimates = NULL, SE_estimates = NULL,
                            pvalues = NULL)
    }
    if ("mult_reg" %in% estimates) {
        results_mult_reg <- list(point_estimates = NULL, SE_estimates = NULL,
                                 pvalues = NULL)
    }
    if ("res_reg" %in% estimates) {
        results_res_reg <- list(point_estimates = NULL, SE_estimates = NULL,
                            pvalues = NULL)
    }
    if ("ee" %in% estimates) {
        results_ee <- list(point_estimates = NULL, SE_estimates = NULL,
                           wald_test_stat = NULL, pvalues = NULL)
    }
    k <- dim(X)[2]
    for (i in 1:k) {
        Xi <- X[, i]
        if (setting == "GLM") {
            data_help <- data.frame(Y = Y, X = Xi, K = K, L = L)
            data_help <- data_help[stats::complete.cases(data_help), ]

            if ("sem" %in% estimates) {
                sem_help <- sem_appl(Y = data_help$Y, X = data_help$X,
                                     K = data_help$K, L = data_help$L)
                results_sem$point_estimates[i] <- sem_help$point_estimates[5]
                results_sem$SE_estimates[i] <- sem_help$SE_estimates[5]
                results_sem$pvalues[i] <- sem_help$pvalues[5]
                names(results_sem$point_estimates)[i] <-
                  names(results_sem$SE_estimates)[i] <-
                  names(results_sem$pvalues)[i] <- names(X)[i]
            }
            if ("mult_reg" %in% estimates) {
                mult_reg_help <- mult_reg(setting = setting, Y = data_help$Y,
                                          X = data_help$X, K = data_help$K,
                                          L = data_help$L)
                results_mult_reg$point_estimates[i] <- mult_reg_help$point_estimates[3]
                results_mult_reg$SE_estimates[i] <- mult_reg_help$SE_estimates[3]
                results_mult_reg$pvalues[i] <- mult_reg_help$pvalues[3]
                names(results_mult_reg$point_estimates)[i] <-
                  names(results_mult_reg$SE_estimates)[i] <-
                  names(results_mult_reg$pvalues)[i] <- names(X)[i]
            }
            if ("res_reg" %in% estimates) {
                res_reg_help <- res_reg(Y = data_help$Y, X = data_help$X,
                                        K = data_help$K, L = data_help$L)
                results_res_reg$point_estimates[i] <- res_reg_help$point_estimates[5]
                results_res_reg$SE_estimates[i] <- res_reg_help$SE_estimates[5]
                results_res_reg$pvalues[i] <- res_reg_help$pvalues[5]
                names(results_res_reg$point_estimates)[i] <-
                  names(results_res_reg$SE_estimates)[i] <-
                  names(results_res_reg$pvalues)[i] <- names(X)[i]
            }
            if ("ee" %in% estimates) {
                point_estimates_ee <- get_estimates(setting = setting,
                                                    Y = data_help$Y,
                                                    X = data_help$X,
                                                    K = data_help$K,
                                                    L = data_help$L)
                if (ee_se == "sandwich") {
                    # Obtain estimating functions expressions
                    estfunct <- est_funct_expr(setting = "GLM")
                    # Obtain matrices with all first and second derivatives
                    derivobj <- deriv_obj(setting = setting,
                                          logL1 = estfunct$logL1,
                                          logL2 = estfunct$logL2,
                                          Y = data_help$Y,
                                          X = data_help$X,
                                          K = data_help$K,
                                          L = data_help$L,
                                          estimates = point_estimates_ee)
                    # Obtain score and hessian matrices
                    results_scores <- scores(derivobj)
                    results_hessian <- hessian(derivobj)
                    # Obtain sandwich standard error estimates of the parameters
                    se_estimates_ee <- sandwich_se(setting = setting,
                                                   scores = results_scores,
                                                   hessian = results_hessian)
                }
                if (ee_se == "bootstrap") {
                    se_estimates_ee <- bootstrap_se(setting = setting,
                                                    BS_rep = BS_rep,
                                                    Y = data_help$Y,
                                                    X = data_help$X,
                                                    K = data_help$K,
                                                    L = data_help$L)
                }
                if (ee_se == "naive") {
                    se_estimates_ee <- naive_se(setting = setting,
                                                Y = data_help$Y,
                                                X = data_help$X,
                                                K = data_help$K,
                                                L = data_help$L)
                }
                results_ee$point_estimates[i] <- point_estimates_ee[7]
                results_ee$SE_estimates[i] <- se_estimates_ee[7]
                results_ee$wald_test_stat[i] <- point_estimates_ee[7]/
                                                  se_estimates_ee[7]
                results_ee$pvalues[i] <- 2 * stats::pnorm(-abs(point_estimates_ee[7]/
                                                        se_estimates_ee[7]))
                names(results_ee$point_estimates)[i] <-
                  names(results_ee$SE_estimates)[i] <-
                  names(results_ee$wald_test_stat)[i] <-
                  names(results_ee$pvalues)[i] <- names(X)[i]
            }
        }
        if (setting == "AFT") {
            data_help <- data.frame(Y = Y, X = Xi, K = K, L = L, C = C)
            data_help <- data_help[stats::complete.cases(data_help), ]
            if ("mult_reg" %in% estimates) {
                mult_reg_help <- mult_reg(setting = setting, Y = data_help$Y,
                                          X = data_help$X, K = data_help$K,
                                          L = data_help$L, C = data_help$C)
                results_mult_reg$point_estimates[i] <- mult_reg_help$point_estimates[3]
                results_mult_reg$SE_estimates[i] <- mult_reg_help$SE_estimates[3]
                results_mult_reg$pvalues[i] <- mult_reg_help$pvalues[3]
                names(results_mult_reg$point_estimates)[i] <-
                  names(results_mult_reg$SE_estimates)[i] <-
                  names(results_mult_reg$pvalues)[i] <- names(X)[i]
            }
            if ("ee" %in% estimates) {
                point_estimates_ee <- get_estimates(setting = setting,
                                                    Y = data_help$Y,
                                                    X = data_help$X,
                                                    K = data_help$K,
                                                    L = data_help$L,
                                                    C = data_help$C)
                if (ee_se == "sandwich") {
                    # Obtain estimating functions expressions
                    estfunct <- est_funct_expr(setting = setting)
                    # Obtain matrices with all first and second derivatives
                    derivobj <- deriv_obj(setting = setting,
                                          logL1 = estfunct$logL1,
                                          logL2 = estfunct$logL2,
                                          Y = data_help$Y,
                                          X = data_help$X,
                                          K = data_help$K,
                                          L = data_help$L,
                                          C = data_help$C,
                                          estimates = point_estimates_ee)
                    # Obtain score and hessian matrices
                    results_scores <- scores(derivobj)
                    results_hessian <- hessian(derivobj)
                    # Obtain sandwich standard error estimates of the parameters
                    se_estimates_ee <- sandwich_se(setting = setting,
                                                   scores = results_scores,
                                                   hessian = results_hessian)
                }
                if (ee_se == "bootstrap") {
                    se_estimates_ee <- bootstrap_se(setting = setting,
                                                    BS_rep = BS_rep,
                                                    Y = data_help$Y,
                                                    X = data_help$X,
                                                    K = data_help$K,
                                                    L = data_help$L,
                                                    C = data_help$C)
                }
                if (ee_se == "naive") {
                    se_estimates_ee <- naive_se(setting = setting,
                                                Y = data_help$Y,
                                                X = data_help$X,
                                                K = data_help$K,
                                                L = data_help$L,
                                                C = data_help$C)
                }
                results_ee$point_estimates[i] <- point_estimates_ee[7]
                results_ee$SE_estimates[i] <- se_estimates_ee[7]
                results_ee$wald_test_stat[i] <- point_estimates_ee[7]/
                                                  se_estimates_ee[7]
                results_ee$pvalues[i] <- 2 * stats::pnorm(-abs(point_estimates_ee[7]/
                                                        se_estimates_ee[7]))
                names(results_ee$point_estimates)[i] <-
                  names(results_ee$SE_estimates)[i] <-
                  names(results_ee$wald_test_stat)[i] <-
                  names(results_ee$pvalues)[i] <- names(X)[i]
            }
        }
    }
    output <- list()
    if ("ee" %in% estimates) {
        output$results_ee <- results_ee
    }
    if ("mult_reg" %in% estimates) {
        output$results_mult_reg <- results_mult_reg
    }
    if ("res_reg" %in% estimates) {
        output$results_res_reg <- results_res_reg
    }
    if ("sem" %in% estimates) {
        output$results_sem <- results_sem
    }
    class(output) <- "ciee"
    return(output)
}
