#' Bootstrap standard error estimates
#'
#' Function to obtain bootstrap standard error estimates for the parameter
#' estimates of the \code{\link{get_estimates}} function, under the generalized
#' linear model (GLM) or accelerated failure time (AFT) setting for the analysis
#' of a normally-distributed or censored time-to-event primary outcome.
#'
#' Under the GLM setting for the analysis of a normally-distributed primary
#' outcome Y, bootstrap standard error estimates are obtained for the estimates
#' of the parameters
#' \eqn{\alpha_0, \alpha_1, \alpha_2, \alpha_3, \sigma_1^2, \alpha_4, \alpha_{XY}, \sigma_2^2}{\alpha0, \alpha1, \alpha2, \alpha3, \sigma1^2, \alpha4, \alphaXY, \sigma2^2}
#' in the models
#' \deqn{Y = \alpha_0 + \alpha_1 \cdot K + \alpha_2 \cdot X + \alpha_3 \cdot L + \epsilon_1, \epsilon_1 \sim N(0,\sigma_1^2)}{Y = \alpha0 + \alpha1*K + \alpha2*X + \alpha3*L + \epsilon1, \epsilon1 ~ N(0,\sigma1^2)}
#' \deqn{Y^* = Y - \overline{Y} - \alpha_1 \cdot (K-\overline{K})}{Y* = Y - mean(Y) - \alpha1*(K-mean(K))}
#' \deqn{Y^* = \alpha_0 + \alpha_{XY} \cdot X + \epsilon_2, \epsilon_2 \sim N(0,\sigma_2^2),}{Y* = \alpha0 + \alphaXY*X + \epsilon2, \epsilon2 ~ N(0,\sigma2^2),}
#' accounting for the additional variability from the 2-stage approach.
#'
#' Under the AFT setting for the analysis of a censored time-to-event primary
#' outcome, bootstrap standard error estimates are similarly obtained of the
#' parameter estimates of
#' \eqn{\alpha_0, \alpha_1, \alpha_2, \alpha_3, \sigma_1, \alpha_4, \alpha_{XY}, \sigma_2^2}{\alpha0, \alpha1, \alpha2, \alpha3, \sigma1, \alpha4, \alphaXY, \sigma2^2}
#'
#' @param setting String with value \code{"GLM"} or \code{"AFT"} indicating
#'                whether standard error estimates are obtained for a
#'                normally-distributed (\code{"GLM"}) or censored time-to-event
#'                (\code{"AFT"}) primary outcome \code{Y}.
#' @param BS_rep Integer indicating the number of bootstrap samples that are drawn.
#' @param Y Numeric input vector for the primary outcome.
#' @param X Numeric input vector for the exposure variable.
#' @param K Numeric input vector for the intermediate outcome.
#' @param L Numeric input vector for the observed confounding factor.
#' @param C Numeric input vector for the censoring indicator under the AFT setting
#'          (must be coded 0 = censored, 1 = uncensored).
#'
#' @return Returns a vector with the bootstrap standard error estimates
#'         of the parameter estimates.
#'
#' @examples
#'
#' dat <- generate_data(setting = "GLM", n = 100)
#'
#' # For illustration use here only 100 bootstrap samples, recommended is using 1000
#' bootstrap_se(setting = "GLM", BS_rep = 100, Y = dat$Y, X = dat$X,
#'              K = dat$K, L = dat$L)
#'
#' @export
#'

bootstrap_se <- function(setting = "GLM", BS_rep = 1000, Y = NULL, X = NULL,
                         K = NULL, L = NULL, C = NULL) {
    if (is.null(setting)) {
        stop("setting has to be supplied.")
    }
    if (is.null(Y) | is.null(X) | is.null(K) | is.null(L)) {
        stop("Data of one or more variables are not supplied.")
    }
    if (setting == "AFT" & is.null(C)) {
        stop("C has to be supplied for the AFT setting.")
    }
    n <- length(Y)

    if (setting == "GLM") {
        alpha_0_SE_BS_help <- alpha_1_SE_BS_help <-
          alpha_2_SE_BS_help <- alpha_3_SE_BS_help <-
          sigma_1_sq_SE_BS_help <- alpha_4_SE_BS_help <-
          alpha_XY_SE_BS_help <- sigma_2_sq_SE_BS_help <- NULL
    }
    if (setting == "AFT") {
        alpha_0_SE_BS_help <- alpha_1_SE_BS_help <-
          alpha_2_SE_BS_help <- alpha_3_SE_BS_help <-
          sigma_1_SE_BS_help <- alpha_4_SE_BS_help <-
          alpha_XY_SE_BS_help <- sigma_2_sq_SE_BS_help <- NULL
    }
    for (rep in 1:BS_rep) {
        id <- sample(1:n, n, replace = TRUE)
        if (setting == "GLM") {
            data_id <- data.frame(X = X[id], L = L[id], K = K[id], Y = Y[id])
            estimates <- get_estimates(setting = setting, Y = data_id$Y, X = data_id$X,
                                       K = data_id$K, L = data_id$L)
        }
        if (setting == "AFT") {
            data_id <- data.frame(X = X[id], L = L[id], K = K[id], Y = Y[id],
                                  T = T[id], C = C[id])
            estimates <- get_estimates(setting = setting, Y = data_id$Y, X = data_id$X,
                                       K = data_id$K, L = data_id$L, C = data_id$C)
        }
        alpha_0_SE_BS_help[rep] <- estimates[names(estimates) == "alpha_0"]
        alpha_1_SE_BS_help[rep] <- estimates[names(estimates) == "alpha_1"]
        alpha_2_SE_BS_help[rep] <- estimates[names(estimates) == "alpha_2"]
        alpha_3_SE_BS_help[rep] <- estimates[names(estimates) == "alpha_3"]
        if (setting == "GLM") {
            sigma_1_sq_SE_BS_help[rep] <- estimates[names(estimates) == "sigma_1_sq"]
        }
        if (setting == "AFT") {
            sigma_1_SE_BS_help[rep] <- estimates[names(estimates) == "sigma_1"]
        }
        alpha_4_SE_BS_help[rep] <- estimates[names(estimates) == "alpha_4"]
        alpha_XY_SE_BS_help[rep] <- estimates[names(estimates) == "alpha_XY"]
        sigma_2_sq_SE_BS_help[rep] <- estimates[names(estimates) == "sigma_2_sq"]
    }
    if (setting == "GLM") {
        theta_bootstrap_se <- c(stats::sd(alpha_0_SE_BS_help,na.rm=T),
                                stats::sd(alpha_1_SE_BS_help,na.rm=T),
                                stats::sd(alpha_2_SE_BS_help,na.rm=T),
                                stats::sd(alpha_3_SE_BS_help,na.rm=T),
                                stats::sd(sigma_1_sq_SE_BS_help,na.rm=T),
                                stats::sd(alpha_4_SE_BS_help,na.rm=T),
                                stats::sd(alpha_XY_SE_BS_help,na.rm=T),
                                stats::sd(sigma_2_sq_SE_BS_help,na.rm=T))
        names(theta_bootstrap_se) <- c("alpha_0", "alpha_1", "alpha_2", "alpha_3",
                                       "sigma_1_sq", "alpha_4", "alpha_XY", "sigma_2_sq")
    }
    if (setting == "AFT") {
        theta_bootstrap_se <- c(stats::sd(alpha_0_SE_BS_help,na.rm=T),
                                stats::sd(alpha_1_SE_BS_help,na.rm=T),
                                stats::sd(alpha_2_SE_BS_help,na.rm=T),
                                stats::sd(alpha_3_SE_BS_help,na.rm=T),
                                stats::sd(sigma_1_SE_BS_help,na.rm=T),
                                stats::sd(alpha_4_SE_BS_help,na.rm=T),
                                stats::sd(alpha_XY_SE_BS_help,na.rm=T),
                                stats::sd(sigma_2_sq_SE_BS_help,na.rm=T))
        names(theta_bootstrap_se) <- c("alpha_0", "alpha_1", "alpha_2", "alpha_3",
                                       "sigma_1", "alpha_4", "alpha_XY", "sigma_2_sq")
    }
    return(theta_bootstrap_se)
}
