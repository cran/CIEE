#' Sandwich standard error estimates
#'
#' Function to obtain consistent and robust sandwich standard error estimates
#' based on estimating equations, for the parameter estimates of the
#' \code{\link{get_estimates}} function, under the GLM or AFT setting
#' for the analysis of a normally-distributed or censored time-to-event primary
#' outcome.
#'
#' Under the GLM setting for the analysis of a normally-distributed primary
#' outcome Y, robust sandwich standard error estimates are obtained for the
#' estimates of the parameters
#' \eqn{\alpha_0, \alpha_1, \alpha_2, \alpha_3, \sigma_1^2, \alpha_4, \alpha_{XY}, \sigma_2^2}{\alpha0, \alpha1, \alpha2, \alpha3, \sigma1^2, \alpha4, \alphaXY, \sigma2^2}
#' in the model
#' \deqn{Y = \alpha_0 + \alpha_1 \cdot K + \alpha_2 \cdot X + \alpha_3 \cdot L + \epsilon_1, \epsilon_1 \sim N(0,\sigma_1^2)}{Y = \alpha0 + \alpha1*K + \alpha2*X + \alpha3*L + \epsilon1, \epsilon1 ~ N(0,\sigma1^2)}
#' \deqn{Y^* = Y - \overline{Y} - \alpha_1 \cdot (K-\overline{K})}{Y* = Y - mean(Y) - \alpha1*(K-mean(K))}
#' \deqn{Y^* = \alpha_0 + \alpha_{XY} \cdot X + \epsilon_2, \epsilon_2 \sim N(0,\sigma_2^2)}{Y* = \alpha0 + \alphaXY*X + \epsilon2, \epsilon2 ~ N(0,\sigma2^2)}
#' by using the score and hessian matrices of the parameters.
#'
#' Under the AFT setting for the analysis of a censored time-to-event primary
#' outcome, robust sandwich standard error estimates are  similarly obtained of
#' the parameter estimates of
#' \eqn{\alpha_0, \alpha_1, \alpha_2, \alpha_3, \sigma_1, \alpha_4, \alpha_{XY}, \sigma_2^2}{\alpha0, \alpha1, \alpha2, \alpha3, \sigma1, \alpha4, \alphaXY, \sigma2^2}.
#' For more details and the underlying model, see the vignette.
#'
#' @param setting String with value \code{"GLM"} or \code{"AFT"} indicating
#'                whether standard error estimates are obtained for a
#'                normally-distributed (\code{"GLM"}) or censored time-to-event
#'                (\code{"AFT"}) primary outcome \code{Y}.
#' @param scores Score matrix of the parameters, which can be obtained using the
#'               \code{\link{scores}} function.
#' @param hessian Hessian matrix of the parameters, which can be obtained using the
#'               \code{\link{hessian}} function.
#'
#' @return Returns a vector with the CIEE sandwich standard error estimates
#'         of the parameter estimates.
#'
#' @examples
#'
#' # Generate data including Y, K, L, X under the GLM setting
#' dat <- generate_data(setting = "GLM")
#'
#' # Obtain estimating functions expressions
#' estfunct <- est_funct_expr(setting = "GLM")
#'
#' # Obtain point estimates of the parameters
#' estimates <- get_estimates(setting = "GLM", Y = dat$Y, X = dat$X,
#'                            K = dat$K, L = dat$L)
#'
#' # Obtain matrices with all first and second derivatives
#' derivobj <- deriv_obj(setting = "GLM", logL1 = estfunct$logL1,
#'                       logL2 = estfunct$logL2, Y = dat$Y, X = dat$X,
#'                       K = dat$K, L = dat$L, estimates = estimates)
#'
#' # Obtain score and hessian matrices
#' results_scores <- scores(derivobj)
#' results_hessian <- hessian(derivobj)
#'
#' # Obtain sandwich standard error estimates of the parameters
#' sandwich_se(scores = results_scores, hessian = results_hessian)
#'
#' @export
#'

sandwich_se <- function(setting = "GLM", scores = NULL, hessian = NULL) {
    if (is.null(scores) | is.null(hessian)) {
        stop("scores and hessian have to be supplied.")
    }
    n <- dim(scores)[1]

    ### A_n matrix ###
    A_n <- matrix(, nrow = 8, ncol = 8)
    for (A_n_i in 1:8) {
        for (A_n_j in 1:8) {
            A_n[A_n_i, A_n_j] <- sum(hessian[, A_n_i, A_n_j])
        }
    }
    A_n <- -(1/n) * A_n

    ### B_n matrix ###
    B_n <- matrix(, nrow = 8, ncol = 8)
    for (B_n_i in 1:8) {
        for (B_n_j in 1:8) {
            B_n[B_n_i, B_n_j] <- sum(scores[, B_n_i] * scores[, B_n_j])
        }
    }
    B_n <- (1/n) * B_n

    ### C_n matrix ###
    C_n <- solve(A_n) %*% B_n %*% t(solve(A_n))

    ### Variance estimates of coefficients ###
    theta_EE_se <- (1/n) * diag(C_n)
    theta_EE_se <- sqrt(theta_EE_se)
    if (setting == "GLM") {
        names(theta_EE_se) <- c("alpha_0", "alpha_1", "alpha_2", "alpha_3",
                                "sigma_1_sq", "alpha_4", "alpha_XY", "sigma_2_sq")
    }
    if (setting == "AFT") {
        names(theta_EE_se) <- c("alpha_0", "alpha_1", "alpha_2", "alpha_3",
                                "sigma_1", "alpha_4", "alpha_XY", "sigma_2_sq")
    }
    return(theta_EE_se)
}
