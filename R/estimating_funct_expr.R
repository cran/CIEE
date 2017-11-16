#' Estimating functions.
#'
#' Function to compute \code{logL1} and \code{logL2} under the GLM and AFT setting
#' for the analysis of a normally-distributed and of a censored time-to-event
#' primary outcome. \code{logL1} and \code{logL2} are functions which underlie
#' the estimating functions of CIEE for the derivation of point estimates and
#' standard error estimates. \code{\link{est_funct_expr}} computes their
#' expression, which is then further used in the functions \code{\link{deriv_obj}},
#' \code{\link{ciee}} and \code{\link{ciee_loop}}.
#'
#' Under the GLM setting for the analysis of a normally-distributed primary
#' outcome \code{Y}, the goal is to obtain estimates for the pararameters
#' \eqn{\alpha_0, \alpha_1, \alpha_2, \alpha_3, \sigma_1^2, \alpha_4, \alpha_{XY}, \sigma_2^2}{\alpha0, \alpha1, \alpha2, \alpha3, \sigma1^2, \alpha4, \alphaXY, \sigma2^2}
#' under the model
#' \deqn{Y = \alpha_0 + \alpha_1 \cdot K + \alpha_2 \cdot X + \alpha_3 \cdot L + \epsilon_1, \epsilon_1 \sim N(0,\sigma_1^2)}{Y = \alpha0 + \alpha1*K + \alpha2*X + \alpha3*L + \epsilon1, \epsilon1 ~ N(0,\sigma1^2)}
#' \deqn{Y^* = Y - \overline{Y} - \alpha_1 \cdot (K-\overline{K})}{Y* = Y - mean(Y) - \alpha1*(K-mean(K))}
#' \deqn{Y^* = \alpha_0 + \alpha_{XY} \cdot X + \epsilon_2, \epsilon_2 \sim N(0,\sigma_2^2)}{Y* = \alpha0 + \alphaXY*X + \epsilon2, \epsilon2 ~ N(0,\sigma2^2).}
#' \code{logL1} underlies the estimating functions for the derivation of the
#' first 5 parameters
#' \eqn{\alpha_0, \alpha_1, \alpha_2, \alpha_3, \sigma_1^2}{\alpha0, \alpha1, \alpha2, \alpha3, \sigma1^2}
#' and
#' \code{logL2} underlies the estimating functions for the derivation of the
#' last 3 parameters
#' \eqn{\alpha_4, \alpha_{XY}, \sigma_2^2}{\alpha4, \alphaXY, \sigma2^2}.
#'
#' Under the AFT setting for the analysis of a censored time-to-event primary
#' outcome \code{Y}, the goal is to obtain estimates of the parameters
#' \eqn{\alpha_0, \alpha_1, \alpha_2, \alpha_3, \sigma_1, \alpha_4, \alpha_{XY}, \sigma_2^2}{\alpha0, \alpha1, \alpha2, \alpha3, \sigma1, \alpha4, \alphaXY, \sigma2^2}.
#' Here, \code{logL1} similarly underlies the estimating functions
#' for the derivation of the first 5 parameters and \code{logL2} underlies the
#' estimating functions for the derivation of the last 3 parameters.
#'
#' \code{logL1}, \code{logL2} equal the log-likelihood functions (\code{logL2}
#' given that \eqn{\alpha_1}{\alpha1} is known). For more details and the underlying model,
#' see the vignette.
#'
#' @param setting String with value \code{"GLM"} or \code{"AFT"} indicating whether
#'                the expression of \code{logL1} and \code{logL2} is computed
#'                under the GLM or AFT setting.
#' @return Returns a list containing the expression of the functions \code{logL1}
#'         and \code{logL2}.
#'
#' @examples
#'
#' est_funct_expr(setting = "GLM")
#' est_funct_expr(setting = "AFT")
#'
#' @export
#'

est_funct_expr <- function(setting = "GLM") {
    if (is.null(setting)) {
        stop("setting has to be supplied.")
    }
    if (setting == "GLM") {
        logL1 <- expression(log((1/sqrt(sigma1sq)) * dnorm((y_i - alpha0 -
                            alpha1 * k_i - alpha2 * x_i - alpha3 * l_i)/
                            sqrt(sigma1sq), mean = 0, sd = 1)))
        logL2 <- expression(log((1/sqrt(sigma2sq)) * dnorm((y_i - y_bar -
                            alpha1 * (k_i - k_bar) - alpha4 - alphaXY * x_i)/
                            sqrt(sigma2sq), mean = 0, sd = 1)))
    }
    if (setting == "AFT") {
        logL1 <- expression(-c_i * log(sigma1) + c_i * log(dnorm((y_i -
                            alpha0 - alpha1 * k_i - alpha2 * x_i -
                            alpha3 * l_i)/sigma1, mean = 0, sd = 1)) +
                            (1 - c_i) * log(1 - pnorm((y_i - alpha0 -
                            alpha1 * k_i - alpha2 * x_i - alpha3 * l_i)/
                            sigma1, mean = 0, sd = 1)))
        logL2 <- expression(log((1/sqrt(sigma2sq)) * dnorm(((c_i * y_i +
                            (1 - c_i) * ((alpha0 + alpha1 * k_i +
                            alpha2 * x_i + alpha3 * l_i) + (sigma1 *
                            dnorm((y_i - alpha0 - alpha1 * k_i -
                            alpha2 * x_i - alpha3 * l_i)/sigma1, mean = 0,
                            sd = 1)/(1 - pnorm((y_i - alpha0 - alpha1 * k_i -
                            alpha2 * x_i - alpha3 * l_i)/sigma1, mean = 0,
                            sd = 1))))) - y_adj_bar - alpha1 * (k_i - k_bar) -
                            alpha4 - alphaXY * x_i)/sqrt(sigma2sq),
                            mean = 0, sd = 1)))
    }
    return(list(logL1 = logL1, logL2 = logL2))
}
