#' CIEE parameter point estimates
#'
#' Function to perform CIEE to obtain point estimates under the GLM or AFT
#' setting for the analysis of a normally-distributed or censored time-to-event
#' primary outcome.
#'
#' Under the GLM setting for the analysis of a normally-distributed primary
#' outcome Y, estimates of the parameters
#' \eqn{\alpha_0, \alpha_1, \alpha_2, \alpha_3, \sigma_1^2, \alpha_4, \alpha_{XY}, \sigma_2^2}{\alpha0, \alpha1, \alpha2, \alpha3, \sigma1^2, \alpha4, \alphaXY, \sigma2^2}
#' are obtained by constructing estimating equations for the models
#' \deqn{Y = \alpha_0 + \alpha_1 \cdot K + \alpha_2 \cdot X + \alpha_3 \cdot L + \epsilon_1, \epsilon_1 \sim N(0,\sigma_1^2)}{Y = \alpha0 + \alpha1*K + \alpha2*X + \alpha3*L + \epsilon1, \epsilon1 ~ N(0,\sigma1^2)}
#' \deqn{Y^* = Y - \overline{Y} - \alpha_1 \cdot (K-\overline{K})}{Y* = Y - mean(Y) - \alpha1*(K-mean(K))}
#' \deqn{Y^* = \alpha_0 + \alpha_{XY} \cdot X + \epsilon_2, \epsilon_2 \sim N(0,\sigma_2^2).}{Y* = \alpha0 + \alphaXY*X + \epsilon2, \epsilon2 ~ N(0,\sigma2^2).}
#' Under the AFT setting for the analysis of a censored time-to-event primary
#' outcome, estimates of the parameters
#' \eqn{\alpha_0, \alpha_1, \alpha_2, \alpha_3, \sigma_1, \alpha_4, \alpha_{XY}, \sigma_2^2}{\alpha0, \alpha1, \alpha2, \alpha3, \sigma1, \alpha4, \alphaXY, \sigma2^2}
#' are obtained by constructing
#' similar estimating equations based on a censored regression model and adding
#' an additional computation to estimate the true underlying survival times.
#' In addition to the parameter estimates, the mean of the estimated true
#' survival times is computed and returned in the output. For more details and
#' the underlying model, see the vignette.
#'
#' For both settings, the point estimates based on estimating equations equal
#' least squares (and maximum likelihood) estimates, and are obtained using
#' the \code{\link[stats]{lm}} and \code{\link[survival]{survreg}}
#' functions for computational purposes.
#'
#' @param setting String with value \code{"GLM"} or \code{"AFT"} indicating
#'                whether CIEE point estimates are obtained for a
#'                normally-distributed (\code{"GLM"}) or censored time-to-event
#'                (\code{"AFT"}) primary outcome \code{Y}.
#' @param Y Numeric input vector for the primary outcome.
#' @param X Numeric input vector for the exposure variable.
#' @param K Numeric input vector for the intermediate outcome.
#' @param L Numeric input vector for the observed confounding factor.
#' @param C Numeric input vector for the censoring indicator under the AFT setting
#'          (must be coded 0 = censored, 1 = uncensored).
#'
#' @return Returns a list with point estimates of the parameters. Under the
#'         AFT setting, the mean of the estimated true survival times is also
#'         computed and returned.
#'
#' @examples
#'
#' dat_GLM <- generate_data(setting = "GLM")
#' get_estimates(setting = "GLM", Y = dat_GLM$Y, X = dat_GLM$X, K = dat_GLM$K,
#'               L = dat_GLM$L)
#'
#' dat_AFT <- generate_data(setting = "AFT", a = 0.2, b = 4.75)
#' get_estimates(setting = "AFT", Y = dat_AFT$Y, X = dat_AFT$X, K = dat_AFT$K,
#'               L = dat_AFT$L, C = dat_AFT$C)
#'
#' @export
#'

get_estimates <- function(setting = "GLM", Y = NULL, X = NULL, K = NULL,
                          L = NULL, C = NULL) {
    if (!requireNamespace("survival", quietly = TRUE)) {
        stop("Pkg needed for this function to work. Please install it.", call. = FALSE)
    }
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
        data_help <- data.frame(Y = Y, X = X, K = K, L = L)
        data_help <- data_help[stats::complete.cases(data_help), ]
        ######### Stage 1 #########
        fit_stage_1 <- stats::lm(data_help$Y ~ data_help$K + data_help$X + data_help$L)

        alpha_0_out <- summary(fit_stage_1)$coefficients[1, 1]
        alpha_1_out <- summary(fit_stage_1)$coefficients[2, 1]
        alpha_2_out <- summary(fit_stage_1)$coefficients[3, 1]
        alpha_3_out <- summary(fit_stage_1)$coefficients[4, 1]
        sigma_1_sq_out <- (n - 4)/n * summary(fit_stage_1)$sigma^2

        ######### Stage 2 #########
        Y_tilde <- data_help$Y - mean(data_help$Y) -
                   alpha_1_out * (data_help$K - mean(data_help$K))
        fit_stage_2 <- stats::lm(Y_tilde ~ data_help$X)

        alpha_4_out <- summary(fit_stage_2)$coefficients[1, 1]
        alpha_XY_out <- summary(fit_stage_2)$coefficients[2, 1]
        sigma_2_sq_out <- (n - 2)/n * summary(fit_stage_2)$sigma^2

        point_estimates <- c(alpha_0_out, alpha_1_out, alpha_2_out, alpha_3_out,
                             sigma_1_sq_out, alpha_4_out, alpha_XY_out, sigma_2_sq_out)
        names(point_estimates) <- c("alpha_0", "alpha_1", "alpha_2", "alpha_3",
                                    "sigma_1_sq", "alpha_4", "alpha_XY", "sigma_2_sq")
    }
    if (setting == "AFT") {
        data_help <- data.frame(Y = Y, X = X, K = K, L = L, C = C)
        data_help <- data_help[stats::complete.cases(data_help), ]
        ######### Stage 1 #########
        fit_stage_1 <- survival::survreg(survival::Surv(data_help$Y, data_help$C) ~
                                         data_help$K + data_help$X + data_help$L,
                                         dist = "gaussian")

        alpha_0_out <- summary(fit_stage_1)$table[1, 1]
        alpha_1_out <- summary(fit_stage_1)$table[2, 1]
        alpha_2_out <- summary(fit_stage_1)$table[3, 1]
        alpha_3_out <- summary(fit_stage_1)$table[4, 1]
        sigma_1_out <- fit_stage_1$scale

        ######### Stage 2 #########
        mu <- fit_stage_1$linear.predictors
        Y_adj <- data_help$C * data_help$Y + (1 - data_help$C) * (mu +
                 (sigma_1_out * stats::dnorm((data_help$Y - mu)/sigma_1_out,
                 mean = 0, sd = 1)/(1 - stats::pnorm((data_help$Y - mu)/sigma_1_out,
                 mean = 0, sd = 1))))
        Y_tilde <- Y_adj - mean(Y_adj) -
                   alpha_1_out * (data_help$K - mean(data_help$K))
        fit_stage_2 <- stats::lm(Y_tilde ~ data_help$X)

        alpha_4_out <- summary(fit_stage_2)$coefficients[1, 1]
        alpha_XY_out <- summary(fit_stage_2)$coefficients[2, 1]
        sigma_2_sq_out <- (n - 2)/n * summary(fit_stage_2)$sigma^2

        point_estimates <- c(alpha_0_out, alpha_1_out, alpha_2_out, alpha_3_out,
                             sigma_1_out, alpha_4_out, alpha_XY_out, sigma_2_sq_out,
                             mean(Y_adj))
        names(point_estimates) <- c("alpha_0", "alpha_1", "alpha_2", "alpha_3",
                                    "sigma_1", "alpha_4", "alpha_XY", "sigma_2_sq",
                                    "y_adj_bar")
    }
    return(point_estimates)
}
