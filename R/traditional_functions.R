#' Traditional regression approaches.
#'
#' Functions to fit traditional regression approaches for a quantitative
#' normally-distributed primary outcome (\code{setting} = \code{"GLM"})
#' and a censoredtime-to-event primary outcome (\code{setting} = \code{"AFT"}).
#' \code{\link{mult_reg}} fits the multiple regression approach and
#' \code{\link{res_reg}} computes the regression of residuals approach.
#'
#' In more detail, for a quantitative normally-distributed primary outcome
#' \code{Y}, \code{\link{mult_reg}} fits the model
#' \deqn{Y = \alpha_0 + \alpha_1 \cdot K + \alpha_{XY} \cdot X + \alpha_2 \cdot L + \epsilon}{Y = \alpha0 + \alpha1*K + \alphaXY*X + \alpha2*L + \epsilon}
#' and obtains point and standard error estimates for the parameters
#' \eqn{\alpha_0, \alpha_1, \alpha_{XY}, \alpha_2}{\alpha0, \alpha1, \alphaXY, \alpha2}.
#' \code{\link{res_reg}} obtains point and standard
#' error estimates for the parameters
#' \eqn{\alpha_0, \alpha_1, \alpha_2, \alpha_3, \alpha_{XY}}{\alpha0, \alpha1, \alpha2, \alpha3, \alphaXY}
#' by fitting the models
#' \deqn{Y = \alpha_0 + \alpha_1 \cdot K + \alpha_2 \cdot L + \epsilon_1}{Y = \alpha0 + \alpha1*K + \alpha2*L + \epsilon1,}
#' \deqn{\widehat{\epsilon}_1 = \alpha_3 + \alpha_{XY} \cdot X + \epsilon_2}{hat(\epsilon1) = \alpha3 + \alphaXY*X + \epsilon2.}
#' Both functions use the \code{\link[stats]{lm}} function and also report the
#' provided p-values from t-tests that each parameter equals 0.
#' For the analysis of a censored time-to-event primary outcome \code{Y},
#' only the multiple regression approach is implemented. Here,
#' \code{\link{mult_reg}} fits the according censored regression model to obtain
#' coefficient and standard error estimates as well as p-values from large-sample
#' Wald-type tests by using the \code{\link[survival]{survreg}} function.
#' See the vignette for more details.
#'
#' @param setting String with value \code{"GLM"} or \code{"AFT"} indicating
#'                whether the approaches are fitted for a normally-distributed
#'                primary outcome \code{Y} (\code{"GLM"}) or a censored
#'                time-to-event primary outcome \code{Y} (\code{"AFT"}). Under
#'                the \code{"AFT"} setting, only \code{mult_reg} is
#'                available.
#' @param Y Numeric input vector of the primary outcome.
#' @param X Numeric input vector of the exposure variable.
#' @param K Numeric input vector of the intermediate outcome.
#' @param L Numeric input vector of the observed confounding factor.
#' @param C Numeric input vector of the censoring indicator under the AFT setting
#'          (must be coded 0 = censored, 1 = uncensored).
#'
#' @return Returns a list with point estimates of the parameters
#'         \code{point_estimates}, standard error estimates \code{SE_estimates}
#'         and p-values \code{pvalues}.
#'
#' @examples
#'
#' dat_GLM <- generate_data(setting = "GLM")
#' mult_reg(setting = "GLM", Y = dat_GLM$Y, X = dat_GLM$X, K = dat_GLM$K,
#'          L = dat_GLM$L)
#' res_reg(Y = dat_GLM$Y, X = dat_GLM$X, K = dat_GLM$K, L = dat_GLM$L)
#'
#' dat_AFT <- generate_data(setting = "AFT", a = 0.2, b = 4.75)
#' mult_reg(setting = "AFT", Y = dat_AFT$Y, X = dat_AFT$X, K = dat_AFT$K,
#'          L = dat_AFT$L, C = dat_AFT$C)
#'

#' @name traditional_regression_functions
NULL

#' @rdname traditional_regression_functions
#' @export

mult_reg <- function(setting = "GLM", Y = NULL, X = NULL, K = NULL,
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
    if (setting == "GLM") {
        fit_mult_reg <- stats::lm(Y ~ K + X + L)
        point_estimates <- c(summary(fit_mult_reg)$coefficients[1, 1],
                             summary(fit_mult_reg)$coefficients[2, 1],
                             summary(fit_mult_reg)$coefficients[3, 1],
                             summary(fit_mult_reg)$coefficients[4, 1])
        SE_estimates <- c(summary(fit_mult_reg)$coefficients[1, 2],
                          summary(fit_mult_reg)$coefficients[2, 2],
                          summary(fit_mult_reg)$coefficients[3, 2],
                          summary(fit_mult_reg)$coefficients[4, 2])
        pvalues <- c(summary(fit_mult_reg)$coefficients[1, 4],
                     summary(fit_mult_reg)$coefficients[2, 4],
                     summary(fit_mult_reg)$coefficients[3, 4],
                     summary(fit_mult_reg)$coefficients[4, 4])
        names(point_estimates) <- names(SE_estimates) <- names(pvalues) <-
          c("alpha_0", "alpha_1", "alpha_XY", "alpha_2")
    }
    if (setting == "AFT") {
        fit_mult_reg <- survival::survreg(survival::Surv(Y, C) ~ K + X + L,
                                          dist = "gaussian")
        point_estimates <- c(summary(fit_mult_reg)$table[1, 1],
                             summary(fit_mult_reg)$table[2, 1],
                             summary(fit_mult_reg)$table[3, 1],
                             summary(fit_mult_reg)$table[4, 1])
        SE_estimates <- c(summary(fit_mult_reg)$table[1, 2],
                          summary(fit_mult_reg)$table[2, 2],
                          summary(fit_mult_reg)$table[3, 2],
                          summary(fit_mult_reg)$table[4, 2])
        pvalues <- c(summary(fit_mult_reg)$table[1, 4],
                     summary(fit_mult_reg)$table[2, 4],
                     summary(fit_mult_reg)$table[3, 4],
                     summary(fit_mult_reg)$table[4, 4])
        names(point_estimates) <- names(SE_estimates) <- names(pvalues) <-
          c("alpha_0", "alpha_1", "alpha_XY", "alpha_2")
    }
    return(list(point_estimates = point_estimates,
                SE_estimates = SE_estimates, pvalues = pvalues))
}

#' @rdname traditional_regression_functions
#' @export

res_reg <- function(Y = NULL, X = NULL, K = NULL, L = NULL) {
    if (is.null(Y) | is.null(X) | is.null(K) | is.null(L)) {
        stop("Data of one or more variables is not supplied.")
    }
    data_help <- data.frame(Y = Y, X = X, K = K, L = L)
    data_help <- data_help[stats::complete.cases(data_help), ]
    fit_res_reg_1 <- stats::lm(data_help$Y ~ data_help$K + data_help$L)
    res <- fit_res_reg_1$residuals
    fit_res_reg_2 <- stats::lm(res ~ data_help$X)

    point_estimates <- c(summary(fit_res_reg_1)$coefficients[1, 1],
                         summary(fit_res_reg_1)$coefficients[2, 1],
                         summary(fit_res_reg_1)$coefficients[3, 1],
                         summary(fit_res_reg_2)$coefficients[1, 1],
                         summary(fit_res_reg_2)$coefficients[2, 1])
    SE_estimates <- c(summary(fit_res_reg_1)$coefficients[1, 2],
                      summary(fit_res_reg_1)$coefficients[2, 2],
                      summary(fit_res_reg_1)$coefficients[3, 2],
                      summary(fit_res_reg_2)$coefficients[1, 2],
                      summary(fit_res_reg_2)$coefficients[2, 2])
    pvalues <- c(summary(fit_res_reg_1)$coefficients[1, 4],
                 summary(fit_res_reg_1)$coefficients[2, 4],
                 summary(fit_res_reg_1)$coefficients[3, 4],
                 summary(fit_res_reg_2)$coefficients[1, 4],
                 summary(fit_res_reg_2)$coefficients[2, 4])
    names(point_estimates) <- names(SE_estimates) <- names(pvalues) <-
      c("alpha_0", "alpha_1", "alpha_2", "alpha_3", "alpha_XY")
    return(list(point_estimates = point_estimates,
                SE_estimates = SE_estimates, pvalues = pvalues))
}
