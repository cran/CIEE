#' Structural equation modeling approach
#'
#' Function which uses the \code{\link[lavaan]{sem}} function in the
#' \code{lavaan} package to fit the model
#' \deqn{L = \alpha_0 + \alpha_1 \cdot X + \epsilon_1, \epsilon_1 \sim N(0,\sigma_1^2)}{L = \alpha0 + \alpha1*X + \epsilon1, \epsilon1 ~ N(0,\sigma1^2)}
#' \deqn{K = \alpha_2 + \alpha_3 \cdot X + \alpha_4 \cdot L + \epsilon_2, \epsilon_2 \sim~ N(0,\sigma_2^2)}{K = \alpha2 + \alpha3*X + \alpha4*L + \epsilon2, \epsilon2 ~ N(0,\sigma2^2)}
#' \deqn{Y = \alpha_5 + \alpha_6 \cdot K + \alpha_{XY} \cdot X + \epsilon_3, \epsilon_3 \sim N(0,\sigma_3^2)}{Y = \alpha5 + \alpha6*K + \alphaXY*X + \epsilon3, \epsilon3 ~ N(0,\sigma3^2)}
#' in order to obtain point and standard error estimates
#' of the parameters
#' \eqn{\alpha_1, \alpha_3, \alpha_4, \alpha_6, \alpha_{XY}}{\alpha1, \alpha3, \alpha4, \alpha6, \alphaXY}
#' for the GLM setting.
#' See the vignette for more details.
#'
#' @param Y Numeric input vector for the primary outcome.
#' @param X Numeric input vector for the exposure variable.
#' @param K Numeric input vector for the intermediate outcome.
#' @param L Numeric input vector for the observed confounding factor.
#'
#' @return Returns a list with point estimates of the parameters
#'         (\code{point_estimates}), standard error estimates
#'         (\code{SE_estimates}) and p-values from large-sample
#'         Wald-type tests (\code{pvalues}).
#'
#' @examples
#'
#' dat <- generate_data(setting = "GLM")
#' sem_appl(Y = dat$Y, X = dat$X, K = dat$K, L = dat$L)
#'
#' @export
#'

sem_appl <- function(Y = NULL, X = NULL, K = NULL, L = NULL) {
    if (!requireNamespace("lavaan", quietly = TRUE)) {
        stop("Pkg needed for this function to work. Please install it.", call. = FALSE)
    }
    if (is.null(Y) | is.null(X) | is.null(K) | is.null(L)) {
        stop("Data of one or more variables are not supplied.")
    }
    data_help <- data.frame(Y = Y, X = X, K = K, L = L)
    data_help <- data_help[stats::complete.cases(data_help), ]
    model <- "
              L ~ X
              K ~ X + L
              Y ~ K + X
             "
    fit <- lavaan::sem(model, data = data_help)

    point_estimates <- c(fit@Fit@est[1], fit@Fit@est[2], fit@Fit@est[3],
                         fit@Fit@est[4], fit@Fit@est[5])
    SE_estimates <- c(fit@Fit@se[1], fit@Fit@se[2], fit@Fit@se[3],
                      fit@Fit@se[4], fit@Fit@se[5])
    pvalues <- 2 * stats::pnorm(-abs(point_estimates/SE_estimates))
    names(point_estimates) <- names(SE_estimates) <- names(pvalues) <-
      c("alpha_1", "alpha_3", "alpha_4", "alpha_6", "alpha_XY")
    return(list(point_estimates = point_estimates,
                SE_estimates = SE_estimates, pvalues = pvalues))
}
