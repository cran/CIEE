#' Data generation function
#'
#' Function to generate data with \code{n} observations of a primary
#' outcome \code{Y}, secondary outcome \code{K}, exposure \code{X}, and
#' measured as well as unmeasured confounders \code{L} and \code{U}, where
#' the primary outcome is a quantitative normally-distributed variable
#' (\code{setting} = \code{"GLM"}) or censored time-to-event outcome under
#' an accelerated failure time (AFT) model (\code{setting} = \code{"AFT"}).
#' Under the AFT setting, the observed time-to-event variable \code{T=exp(Y)}
#' as well as the censoring indicator \code{C} are also computed. \code{X}
#' is generated as a genetic exposure variable in the form of a single
#' nucleotide variant (SNV) in 0-1-2 additive coding with minor allele
#' frequency \code{maf}. \code{X} can be generated independently of \code{U}
#' (\code{X_orth_U} = \code{TRUE}) or dependent on \code{U}
#' (\code{X_orth_U} = \code{FALSE}). For more details regarding the underlying
#' model, see the vignette.
#'
#' @param setting String with value \code{"GLM"} or \code{"AFT"} indicating
#'                whether the primary outcome is generated as a
#'                normally-distributed quantitative outcome (\code{"GLM"}) or
#'                censored time-to-event outcome (\code{"AFT"}).
#' @param n Numeric. Sample size.
#' @param maf Numeric. Minor allele frequency of the genetic exposure variable.
#' @param cens Numeric. Desired percentage of censored individuals and has to be
#'             specified under the AFT setting. Note that the actual censoring
#'             rate is generated through specification of the parameters
#'             \code{a} and \code{b}, and \code{cens} is mostly used as a check
#'             whether the desired censoring rate is obtained through \code{a}
#'             and \code{b} (otherwise, a warning is issued).
#' @param a Integer for generating the desired censoring rate under the AFT
#'          setting. Has to be specified under the AFT setting.
#' @param b Integer for generating the desired censoring rate under the AFT
#'          setting. Has to be specified under the AFT setting.
#' @param aUL Numeric. Size of the effect of \code{U} on \code{L}.
#' @param aXL Numeric. Size of the effect of \code{X} on \code{L}.
#' @param aXK Numeric. Size of the effect of \code{X} on \code{K}.
#' @param aLK Numeric. Size of the effect of \code{L} on \code{K}.
#' @param aUY Numeric. Size of the effect of \code{U} on \code{Y}.
#' @param aKY Numeric. Size of the effect of \code{K} on \code{Y}.
#' @param aXY Numeric. Size of the effect of \code{X} on \code{Y}.
#' @param aLY Numeric. Size of the effect of \code{L} on \code{Y}.
#' @param mu_U Numeric. Expected value of \code{U}.
#' @param sd_U Numeric. Standard deviation of \code{U}.
#' @param X_orth_U Logical. Indicator whether \code{X} should be generated
#'                 independently of \code{U} (\code{X_orth_U} = \code{TRUE})
#'                 or dependent on \code{U} (\code{X_orth_U} = \code{FALSE}).
#' @param mu_X Numeric. Expected value of \code{X}.
#' @param sd_X Numeric. Standard deviation of \code{X}.
#' @param mu_L Numeric. Expected value of \code{L}.
#' @param sd_L Numeric. Standard deviation of \code{L}.
#' @param mu_K Numeric. Expected value of \code{K}.
#' @param sd_K Numeric. Standard deviation of \code{K}.
#' @param mu_Y Numeric. Expected value of \code{Y}.
#' @param sd_Y Numeric. Standard deviation of \code{Y}.
#'
#' @return A dataframe containing \code{n} observations of the variables \code{Y},
#'         \code{K}, \code{X}, \code{L}, \code{U}. Under the AFT setting,
#'         \code{T=exp(Y)} and the censoring indicator \code{C} (0 = censored,
#'         1 = uncensored) are also computed.
#'
#' @examples
#' # Generate data under the GLM setting with default values
#' dat_GLM <- generate_data()
#' head(dat_GLM)
#'
#' # Generate data under the AFT setting with default values
#' dat_AFT <- generate_data(setting = "AFT", a = 0.2, b = 4.75)
#' head(dat_AFT)
#'
#' @export
#'

generate_data <- function(setting = "GLM", n = 1000, maf = 0.2, cens = 0.3,
                          a = NULL, b = NULL, aXK = 0.2, aXY = 0.1, aXL = 0,
                          aKY = 0.3, aLK = 0, aLY = 0, aUY = 0, aUL = 0,
                          mu_X = NULL, sd_X = NULL, X_orth_U = TRUE, mu_U = 0,
                          sd_U = 1, mu_K = 0, sd_K = 1, mu_L = 0, sd_L = 1,
                          mu_Y = 0, sd_Y = 1) {
    U_out <- stats::rnorm(n, mean = mu_U, sd = sd_U)
    if (setting == "AFT" & (is.null(a) | is.null(b))) {
        stop("a and b have to be specified under the AFT setting.")
    }
    if (X_orth_U == TRUE) {
        X_out <- stats::rbinom(n, size = 2, prob = maf)
    }
    if (X_orth_U == FALSE) {
        X_out <- stats::pnorm(U_out, mean = mu_X, sd = sd_X)
        p <- 1 - maf
        for (j in 1:length(X_out)) {
            if (X_out[j] < p^2) {
                X_out[j] <- 0
                next
            }
            if (X_out[j] >= p^2 & X_out[j] < p^2 + 2 * p * (1 - p)) {
                X_out[j] <- 1
                next
            }
            if (X_out[j] >= p^2 + 2 * p * (1 - p)) {
                X_out[j] <- 2
                next
            }
        }
    }
    L_out <- aUL * U_out + aXL * X_out + stats::rnorm(n, mean = mu_L, sd = sd_L)
    K_out <- aXK * X_out + aLK * L_out + stats::rnorm(n, mean = mu_K, sd = sd_K)
    Y_out <- aUY * U_out + aKY * K_out + aXY * X_out + aLY * L_out +
             stats::rnorm(n, mean = mu_Y, sd = sd_Y)
    data <- data.frame(Y = Y_out, K = K_out, X = X_out, L = L_out, U = U_out)
    if (setting == "AFT") {
        T_help <- exp(Y_out)
        ### Create censoring indicator and censored times no censoring
        if (cens == 0) {
            T_out <- T_help
            C_out <- rep(1, n)  # C_out==0 is censored, C_out==1 is uncensored
        }
        if (!cens == 0) {
            # there is censoring; cens is the percentage of censored data
            T_cens <- stats::runif(n, min = a, max = b)  # a, b for desired censoring rate
            C_out <- as.numeric(T_help < T_cens)  # C==0 censored, C==1 uncensored
            T_out <- pmin(T_help, T_cens)
            Y_out <- log(T_out)
        }
        cens_out <- sum(abs(C_out - 1))/n
        data <- data.frame(Y = Y_out, K = K_out, X = X_out, L = L_out, U = U_out,
                           T = T_out, C = C_out)
        print(paste("The empirical censoring rate obtained through the specified parameters a=", a, " and b=", b, " is ", cens_out, ".", sep = ""))
        if (abs(cens_out - cens) > 0.1) {
            warning(paste("This obtained empirical censoring rate is quite different from the desired censoring rate cens=", cens, ". Please check and adapt values for a and b.",
                sep = ""))
        }
    }
    return(data)
}
