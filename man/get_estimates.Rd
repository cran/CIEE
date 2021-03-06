% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/point_estimates_ee.R
\name{get_estimates}
\alias{get_estimates}
\title{CIEE parameter point estimates}
\usage{
get_estimates(setting = "GLM", Y = NULL, X = NULL, K = NULL, L = NULL,
  C = NULL)
}
\arguments{
\item{setting}{String with value \code{"GLM"} or \code{"AFT"} indicating
whether CIEE point estimates are obtained for a
normally-distributed (\code{"GLM"}) or censored time-to-event
(\code{"AFT"}) primary outcome \code{Y}.}

\item{Y}{Numeric input vector for the primary outcome.}

\item{X}{Numeric input vector for the exposure variable.}

\item{K}{Numeric input vector for the intermediate outcome.}

\item{L}{Numeric input vector for the observed confounding factor.}

\item{C}{Numeric input vector for the censoring indicator under the AFT setting
(must be coded 0 = censored, 1 = uncensored).}
}
\value{
Returns a list with point estimates of the parameters. Under the
        AFT setting, the mean of the estimated true survival times is also
        computed and returned.
}
\description{
Function to perform CIEE to obtain point estimates under the GLM or AFT
setting for the analysis of a normally-distributed or censored time-to-event
primary outcome.
}
\details{
Under the GLM setting for the analysis of a normally-distributed primary
outcome Y, estimates of the parameters
\eqn{\alpha_0, \alpha_1, \alpha_2, \alpha_3, \sigma_1^2, \alpha_4, \alpha_{XY}, \sigma_2^2}{\alpha0, \alpha1, \alpha2, \alpha3, \sigma1^2, \alpha4, \alphaXY, \sigma2^2}
are obtained by constructing estimating equations for the models
\deqn{Y = \alpha_0 + \alpha_1 \cdot K + \alpha_2 \cdot X + \alpha_3 \cdot L + \epsilon_1, \epsilon_1 \sim N(0,\sigma_1^2)}{Y = \alpha0 + \alpha1*K + \alpha2*X + \alpha3*L + \epsilon1, \epsilon1 ~ N(0,\sigma1^2)}
\deqn{Y^* = Y - \overline{Y} - \alpha_1 \cdot (K-\overline{K})}{Y* = Y - mean(Y) - \alpha1*(K-mean(K))}
\deqn{Y^* = \alpha_0 + \alpha_{XY} \cdot X + \epsilon_2, \epsilon_2 \sim N(0,\sigma_2^2).}{Y* = \alpha0 + \alphaXY*X + \epsilon2, \epsilon2 ~ N(0,\sigma2^2).}
Under the AFT setting for the analysis of a censored time-to-event primary
outcome, estimates of the parameters
\eqn{\alpha_0, \alpha_1, \alpha_2, \alpha_3, \sigma_1, \alpha_4, \alpha_{XY}, \sigma_2^2}{\alpha0, \alpha1, \alpha2, \alpha3, \sigma1, \alpha4, \alphaXY, \sigma2^2}
are obtained by constructing
similar estimating equations based on a censored regression model and adding
an additional computation to estimate the true underlying survival times.
In addition to the parameter estimates, the mean of the estimated true
survival times is computed and returned in the output. For more details and
the underlying model, see the vignette.

For both settings, the point estimates based on estimating equations equal
least squares (and maximum likelihood) estimates, and are obtained using
the \code{\link[stats]{lm}} and \code{\link[survival]{survreg}}
functions for computational purposes.
}
\examples{

dat_GLM <- generate_data(setting = "GLM")
get_estimates(setting = "GLM", Y = dat_GLM$Y, X = dat_GLM$X, K = dat_GLM$K,
              L = dat_GLM$L)

dat_AFT <- generate_data(setting = "AFT", a = 0.2, b = 4.75)
get_estimates(setting = "AFT", Y = dat_AFT$Y, X = dat_AFT$X, K = dat_AFT$K,
              L = dat_AFT$L, C = dat_AFT$C)

}
