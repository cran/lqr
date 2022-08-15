\name{best.lqr}
\alias{best.lqr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Best Fit in Robust Linear Quantile Regression
}
\description{
It finds the best fit distribution in robust linear quantile regression model. It adjusts the Normal, Student's t, Laplace, Slash and Contaminated Normal models. It shows a summary table with the likelihood-based criterion, envelopes plots and the histogram of the residuals with fitted densities for all models. Estimates and full inference are provided for the best model.
}
\usage{
best.lqr(formula,data = NULL,subset = NULL,
              p = 0.5, precision = 10^-6,
              criterion = "AIC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
}
\item{data}{
an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}.
}
\item{subset}{
an optional string specifying a subset of observations to be used in the fitting process. Be aware of the use of double quotes in a proper way when necessary, e.g., in \code{"(sex=='F')"}.
}  
\item{p}{
An unique quantile or a set of quantiles related to the quantile regression.
}
  \item{precision}{
The convergence maximum error permitted. By default is 10^-6.
}
  \item{criterion}{
Likelihood-based criterion to be used for choosen the best model. It could be \code{AIC}, \code{BIC}, \code{HQ} or \code{loglik} (log-likelihood). By default \code{AIC} criterion will be used.
}
}
\details{
The \code{best.fit()} function finds the best model only for one quantile. For fitting a grid of quantiles \code{lqr()} might be used but the distribution must be provided.
}
\value{

For the best model:

\item{iter}{number of iterations.}
\item{criteria}{attained criteria value.}
\item{beta}{fixed effects estimates.}
\item{sigma}{scale parameter estimate for the error term.}
\item{nu}{Estimate of \code{nu} parameter detailed above.}
\item{gamma}{Estimate of \code{gamma} parameter detailed above.}
\item{SE}{Standard Error estimates.}
\item{table}{Table containing the inference for the fixed effects parameters.}
\item{loglik}{Log-likelihood value.}
\item{AIC}{Akaike information criterion.}
\item{BIC}{Bayesian information criterion.}
\item{HQ}{Hannan-Quinn information criterion.}
\item{fitted.values}{vector containing the fitted values.}
\item{residuals}{vector containing the residuals.}
}

\references{
Galarza, C., Lachos, V. H., Cabral, C. R. B., & Castro, C. L. (2017). Robust quantile regression using a generalized class of skewed distributions. Stat,
6(1), 113-130.

Wichitaksorn, N., Choy, S. T., & Gerlach, R. (2014). A generalized class of skew distributions and associated robust quantile regression models. Canadian Journal of Statistics, 42(4), 579-596.
}

\author{Christian E. Galarza <\email{cgalarza88@gmail.com}>, {Luis Benites} <\email{lsanchez@ime.usp.br}> and
Victor H. Lachos <\email{hlachos@ime.unicamp.br}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}

\keyword{package}
\keyword{quantile}
\keyword{regression}
\keyword{skew}

\seealso{\code{\link{lqr}},\code{\link{Log.lqr}},\code{\link{Log.best.lqr}},\code{\link{dSKD}}%,\code{\link[qrLMM]{QRLMM}},
%\code{\link[qrNLMM]{QRNLMM}}
}


\examples{
\donttest{

data(crabs,package = "MASS")

#Finding the best model for the 3rd quartile based on BIC
best.lqr(BD~FL,data = crabs, p = 0.75, criterion = "BIC")
 }
}
