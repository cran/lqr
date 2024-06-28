\name{lqr}
\alias{lqr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust Linear Quantile Regression
}
\description{
It fits a robust linear quantile regression model using a new family of zero-quantile distributions for the error term. This family of distribution includes skewed versions of the Normal, Student's t, Laplace, Slash and Contaminated Normal distribution. It provides estimates and full inference. It also provides envelopes plots for assessing the fit and confidences bands when several quantiles are provided simultaneously.
}
\usage{
lqr(formula,data = NULL,subset = NULL,
               p=0.5,dist = "normal",
               nu=NULL,gamma=NULL,
               precision = 10^-6,envelope=FALSE,
               CI=0.95,silent = FALSE
)

#lqr(y~x, data, p = 0.5, dist = "normal")
#lqr(y~x, data, p = 0.5, dist = "t")
#lqr(y~x, data, p = 0.5, dist = "laplace")
#lqr(y~x, data, p = 0.5, dist = "slash")
#lqr(y~x, data, p = 0.5, dist = "cont")

#lqr(y~x, p = c(0.25,0.50,0.75), dist = "normal")
}

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
  \item{dist}{
represents the distribution to be used for the error term. The values are \code{normal} for Normal distribution, \code{t} for Student's t distribution, \code{laplace} for Laplace distribution, \code{slash} for Slash distribution and \code{cont} for the Contaminated normal distribution.
}
\item{nu}{
It represents the degrees of freedom when \code{dist = t}. For the Slash distribution (\code{dist = slash}) it is a shape parameter \eqn{\nu>0}. For the Contaminated Normal distribution, \eqn{\nu} is the parameter that represents the percentage of outliers. When is not provided, we use the MLE.
}
  \item{gamma}{
It represents a scale factor for the contaminated normal distribution. When is not provided, we use the MLE.
}
  \item{precision}{
The convergence maximum error permitted. By default is 10^-6.
}
  \item{envelope}{
if \code{TRUE}, it will show a confidence envelope for a curve based on bootstrap replicates.
By default it is \code{TRUE} when just one quantile is provided.
If a grid of quantiles is provided it will be \code{FALSE} by default.
}
  \item{CI}{
Confidence to be used for the Confidence Interval when a grid of quantiles is provided. Default = 0.95.
}
  \item{silent}{
if \code{FALSE} (by default), the function prints some output.
}
}
\details{
When a grid of quantiles is provided, a graphical summary with point estimates and Confidence Intervals for model parameters is shown.
}
\value{
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
\author{Christian E. Galarza <\email{cgalarza88@gmail.com}>, Luis Benites <\email{lsanchez@ime.usp.br}> and
Victor H. Lachos <\email{hlachos@ime.unicamp.br}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}
\note{
If a grid of quantiles is provided, the result will be a list of the same dimension where each element corresponds to each quantile as detailed above.
}

\examples{
\donttest{
#Example 1

##Load the data
data(ais)
attach(ais)

## Fitting a median regression with Normal errors (by default)

modelF = lqr(BMI~LBM,data = ais,subset = "(Sex==1)")
modelM = lqr(BMI~LBM,data = ais,subset = "(Sex==0)")

plot(LBM,BMI,col=Sex*2+1,
     xlab="Lean Body Mass",
     ylab="Body4 Mass Index",
     main="Quantile Regression")

abline(a = modelF$beta[1],b = modelF$beta[2],lwd=2,col=3)
abline(a = modelM$beta[1],b = modelM$beta[2],lwd=2,col=1)
legend(x = "topleft",legend = c("Male","Female"),lwd = 2,col = c(1,3))

#COMPARING SOME MODELS for median regression
modelN  = lqr(BMI~LBM,dist = "normal")
modelT  = lqr(BMI~LBM,dist = "t")
modelL  = lqr(BMI~LBM,dist = "laplace")

#Comparing AIC criteria
modelN$AIC;modelT$AIC;modelL$AIC

#This could be automatically done using best.lqr()
best.model = best.lqr(BMI~LBM,data = ais,
                      p = 0.75, #third quartile
                      criterion = "AIC")

#Let's use a grid of quantiles (no output)
modelfull = lqr(BMI~LBM,data = ais,
                p = seq(from = 0.10,to = 0.90,by = 0.05),
                dist = "normal",silent = TRUE)

#Plotting quantiles 0.10,0.25,0.50,0.75 and 0.90

if(TRUE){
  plot(LBM,BMI,xlab = "Lean Body Mass"
       ,ylab = "Body Mass Index", main = "Quantile Regression",pch=16)
  
  colvec = c(2,2,3,3,4)
  imodel = c(1,17,4,14,9)
  for(i in 1:5){
    abline(a = modelfull[[imodel[i]]]$beta[1],
           b = modelfull[[imodel[i]]]$beta[2],
           lwd=2,col=colvec[i])  
  }
  legend(x = "topleft",
         legend = rev(c("0.10","0.25","0.50","0.75","0.90")),
         lwd = 2,col = c(2,3,4,3,2))
}

#Example 2
##Load the data

data(crabs,package = "MASS")
attach(crabs)

## Fitting a median regression with Normal errors (by default) #Note the double quotes
crabsF = lqr(BD~FL,data = crabs,subset = "(sex=='F')")
crabsM = lqr(BD~FL,data = crabs,subset = "(sex=='M')")

if(TRUE){
  plot(FL,BD,col=as.numeric(sex)+1,
       xlab="Frontal lobe size",ylab="Body depth",main="Quantile Regression")
  abline(a = crabsF$beta[1],b = crabsF$beta[2],lwd=2,col=2)
  abline(a = crabsM$beta[1],b = crabsM$beta[2],lwd=2,col=3)
  legend(x = "topleft",legend = c("Male","Female"),
         lwd = 2,col = c(3,2))
}

#Median regression for different distributions

modelN  = lqr(BD~FL,dist = "normal")
modelT  = lqr(BD~FL,dist = "t")
modelL  = lqr(BD~FL,dist = "laplace")
modelS  = lqr(BD~FL,dist = "slash")
modelC  = lqr(BD~FL,dist = "cont" )

#Comparing AIC criterias
modelN$AIC;modelT$AIC;modelL$AIC;modelS$AIC;modelC$AIC

# best model based on BIC
best.lqr(BD~FL,criterion = "BIC")

#Let's use a grid of quantiles for the Student's t distribution
modelfull = lqr(BD~FL,data = crabs,
                p = seq(from = 0.10,to = 0.90,by = 0.05),
                dist = "t") # silent = FALSE

#Plotting quantiles 0.10,0.25,0.50,0.75 and 0.90
if(TRUE){
  plot(FL,BD,xlab = "Frontal lobe size"
       ,ylab = "Body depth", main = "Quantile Regression",pch=16)
  colvec = c(2,2,3,3,4)
  imodel = c(1,17,4,14,9)
  for(i in 1:5){
    abline(a = modelfull[[imodel[i]]]$beta[1],
           b = modelfull[[imodel[i]]]$beta[2],
           lwd=2,col=colvec[i])  
  }
  legend(x = "topleft",
         legend = rev(c("0.10","0.25","0.50","0.75","0.90")),
         lwd = 2,col = c(2,3,4,3,2))
}
}
}

\references{
Galarza, C., Lachos, V. H., Cabral, C. R. B., & Castro, C. L. (2017). Robust quantile regression using a generalized class of skewed distributions. Stat,
6(1), 113-130.

Wichitaksorn, N., Choy, S. T., & Gerlach, R. (2014). A generalized class of skew distributions and associated robust quantile regression models. Canadian Journal of Statistics, 42(4), 579-596.
}

\keyword{package}
\keyword{quantile}
\keyword{regression}
\keyword{skew}

\seealso{\code{\link{cens.lqr}},\code{\link{best.lqr}},\code{\link{Log.lqr}},
\code{\link{Log.best.lqr}},\code{\link{dSKD}}%,\code{\link[qrLMM]{QRLMM}},
%\code{\link[qrNLMM]{QRNLMM}}
}
