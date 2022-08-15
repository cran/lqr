\name{Log.lqr}
\alias{Log.lqr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust Logistic Linear Quantile Regression
}
\description{
It performs the logistic transformation in Galarza et.al.(2020) (see references) for estimating quantiles for a bounded response. Once the response is transformed, it uses the \code{lqr} function.
}
\usage{
Log.lqr(formula,data = NULL,subset = NULL,
                   p=0.5,a=0,b=1,
                   dist = "normal",
                   nu=NULL,
                   gamma=NULL,
                   precision = 10^-6,
                   epsilon = 0.001,
                   CI=0.95,
                   silent = FALSE)
}

\arguments{
We will detail first the only three arguments that differ from \code{lqr} function.

  \item{a}{
lower bound for the response (default = 0)
}
 \item{b}{
upper bound for the response (default = 1)
}
\item{epsilon}{
a small quantity \eqn{\epsilon>0} that ensures that the logistic transform is defined for all values of the response.
}
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
  \item{CI}{
Confidence to be used for the Confidence Interval when a grid of quantiles is provided. Default = 0.95.
}
\item{silent}{
if \code{FALSE} (by default), the function prints some output.
}
}
\details{
We follow the transformation in Bottai et.al. (2009) defined as

\deqn{h(y)=logit(y)=log(\frac{y-a}{b-y})}

that implies

\deqn{Q_{y}(p)=\frac{b\,exp(X\beta) + a}{1 + exp(X\beta)}}

where \eqn{Q_{y}(p)} represents the conditional quantile of the response. Once estimates for the regression coefficients \eqn{\beta_p} are obtained, inference on \eqn{Q_{y}(p)} can then be made through the inverse transform above. This equation (as function) is provided in the output. See example.

The interpretation of the regression coefficients is analogous to the interpretation of the coefficients of a logistic regression for binary outcomes.

For example, let \eqn{x_1} be the gender (male = 0, female=1). Then \eqn{exp(\beta_{0.5,1})} represents the odds ratio of median score in males vs females, where the odds are defined using the score instead of a probability, \eqn{(y-a)/(b-y)}. When the covariate is continous, the respective \eqn{\beta} coeficient can be interpretated as the increment (or decrement) over the log(odd ratio) when the covariate increases one unit.


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
\author{Christian E. Galarza <\email{cgalarza88@gmail.com}>, {Luis Benites} <\email{lsanchez@ime.usp.br}> and
Victor H. Lachos <\email{hlachos@ime.unicamp.br}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}
\note{
When a grid of quantiles is provided, a graphical summary with point estimates and Confidence Intervals for model parameters is shown. Also, the result will be a list of the same dimension where each element corresponds to each quantile as detailed above.
}

\examples{
\donttest{

##Load the data
data(resistance)
attach(resistance)

#EXAMPLE 1.1

#Comparing the resistence to death of two types of tumor-cells.
#The response is a score in [0,4].

boxplot(score~type,ylab="score",xlab="type")

#Student't median logistic quantile regression
res = Log.lqr(score~type,data = resistance,a=0,b=4,dist="t")

# The odds ratio of median score in type B vs type A
exp(res$beta[2])

#Proving that exp(res$beta[2])  is approx median odd ratio
medA  = median(score[type=="A"])
medB  = median(score[type=="B"])
rateA = (medA - 0)/(4 - medA)
rateB = (medB - 0)/(4 - medB)
odd   = rateB/rateA

round(c(exp(res$beta[2]),odd),3)


#EXAMPLE 1.2
############

#Comparing the resistence to death depending of dose.

#descriptive
plot(dose,score,ylim=c(0,4),col="dark gray");abline(h=c(0,4),lty=2)
dosecat<-cut(dose, 6, ordered = TRUE)
boxplot(score~dosecat,ylim=c(0,4))
abline(h=c(0,4),lty=2)

#(Non logistic) Best quantile regression for quantiles
# 0.05, 0.50 and 0.95
p05 = best.lqr(score~poly(dose,3),data = resistance,p = 0.05)
p50 = best.lqr(score~poly(dose,3),data = resistance,p = 0.50)
p95 = best.lqr(score~poly(dose,3),data = resistance,p = 0.95)
res3  = list(p05,p50,p95)

plot(dose,score,ylim=c(-1,5),col="gray");abline(h=c(0,4),lty=2)
lines(sort(dose), p05$fitted.values[order(dose)], col='red', type='l')
lines(sort(dose), p50$fitted.values[order(dose)], col='blue', type='l')
lines(sort(dose), p95$fitted.values[order(dose)], col='red', type='l')

#Using Student's t logistic quantile regression for obtaining preditypeBions inside bounds

logp05 = Log.lqr(score~poly(dose,3),data = resistance,p = 0.05,b = 4,dist = "t") #a = 0 by default
logp50 = Log.lqr(score~poly(dose,3),data = resistance,p = 0.50,b = 4,dist = "t")
logp95 = Log.lqr(score~poly(dose,3),data = resistance,p = 0.95,b = 4,dist = "t")
res4  = list(logp05,logp50,logp95)

#No more predited curves out-of-bounds
plot(dose,score,ylim=c(-1,5),col="gray");abline(h=c(0,4),lty=2)
lines(sort(dose), logp05$fitted.values[order(dose)], col='red', type='l')
lines(sort(dose), logp50$fitted.values[order(dose)], col='blue', type='l')
lines(sort(dose), logp95$fitted.values[order(dose)], col='red', type='l')

#EXAMPLE 1.3
############

#A full model using dose and type for a grid of quantiles

res5 = Log.lqr(formula = score ~ poly(dose,3)*type,data = resistance,
               a = 0,b = 4,
               p = seq(from = 0.05,to = 0.95,by = 0.05),dist = "t",
               silent = TRUE)

#A nice plot

if(TRUE){
  par(mfrow=c(1,2))
  typeB = (resistance$type == "B")
  
  plot(dose,score,
       ylim=c(0,4),
       col=c(8*typeB + 1*!typeB),main="Type A")
  abline(h=c(0,4),lty=2)
  
  lines(sort(dose[!typeB]),
        res5[[2]]$fitted.values[!typeB][order(dose[!typeB])],
        col='red')
  
  lines(sort(dose[!typeB]),
        res5[[5]]$fitted.values[!typeB][order(dose[!typeB])],
        col='green')
  
  lines(sort(dose[!typeB]),
        res5[[10]]$fitted.values[!typeB][order(dose[!typeB])],
        col='blue',lwd=2)
  
  lines(sort(dose[!typeB]),
        res5[[15]]$fitted.values[!typeB][order(dose[!typeB])],
        col='green')
  
  lines(sort(dose[!typeB]),
        res5[[18]]$fitted.values[!typeB][order(dose[!typeB])],
        col='red')
  
  plot(dose,score,
       ylim=c(0,4),
       col=c(1*typeB + 8*!typeB),main="Type B")
  abline(h=c(0,4),lty=2)
  
  lines(sort(dose[typeB]),
        res5[[2]]$fitted.values[typeB][order(dose[typeB])],
        col='red')
  
  lines(sort(dose[typeB]),
        res5[[5]]$fitted.values[typeB][order(dose[typeB])],
        col='green')
  
  lines(sort(dose[typeB]),
        res5[[10]]$fitted.values[typeB][order(dose[typeB])],
        col='blue',lwd=2)
  
  lines(sort(dose[typeB]),
        res5[[15]]$fitted.values[typeB][order(dose[typeB])],
        col='green')
  
  lines(sort(dose[typeB]),
        res5[[18]]$fitted.values[typeB][order(dose[typeB])],
        col='red')
}
}
}

\references{
Galarza, C.M., Zhang P. and Lachos, V.H. (2020). Logistic Quantile Regression for Bounded Outcomes Using
a Family of Heavy-Tailed Distributions. Sankhya B: The Indian Journal of Statistics. \doi{10.1007/s13571-020-00231-0}

Galarza, C., Lachos, V. H., Cabral, C. R. B., & Castro, C. L. (2017). Robust quantile regression using a generalized class of skewed distributions. Stat, 6(1), 113-130.
}

\keyword{bounded}
\keyword{logistic}
\keyword{quantile}
\keyword{regression}
\keyword{skew}

\seealso{\code{\link{Log.best.lqr}},\code{\link{best.lqr}},\code{\link{dSKD}}%,\code{\link[qrLMM]{QRLMM}},
%\code{\link[qrNLMM]{QRNLMM}}
}
